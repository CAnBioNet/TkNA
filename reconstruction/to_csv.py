#!/usr/bin/env python
import argparse
import itertools
import json
import numpy
from pathlib import Path
import xarray

from util import CsvWriter, readDataZip, parseConfigFile
from util.configs import aggregateConfigSpec, singleCellConfigSpec

numpy.seterr(all="raise")

def getArgs():
	parser = argparse.ArgumentParser(description="Write CSVs in the specified formats from data.", add_help=False)
	requiredArgGroup = parser.add_argument_group("required arguments")
	optionalArgGroup = parser.add_argument_group("optional arguments")
	requiredArgGroup.add_argument("--data-file", type=str, dest="dataFile", required=True, help="ZIP file containing data")
	requiredArgGroup.add_argument("--config-file", type=str, dest="configFile", required=True, help="Config file used to generate the data")
	requiredArgGroup.add_argument("--out-dir", type=str, dest="outDir", required=True, help="Directory to write the CSV files to")
	optionalArgGroup.add_argument("--singlecell", "-s", action="store_true", default=False, help="Work with single-cell rather than aggregate data")
	optionalArgGroup.add_argument("-h", "--help", action="help", help="Show this help message and exit")
	args = parser.parse_args()

	args.dataFile = Path(args.dataFile)
	if not args.dataFile.exists():
		raise Exception("Specified data file does not exist.")

	args.configFile = Path(args.configFile)
	if not args.configFile.exists():
		raise Exception("Specified config file does not exist.")

	args.outDir = Path(args.outDir)
	if args.outDir.exists() and not args.outDir.is_dir():
		raise Exception("Specified output directory exists but is not a directory.")

	return args

def computeFoldChanges_(experimentData, foldChangeType, treatments, *args):
	grouped = experimentData.groupby("treatment")
	treatmentData = [grouped[treatment] for treatment in treatments]
	combinedOrganismData = [getattr(treatmentDataset, foldChangeType)("organism") for treatmentDataset in treatmentData]
	return numpy.log2(combinedOrganismData[0] / combinedOrganismData[1])

def setupMeasurableCsv(data, config):
	foldChangeType = config["foldChangeType"]

	csvConfig = CsvWriter.Config("measurable",
		CsvWriter.Coordinate("ID"),
		CsvWriter.Property("Type", "originalData", "measurableType"),
		CsvWriter.Per("Median Value ({})", "medianValue", "experiment"),
		CsvWriter.Per("Mean Value ({})", "meanValue", "experiment"),
		CsvWriter.Per("Comparison p-value ({})", "differencePValues", "experiment"),
		CsvWriter.Per("Median Log2 Fold Change ({})", "medianFoldChanges", "experiment"),
		CsvWriter.Per("Mean Log2 Fold Change ({})", "meanFoldChanges", "experiment"),
		CsvWriter.Per(foldChangeType.capitalize() + " Log2 Fold Change Direction ({})", "foldChangeSigns", "experiment"),
		CsvWriter.Column(foldChangeType.capitalize() + " Log2 Fold Changes Consistent", "consistentFoldChange"),
		CsvWriter.Column("Combined Comparison p-value", "combinedDifferencePValues"),
		CsvWriter.Column("Corrected Comparison p-value", "correctedDifferencePValues")
	)

	data["meanValue"] = None
	data["medianValue"] = None
	data["meanFoldChanges"] = None
	data["medianFoldChanges"] = None
	data["consistentFoldChange"] = None

	dataKeys = csvConfig.getDataKeys()
	dataKeys.append("foldChangeSigns")
	if not all(key in data for key in dataKeys):
		return

	data["meanValue"] = data["originalData"].groupby("experiment").map(lambda a: a.mean(dim="organism"))
	data["medianValue"] = data["originalData"].groupby("experiment").map(lambda a: a.median(dim="organism"))
	data["consistentFoldChange"] = xarray.apply_ufunc(lambda signs: numpy.all(signs == signs[0]), data["foldChangeSigns"], input_core_dims=[["experiment"]], vectorize=True)

	treatments = config["comparisonTreatments"]
	data["meanFoldChanges"] = data["originalData"].groupby("experiment").map(computeFoldChanges_, args=("mean",treatments))
	data["medianFoldChanges"] = data["originalData"].groupby("experiment").map(computeFoldChanges_, args=("median",treatments))

	return csvConfig, data

def writeComparisons(data, config, csvConfig, outDir):
	fileName = "all_comparisons.csv"
	CsvWriter.writeCsv(outDir / fileName, csvConfig, data, data["differencePValues"].coords["measurable"].data)

def writeNodes(data, config, csvConfig, outDir):
	if "filteredData" not in data:
		return

	fileName = "node_comparisons.csv"
	CsvWriter.writeCsv(outDir / fileName, csvConfig, data, data["filteredData"].coords["measurable"].data)

def setupEdgeCsv(data, config):
	correlationMethod = config["correlationMethod"]
	correlationFilterMethod = config["correlationFilterMethod"]
	if correlationFilterMethod == "percentagreement":
		percentAgreementThreshold = config["correlationFilterPercentAgreementThreshold"]
		consistencyDescriptor = "{}% Agreement".format(percentAgreementThreshold * 100)
	else:
		consistencyDescriptor = "All Agree"

	csvConfig = CsvWriter.Config(["measurable1", "measurable2"],
		CsvWriter.CoordinateFormatted("Edge name", "{}<==>{}"),
		CsvWriter.CoordComponent("partner1", 0),
		CsvWriter.CoordComponent("partner2", 1),
		CsvWriter.CoordComponentPropertyFormatted("Edge Type", "originalData", "{}<==>{}", "measurable", "measurableType"),
		CsvWriter.Per(correlationMethod.capitalize() + " corr pval ({})", "correlationPValues", "metatreatment"),
		CsvWriter.Column("Combined corr Pvalue", "combinedCorrelationPValues"),
		CsvWriter.Per("rho Coefficient ({})", "correlationCoefficients", "metatreatment"),
		CsvWriter.Column("Mean rho Coefficient", "combinedCoefficients"),
		CsvWriter.Column("Combined corr FDR", "correctedCorrelationPValues"),
		CsvWriter.CoordComponent("partner1 name", 0),
		CsvWriter.CoordComponentPer("partner1_MedianValue ({})", "medianValue", 0, "measurable", "experiment"),
		CsvWriter.CoordComponentPer("partner1_MeanValue ({})", "meanValue", 0, "measurable", "experiment"),
		CsvWriter.CoordComponentPer("partner1_MedianLog2FoldChange ({})", "medianFoldChanges", 0, "measurable", "experiment"),
		CsvWriter.CoordComponentPer("partner1_MeanLog2FoldChange ({})", "meanFoldChanges", 0, "measurable", "experiment"),
		CsvWriter.CoordComponent("partner2 name", 1),
		CsvWriter.CoordComponentPer("partner2_MedianValue ({})", "medianValue", 1, "measurable", "experiment"),
		CsvWriter.CoordComponentPer("partner2_MeanValue ({})", "meanValue", 1, "measurable", "experiment"),
		CsvWriter.CoordComponentPer("partner2_MedianLog2FoldChange ({})", "medianFoldChanges", 1, "measurable", "experiment"),
		CsvWriter.CoordComponentPer("partner2_MeanLog2FoldChange ({})", "meanFoldChanges", 1, "measurable", "experiment"),
		CsvWriter.Column("Correlations Passed Consistency Filter ({})".format(consistencyDescriptor), "correlationFilter"),
		CsvWriter.Column("All Non-PUC Filters Passed", "nonPucPassed"),
		CsvWriter.Column("combined Coefficient correlation Direction", "combinedCorrelationSigns"),
		CsvWriter.CoordComponentColumn("partner1_FC_direction", "combinedFoldChangeSigns", 0, "measurable"),
		CsvWriter.CoordComponentColumn("partner2_FC_direction", "combinedFoldChangeSigns", 1, "measurable"),
		CsvWriter.Column("IfFoldChangeDirectionMatch", "foldChangeSignProducts"),
		CsvWriter.Column("PUC", "expectedEdgeFilterInt"),
		CsvWriter.Column("Final Network Value (0: No edge, 1: Positive edge, -1: Negative edge)", "edges")
	)

	data["combinedCoefficients"] = None
	data["expectedEdgeFilterInt"] = None
	data["nonPucPassed"] = None

	data["meanValue"] = None
	data["medianValue"] = None
	data["meanFoldChanges"] = None
	data["medianFoldChanges"] = None

	dataKeys = csvConfig.getDataKeys()
	dataKeys.append("originalData")
	dataKeys.append("expectedEdgeFilter")
	if not all(key in data for key in dataKeys):
		return

	data["combinedCoefficients"] = data["correlationCoefficients"].mean(dim="metatreatment")
	data["expectedEdgeFilterInt"] = data["expectedEdgeFilter"].astype(int)
	data["nonPucPassed"] = data["diagonalFilter"] & data["individualCorrelationPValueFilter"] & data["combinedCorrelationPValueFilter"] & data["correctedCorrelationPValueFilter"] & data["correlationFilter"]

	data["meanValue"] = data["originalData"].groupby("experiment").map(lambda a: a.mean(dim="organism"))
	data["medianValue"] = data["originalData"].groupby("experiment").map(lambda a: a.median(dim="organism"))

	treatments = config["comparisonTreatments"]
	data["meanFoldChanges"] = data["originalData"].groupby("experiment").map(computeFoldChanges_, args=("mean",treatments))
	data["medianFoldChanges"] = data["originalData"].groupby("experiment").map(computeFoldChanges_, args=("median",treatments))

	return csvConfig, data

def writeCorrelations(data, config, csvConfig, outDir):
	fileName = "correlations_bw_signif_measurables.csv"
	CsvWriter.writeCsv(outDir / fileName, csvConfig, data, list(itertools.combinations(data["filteredData"].coords["measurable"].data, 2)))

def writeSummary(data, config, csvConfig, outDir):
	includedEdgeIndices = numpy.argwhere((data["edges"] != 0 & ~(data["edges"].isnull())).data)
	includedEdgeEntries = [data["edges"][index[0], index[1]] for index in includedEdgeIndices]
	includedEdges = {frozenset((entry.measurable1.item(), entry.measurable2.item())) for entry in includedEdgeEntries}
	includedEdges = [tuple(edge) for edge in includedEdges]

	fileName = "network_output_comp.csv"
	CsvWriter.writeCsv(outDir / fileName, csvConfig, data, includedEdges)

def writeMeasurableCsvSingleCell(data, config, filePath, nodesOnly):
	combinedDifferencePValuesStacked = data["combinedDifferencePValues"].stack({"measurableAndCellType": ("measurable", "differentialCellType")})
	# Reset index to match the format of the rest of the data
	combinedDifferencePValuesNoNans = combinedDifferencePValuesStacked.dropna("measurableAndCellType").reset_index("measurableAndCellType")
	data["combinedDifferencePValuesStacked"] = combinedDifferencePValuesStacked.reset_index("measurableAndCellType")
	data["correctedDifferencePValuesStacked"] = data["correctedDifferencePValues"].stack({"measurableAndCellType": ("measurable", "differentialCellType")}).reset_index("measurableAndCellType")

	if nodesOnly:
		combinedDifferencePValues = combinedDifferencePValuesStacked.reset_index("measurableAndCellType")
		cellTypeCoords = [data["filteredData"].sel(measurableAndCellType=m).cellType.item() for m in data["filteredData"].coords["measurableAndCellType"].data]
		measurableCoords = [data["filteredData"].sel(measurableAndCellType=m).measurable.item() for m in data["filteredData"].coords["measurableAndCellType"].data]

		coordArr = data["filteredData"]
		indexList = [numpy.argwhere(((combinedDifferencePValues.differentialCellType==cellType) & (combinedDifferencePValues.measurable==measurable)).data).item() for cellType, measurable in zip(cellTypeCoords, measurableCoords)]
		cellTypeDim = "cellType"
	else:
		coordArr = combinedDifferencePValuesNoNans
		indexList = numpy.argwhere(~numpy.isnan(combinedDifferencePValuesStacked.data))
		cellTypeDim = "differentialCellType"

	data["foldChangesStacked"] = data["foldChanges"].stack({"measurableAndCellType": ("measurable", "differentialCellType")}).reset_index("measurableAndCellType")
	data["foldChangeSignsStacked"] = data["foldChangeSigns"].stack({"measurableAndCellType": ("measurable", "differentialCellType")}).reset_index("measurableAndCellType")
	data["foldChangeFilterStacked"] = data["foldChangeFilter"].stack({"measurableAndCellType": ("measurable", "differentialCellType")}).reset_index("measurableAndCellType")

	csvConfig = CsvWriter.Config("measurableAndCellType",
		CsvWriter.CoordinateFunction("Measurable", lambda m: "{}".format(coordArr.coords["measurable"][m].item())),
		CsvWriter.CoordinateFunction("Cell Type", lambda m: "{}".format(coordArr.coords[cellTypeDim][m].item())),
		CsvWriter.CoordinateFunction("Measurable Type", lambda m: "{}".format(coordArr.coords["measurableType"][m].item())),
		CsvWriter.Per("Average Value ({})", "stacked", "organism", coordMap=indexList),
		CsvWriter.Per("Average Log2 Fold Change ({})", "foldChangesStacked", "differentialExperiment", coordMap=indexList),
		CsvWriter.Per("Average Log2 Fold Change Direction ({})", "foldChangeSignsStacked", "differentialExperiment", coordMap=indexList),
		CsvWriter.Column("Average Log2 Fold Changes Consistent", "foldChangeFilterStacked", coordMap=indexList),
		CsvWriter.Per("Corrected Comparison p-value ({})", "correctedDifferencePValuesStacked", "differentialExperiment", coordMap=indexList),
		CsvWriter.Column("Combined Comparison p-value", "combinedDifferencePValuesStacked", coordMap=indexList)
	)

	dataKeys = csvConfig.getDataKeys()
	dataKeys.append("combinedDifferencePValues")
	dataKeys.append("foldChanges")
	dataKeys.append("foldChangeSigns")
	dataKeys.append("foldChangeFilter")
	if not all(key in data for key in dataKeys):
		return

	CsvWriter.writeCsv(filePath, csvConfig, data, coordArr.coords["measurableAndCellType"].data)

def writeComparisonsSingleCell(data, config, outDir):
	writeMeasurableCsvSingleCell(data, config, outDir / "all_comparisons.csv", False)

def writeNodesSingleCell(data, config, outDir):
	fileName = "node_comparisons.csv"
	writeMeasurableCsvSingleCell(data, config, outDir / fileName, True)

def writeEdgeCsvSingleCell(data, config, filePath, finalOnly=False):
	if finalOnly:
		allEdges = list(itertools.combinations(data["filteredData"].coords["measurableAndCellType"].data, 2))
		edgeList = [edge for edge in allEdges if data["edges"].sel(measurableAndCellType1=edge[0], measurableAndCellType2=edge[1]) != 0]
	else:
		edgeList = list(itertools.combinations(data["filteredData"].coords["measurableAndCellType"].data, 2))

	cellTypeCoords = [data["filteredData"].sel(measurableAndCellType=m).cellType.item() for m in data["filteredData"].coords["measurableAndCellType"].data]
	measurableCoords = [data["filteredData"].sel(measurableAndCellType=m).measurable.item() for m in data["filteredData"].coords["measurableAndCellType"].data]
	indexList = [numpy.argwhere(((data["stacked"].cellType==cellType) & (data["stacked"].measurable==measurable)).data).item() for cellType, measurable in zip(cellTypeCoords, measurableCoords)]

	data["foldChangesStacked"] = data["foldChanges"].stack({"measurableAndCellType": ("measurable", "differentialCellType")}).reset_index("measurableAndCellType")
	data["combinedFoldChangeSignsStacked"] = data["combinedFoldChangeSigns"].stack({"measurableAndCellType": ("measurable", "differentialCellType")}).reset_index("measurableAndCellType")

	csvConfig = CsvWriter.Config(["measurableAndCellType1", "measurableAndCellType2"],
		CsvWriter.Property("Measurable 1", "correlationCoefficients", "measurable1"),
		CsvWriter.Property("Measurable 2", "correlationCoefficients", "measurable2"),
		CsvWriter.PropertiesFormatted("Measurables", "correlationCoefficients", "{}<==>{}", ["measurable1", "measurable2"]),
		CsvWriter.Property("Cell Type 1", "correlationCoefficients", "cellType1"),
		CsvWriter.Property("Cell Type 2", "correlationCoefficients", "cellType2"),
		CsvWriter.PropertiesFormatted("Cell Types", "correlationCoefficients", "{}<==>{}", ["cellType1", "cellType2"]),
		CsvWriter.CoordComponentPer("Partner 1 Average Value ({})", "stacked", 0, "measurableAndCellType", "organism", coordMap=indexList),
		CsvWriter.CoordComponentColumn("Partner 1 Average Log2 Fold Change", "foldChangesStacked", 0, "measurableAndCellType", coordMap=indexList),
		CsvWriter.CoordComponentColumn("Partner 1 Fold Change Direction", "combinedFoldChangeSignsStacked", 0, "measurableAndCellType", coordMap=indexList),
		CsvWriter.CoordComponentPer("Partner 2 Average Value ({})", "stacked", 1, "measurableAndCellType", "organism", coordMap=indexList),
		CsvWriter.CoordComponentColumn("Partner 2 Average Log2 Fold Change", "foldChangesStacked", 1, "measurableAndCellType", coordMap=indexList),
		CsvWriter.CoordComponentColumn("Partner 2 Fold Change Direction", "combinedFoldChangeSignsStacked", 1, "measurableAndCellType", coordMap=indexList),
		CsvWriter.Column("Fold Change Signs Match", "foldChangeSignProducts"),
		CsvWriter.Per("Correlation p-value ({})", "correlationPValues", "metatreatment"),
		CsvWriter.Per("Correlation Coefficient ({})", "correlationCoefficients", "metatreatment"),
		CsvWriter.Column("Correlation Coefficients Consistent", "correlationFilter"),
		CsvWriter.Column("Combined Correlation p-value", "combinedCorrelationPValues"),
		CsvWriter.Column("Corrected Correlation p-value", "correctedCorrelationPValues"),
		CsvWriter.Column("PUC", "expectedEdgeFilter"),
		CsvWriter.Column("Final Network Value (0: No edge, 1: Positive edge, -1: Negative edge)", "edges")
	)

	dataKeys = csvConfig.getDataKeys()
	dataKeys.append("foldChanges")
	dataKeys.append("combinedFoldChangeSigns")
	if not all(key in data for key in dataKeys):
		return

	CsvWriter.writeCsv(filePath, csvConfig, data, edgeList)

def writeCorrelationsSingleCell(data, config, outDir):
	fileName = "correlations_bw_signif_measurables.csv"
	writeEdgeCsvSingleCell(data, config, outDir / fileName, False)

def writeSummarySingleCell(data, config, outDir):
	fileName = "network_output_comp.csv"
	writeEdgeCsvSingleCell(data, config, outDir / fileName, True)

def writeConfigValues(config, outDir):
	configValuesFilePath = outDir / "config_values.txt"
	configValuesFile = open(configValuesFilePath, "w")

	def writeThresholds(thresholdKey):
		pValueThresholds = config[thresholdKey]
		arePerTypeThresholds = False
		types = []
		for threshold in pValueThresholds.values():
			if isinstance(threshold, dict):
				arePerTypeThresholds = True
				types = list(threshold.keys())
				break
		if arePerTypeThresholds:
			typeDicts = {}
			for thresholdType, thresholdOrDict in pValueThresholds.items():
				if isinstance(thresholdOrDict, dict):
					for type_, threshold in thresholdOrDict.items():
						if type_ not in typeDicts:
							typeDicts[type_] = {}
						typeDicts[type_][thresholdType] = threshold
				else:
					for type_ in types:
						typeDicts[type_][thresholdType] = thresholdOrDict
			for type_, typeDict in typeDicts.items():
				configValuesFile.write("\t{}: {}\n".format(type_, ", ".join("{}: {}".format(thresholdType, threshold) for thresholdType, threshold in typeDict.items())))
		else:
			for thresholdType, threshold in pValueThresholds.items():
				configValuesFile.write("\t{}: {}\n".format(thresholdType, threshold))

	configValuesFile.write("Comparison thresholds:\n")
	writeThresholds("differencePValueThresholds")

	configValuesFile.write("\n")

	configValuesFile.write("Correlation thresholds:\n")
	writeThresholds("correlationPValueThresholds")

	configValuesFile.close()

if __name__ == "__main__":
	args = getArgs()
	data = readDataZip(args.dataFile)

	configSpec = singleCellConfigSpec if args.singlecell else aggregateConfigSpec
	config = parseConfigFile(configSpec, args.configFile)

	args.outDir.mkdir(parents=True, exist_ok=True)

	if args.singlecell:
		writeComparisonsSingleCell(data, config, args.outDir)
		writeNodesSingleCell(data, config, args.outDir)
		writeCorrelationsSingleCell(data, config, args.outDir)
		writeSummarySingleCell(data, config, args.outDir)
	else:
		nodeCsvConfig, data = setupMeasurableCsv(data, config)
		writeComparisons(data, config, nodeCsvConfig, args.outDir)
		writeNodes(data, config, nodeCsvConfig, args.outDir)

		edgeCsvConfig, data = setupEdgeCsv(data, config)
		writeCorrelations(data, config, edgeCsvConfig, args.outDir)
		writeSummary(data, config, edgeCsvConfig, args.outDir)

	writeConfigValues(config, args.outDir)

