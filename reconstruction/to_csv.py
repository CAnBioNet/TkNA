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

class MissingDataError(Exception):
	def __init__(self, keys):
		self.keys = keys
		self.message = f"Could not find tables matching keys: {keys}"
		super().__init__(self.message)

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

# Used to conditionally include an item in a list, e.g. for columns conditionally included in a CSV table
# Works by returning a collection that is unpacked, allowing for the possibility of
# adding nothing to the containing list if the conditional is false
# Usage: list = [..., *listItemIf(<expression>, <conditional expression>), ...]
def listItemIf(item, conditional):
	if conditional:
		return (item,)
	else:
		return ()

def setupMeasurableCsv(data, config):
	foldChangeType = config["foldChangeType"]

	csvConfig = CsvWriter.Config("measurable",
		CsvWriter.Coordinate("ID"),
		CsvWriter.Property("Type", "originalData", "measurableType"),
		CsvWriter.Per("Median Value (Experiment: {})", "medianValueE", "experiment"),
		CsvWriter.Per("Median Value (Treatment: {})", "medianValueT", "treatment"),
		CsvWriter.Per("Median Value (Experiment: {}, Treatment: {})", "medianValueET", ["experiment", "treatment"]),
		CsvWriter.Per("Mean Value (Experiment: {})", "meanValueE", "experiment"),
		CsvWriter.Per("Mean Value (Treatment: {})", "meanValueT", "treatment"),
		CsvWriter.Per("Mean Value (Experiment: {}, Treatment: {})", "meanValueET", ["experiment", "treatment"]),
		CsvWriter.Per("Comparison p-value ({})", "differencePValues", "experiment"),
		CsvWriter.Per("Median Log2 Fold Change ({})", "medianFoldChanges", "experiment"),
		CsvWriter.Per("Mean Log2 Fold Change ({})", "meanFoldChanges", "experiment"),
		CsvWriter.Per(foldChangeType.capitalize() + " Log2 Fold Change Direction ({})", "foldChangeSigns", "experiment"),
		CsvWriter.Column(foldChangeType.capitalize() + " Log2 Fold Changes Consistent", "consistentFoldChange"),
		CsvWriter.Column("Combined Comparison p-value", "combinedDifferencePValues"),
		CsvWriter.Column("Corrected Comparison p-value", "correctedDifferencePValues")
	)

	# Initialize computed columns so that they are counted as present in the dataset
	data["medianValueE"] = None
	data["medianValueT"] = None
	data["medianValueET"] = None
	data["meanValueE"] = None
	data["meanValueT"] = None
	data["meanValueET"] = None

	data["consistentFoldChange"] = None

	dataKeys = csvConfig.getDataKeys()
	# Add keys for the data the computed columns depend on to ensure its existence
	dataKeys.append("foldChangeSigns")

	# Check that all of the required data is present in the dataset
	missingKeys = [key for key in dataKeys if key not in data]
	if len(missingKeys) > 0:
		raise MissingDataError(missingKeys)

	# Determine computed columns
	data["medianValueE"] = data["originalData"].groupby("experiment").map(lambda a: a.median(dim="organism"))
	data["medianValueT"] = data["originalData"].groupby("treatment").map(lambda a: a.median(dim="organism"))
	data["medianValueET"] = data["originalData"].groupby("experiment").map(lambda a: a.groupby("treatment").map(lambda a: a.median(dim="organism")))

	data["meanValueE"] = data["originalData"].groupby("experiment").map(lambda a: a.mean(dim="organism"))
	data["meanValueT"] = data["originalData"].groupby("treatment").map(lambda a: a.mean(dim="organism"))
	data["meanValueET"] = data["originalData"].groupby("experiment").map(lambda a: a.groupby("treatment").map(lambda a: a.mean(dim="organism")))

	data["consistentFoldChange"] = xarray.apply_ufunc(lambda signs: numpy.all(signs == signs[0]), data["foldChangeSigns"], input_core_dims=[["experiment"]], vectorize=True)

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
		CsvWriter.CoordComponentPer("partner1_MedianValue (Experiment: {})", "medianValueE", 0, "measurable", "experiment"),
		CsvWriter.CoordComponentPer("partner1_MedianValue (Treatment: {})", "medianValueT", 0, "measurable", "treatment"),
		CsvWriter.CoordComponentPer("partner1_MedianValue (Experiment: {}, Treatment: {})", "medianValueET", 0, "measurable", ["experiment", "treatment"]),
		CsvWriter.CoordComponentPer("partner1_MeanValue (Experiment: {})", "meanValueE", 0, "measurable", "experiment"),
		CsvWriter.CoordComponentPer("partner1_MeanValue (Treatment: {})", "meanValueT", 0, "measurable", "treatment"),
		CsvWriter.CoordComponentPer("partner1_MeanValue (Experiment: {}, Treatment: {})", "meanValueET", 0, "measurable", ["experiment", "treatment"]),
		CsvWriter.CoordComponentPer("partner1_MedianLog2FoldChange ({})", "medianFoldChanges", 0, "measurable", "experiment"),
		CsvWriter.CoordComponentPer("partner1_MeanLog2FoldChange ({})", "meanFoldChanges", 0, "measurable", "experiment"),
		CsvWriter.CoordComponent("partner2 name", 1),
		CsvWriter.CoordComponentPer("partner2_MedianValue (Experiment: {})", "medianValueE", 1, "measurable", "experiment"),
		CsvWriter.CoordComponentPer("partner2_MedianValue (Treatment: {})", "medianValueT", 1, "measurable", "treatment"),
		CsvWriter.CoordComponentPer("partner2_MedianValue (Experiment: {}, Treatment: {})", "medianValueET", 1, "measurable", ["experiment", "treatment"]),
		CsvWriter.CoordComponentPer("partner2_MeanValue (Experiment: {})", "meanValueE", 1, "measurable", "experiment"),
		CsvWriter.CoordComponentPer("partner2_MeanValue (Treatment: {})", "meanValueT", 1, "measurable", "treatment"),
		CsvWriter.CoordComponentPer("partner2_MeanValue (Experiment: {}, Treatment: {})", "meanValueET", 1, "measurable", ["experiment", "treatment"]),
		CsvWriter.CoordComponentPer("partner2_MedianLog2FoldChange ({})", "medianFoldChanges", 1, "measurable", "experiment"),
		CsvWriter.CoordComponentPer("partner2_MeanLog2FoldChange ({})", "meanFoldChanges", 1, "measurable", "experiment"),
		CsvWriter.Column("Correlations Passed Consistency Filter ({})".format(consistencyDescriptor), "correlationFilter"),
		*listItemIf(CsvWriter.Column("All Non-PUC Filters Passed", "nonPucPassed"), not config["noPUC"]),
		CsvWriter.Column("combined Coefficient correlation Direction", "combinedCorrelationSigns"),
		CsvWriter.CoordComponentColumn("partner1_FC_direction", "combinedFoldChangeSigns", 0, "measurable"),
		CsvWriter.CoordComponentColumn("partner2_FC_direction", "combinedFoldChangeSigns", 1, "measurable"),
		*listItemIf(CsvWriter.Column("IfFoldChangeDirectionMatch", "foldChangeSignProducts"), not config["noPUC"]),
		*listItemIf(CsvWriter.Column("PUC", "expectedEdgeFilterInt"), not config["noPUC"]),
		CsvWriter.Column("Final Network Value (0: No edge, 1: Positive edge, -1: Negative edge)", "edges")
	)

	# Initialize computed columns so that they are counted as present in the dataset
	data["combinedCoefficients"] = None
	if not config["noPUC"]:
		data["expectedEdgeFilterInt"] = None
		data["nonPucPassed"] = None
	data["medianValueE"] = None
	data["medianValueT"] = None
	data["medianValueET"] = None
	data["meanValueE"] = None
	data["meanValueT"] = None
	data["meanValueET"] = None

	dataKeys = csvConfig.getDataKeys()
	# Add keys for the data the computed columns depend on to ensure its existence
	dataKeys.append("originalData")
	if not config["noPUC"]:
		dataKeys.append("expectedEdgeFilter")

	# Check that all of the required data is present in the dataset
	missingKeys = [key for key in dataKeys if key not in data]
	if len(missingKeys) > 0:
		raise MissingDataError(missingKeys)

	# Determine computed columns
	data["combinedCoefficients"] = data["correlationCoefficients"].mean(dim="metatreatment")
	if not config["noPUC"]:
		data["expectedEdgeFilterInt"] = data["expectedEdgeFilter"].astype(int)
		data["nonPucPassed"] = data["diagonalFilter"] & data["individualCorrelationPValueFilter"] & data["combinedCorrelationPValueFilter"] & data["correctedCorrelationPValueFilter"] & data["correlationFilter"]
		if "correlationCoefficientFilter" in data:
			data["nonPucPassed"] &= data["correlationCoefficientFilter"]

	data["medianValueE"] = data["originalData"].groupby("experiment").map(lambda a: a.median(dim="organism"))
	data["medianValueT"] = data["originalData"].groupby("treatment").map(lambda a: a.median(dim="organism"))
	data["medianValueET"] = data["originalData"].groupby("experiment").map(lambda a: a.groupby("treatment").map(lambda a: a.median(dim="organism")))

	data["meanValueE"] = data["originalData"].groupby("experiment").map(lambda a: a.mean(dim="organism"))
	data["meanValueT"] = data["originalData"].groupby("treatment").map(lambda a: a.mean(dim="organism"))
	data["meanValueET"] = data["originalData"].groupby("experiment").map(lambda a: a.groupby("treatment").map(lambda a: a.mean(dim="organism")))

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
	combinedDifferencePValuesStacked = data["combinedDifferencePValues"].stack({"measurableAndCellType": ("measurable", "cellType")})
	# Reset index to match the format of the rest of the data
	combinedDifferencePValuesNoNans = combinedDifferencePValuesStacked.dropna("measurableAndCellType").reset_index("measurableAndCellType")
	data["combinedDifferencePValuesStacked"] = combinedDifferencePValuesStacked.reset_index("measurableAndCellType")
	data["correctedDifferencePValuesStacked"] = data["correctedDifferencePValues"].stack({"measurableAndCellType": ("measurable", "cellType")}).reset_index("measurableAndCellType")

	if nodesOnly:
		combinedDifferencePValues = combinedDifferencePValuesStacked.reset_index("measurableAndCellType")
		cellTypeCoords = [data["filteredData"].sel(measurableAndCellType=m).cellType.item() for m in data["filteredData"].coords["measurableAndCellType"].data]
		measurableCoords = [data["filteredData"].sel(measurableAndCellType=m).measurable.item() for m in data["filteredData"].coords["measurableAndCellType"].data]

		coordArr = data["filteredData"]
		indexList = [numpy.argwhere(((combinedDifferencePValues.cellType==cellType) & (combinedDifferencePValues.measurable==measurable)).data).item() for cellType, measurable in zip(cellTypeCoords, measurableCoords)]
	else:
		coordArr = combinedDifferencePValuesNoNans
		indexList = numpy.argwhere(~numpy.isnan(combinedDifferencePValuesStacked.data)).flatten()

	data["foldChangesStacked"] = data["foldChanges"].stack({"measurableAndCellType": ("measurable", "cellType")}).reset_index("measurableAndCellType")
	data["foldChangeSignsStacked"] = data["foldChangeSigns"].stack({"measurableAndCellType": ("measurable", "cellType")}).reset_index("measurableAndCellType")
	data["foldChangeFilterStacked"] = data["foldChangeFilter"].stack({"measurableAndCellType": ("measurable", "cellType")}).reset_index("measurableAndCellType")

	csvConfig = CsvWriter.Config("measurableAndCellType",
		CsvWriter.CoordinateFunction("Measurable", lambda m: "{}".format(coordArr.coords["measurable"][m].item())),
		CsvWriter.CoordinateFunction("Cell Type", lambda m: "{}".format(coordArr.coords["cellType"][m].item())),
		CsvWriter.CoordinateFunction("Measurable Type", lambda m: "{}".format(coordArr.coords["measurableType"][m].item())),
		CsvWriter.Per("Average Value ({})", "stacked", "organism", coordMap=indexList),
		CsvWriter.Per("Average Log2 Fold Change ({})", "foldChangesStacked", "experiment", coordMap=indexList),
		CsvWriter.Per("Average Log2 Fold Change Direction ({})", "foldChangeSignsStacked", "experiment", coordMap=indexList),
		CsvWriter.Column("Average Log2 Fold Changes Consistent", "foldChangeFilterStacked", coordMap=indexList),
		CsvWriter.Per("Corrected Comparison p-value ({})", "correctedDifferencePValuesStacked", "experiment", coordMap=indexList),
		CsvWriter.Column("Combined Comparison p-value", "combinedDifferencePValuesStacked", coordMap=indexList)
	)

	dataKeys = csvConfig.getDataKeys()
	dataKeys.append("combinedDifferencePValues")
	dataKeys.append("foldChanges")
	dataKeys.append("foldChangeSigns")
	dataKeys.append("foldChangeFilter")
	missingKeys = [key for key in dataKeys if key not in data]
	if len(missingKeys) > 0:
		raise MissingDataError(missingKeys)

	CsvWriter.writeCsv(filePath, csvConfig, data, coordArr.coords["measurableAndCellType"].data)

def writeComparisonsSingleCell(data, config, outDir):
	fileName = "all_comparisons.csv"
	try:
		writeMeasurableCsvSingleCell(data, config, outDir / fileName, False)
	except MissingDataError as e:
		print(f"WARNING: {e.message}, so {fileName} could not be created")

def writeNodesSingleCell(data, config, outDir):
	fileName = "node_comparisons.csv"
	try:
		writeMeasurableCsvSingleCell(data, config, outDir / fileName, True)
	except MissingDataError as e:
		print(f"WARNING: {e.message}, so {fileName} could not be created")

def writeEdgeCsvSingleCell(data, config, filePath, finalOnly=False):
	if finalOnly:
		allEdges = list(itertools.combinations(data["filteredData"].coords["measurableAndCellType"].data, 2))
		edgeList = [edge for edge in allEdges if data["edges"].sel(measurableAndCellType1=edge[0], measurableAndCellType2=edge[1]) != 0]
	else:
		edgeList = list(itertools.combinations(data["filteredData"].coords["measurableAndCellType"].data, 2))

	cellTypeCoords = [data["filteredData"].sel(measurableAndCellType=m).cellType.item() for m in data["filteredData"].coords["measurableAndCellType"].data]
	measurableCoords = [data["filteredData"].sel(measurableAndCellType=m).measurable.item() for m in data["filteredData"].coords["measurableAndCellType"].data]
	indexList = [numpy.argwhere(((data["stacked"].cellType==cellType) & (data["stacked"].measurable==measurable)).data).item() for cellType, measurable in zip(cellTypeCoords, measurableCoords)]

	data["foldChangesStacked"] = data["foldChanges"].stack({"measurableAndCellType": ("measurable", "cellType")}).reset_index("measurableAndCellType")
	data["combinedFoldChangeSignsStacked"] = data["combinedFoldChangeSigns"].stack({"measurableAndCellType": ("measurable", "cellType")}).reset_index("measurableAndCellType")

	csvConfig = CsvWriter.Config(["measurableAndCellType1", "measurableAndCellType2"],
		CsvWriter.Property("Measurable 1", "correlationCoefficients", "measurable1"),
		CsvWriter.Property("Measurable 2", "correlationCoefficients", "measurable2"),
		CsvWriter.PropertiesFormatted("Measurables", "correlationCoefficients", "{}<==>{}", ["measurable1", "measurable2"]),
		CsvWriter.Property("Cell Type 1", "correlationCoefficients", "cellType1"),
		CsvWriter.Property("Cell Type 2", "correlationCoefficients", "cellType2"),
		CsvWriter.PropertiesFormatted("Cell Types", "correlationCoefficients", "{}<==>{}", ["cellType1", "cellType2"]),
		CsvWriter.CoordComponentPer("Partner 1 Average Value ({})", "stacked", 0, "measurableAndCellType", "organism", coordMap=indexList),
		CsvWriter.CoordComponentPer("Partner 1 Log2 Fold Change {}", "foldChangesStacked", 0, "measurableAndCellType", "experiment", coordMap=indexList),
		CsvWriter.CoordComponentColumn("Partner 1 Fold Change Direction", "combinedFoldChangeSignsStacked", 0, "measurableAndCellType", coordMap=indexList),
		CsvWriter.CoordComponentPer("Partner 2 Average Value ({})", "stacked", 1, "measurableAndCellType", "organism", coordMap=indexList),
		CsvWriter.CoordComponentPer("Partner 2 Log2 Fold Change ({})", "foldChangesStacked", 1, "measurableAndCellType", "experiment", coordMap=indexList),
		CsvWriter.CoordComponentColumn("Partner 2 Fold Change Direction", "combinedFoldChangeSignsStacked", 1, "measurableAndCellType", coordMap=indexList),
		*listItemIf(CsvWriter.Column("Fold Change Signs Match", "foldChangeSignProducts"), not config["noPUC"]),
		CsvWriter.Per("Correlation p-value ({})", "correlationPValues", "metatreatment"),
		CsvWriter.Per("Correlation Coefficient ({})", "correlationCoefficients", "metatreatment"),
		CsvWriter.Column("Correlation Coefficients Consistent", "correlationFilter"),
		CsvWriter.Column("Combined Correlation p-value", "combinedCorrelationPValues"),
		CsvWriter.Column("Corrected Correlation p-value", "correctedCorrelationPValues"),
		*listItemIf(CsvWriter.Column("PUC", "expectedEdgeFilter"), not config["noPUC"]),
		CsvWriter.Column("Final Network Value (0: No edge, 1: Positive edge, -1: Negative edge)", "edges")
	)

	dataKeys = csvConfig.getDataKeys()
	dataKeys.append("foldChanges")
	dataKeys.append("combinedFoldChangeSigns")

	missingKeys = [key for key in dataKeys if key not in data]
	if len(missingKeys) > 0:
		raise MissingDataError(missingKeys)

	CsvWriter.writeCsv(filePath, csvConfig, data, edgeList)

def writeCorrelationsSingleCell(data, config, outDir):
	fileName = "correlations_bw_signif_measurables.csv"
	try:
		writeEdgeCsvSingleCell(data, config, outDir / fileName, False)
	except MissingDataError as e:
		print(f"WARNING: {e.message}, so {fileName} could not be created")

def writeSummarySingleCell(data, config, outDir):
	fileName = "network_output_comp.csv"
	try:
		writeEdgeCsvSingleCell(data, config, outDir / fileName, True)
	except MissingDataError as e:
		print(f"WARNING: {e.message}, so {fileName} could not be created")

def writeConfigValues(config, outDir):
	configValuesFilePath = outDir / "config_values.txt"
	configValuesFile = open(configValuesFilePath, "w")

	def writeThresholds(thresholdKey):
		pValueThresholds = config[thresholdKey]

		if not isinstance(pValueThresholds, dict):
			configValuesFile.write("\tall: {}\n".format(pValueThresholds))
			return

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

	configValuesFile.write("Comparison p-value thresholds:\n")
	writeThresholds("differencePValueThresholds")

	configValuesFile.write("\n")

	configValuesFile.write("Correlation p-value thresholds:\n")
	writeThresholds("correlationPValueThresholds")

	coefficientThresholds = config["correlationCoefficientThresholds"]
	if coefficientThresholds is not None:
		configValuesFile.write("\nCorrelation coefficient thresholds:\n")
		if not isinstance(coefficientThresholds, dict):
			configValuesFile.write("\tall: {}\n".format(coefficientThresholds))
		else:
			for type_, threshold in coefficientThresholds.items():
				configValuesFile.write("\t{}: {}\n".format(type_, threshold))

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
		try:
			nodeCsvConfig, data = setupMeasurableCsv(data, config)
		except MissingDataError as e:
			print(f"WARNING: {e.message}, so node CSV files could not be created")
		else:
			writeComparisons(data, config, nodeCsvConfig, args.outDir)
			writeNodes(data, config, nodeCsvConfig, args.outDir)

		try:
			edgeCsvConfig, data = setupEdgeCsv(data, config)
		except MissingDataError as e:
			print(f"WARNING: {e.message}, so edge CSV files could not be created")
		else:
			writeCorrelations(data, config, edgeCsvConfig, args.outDir)
			writeSummary(data, config, edgeCsvConfig, args.outDir)

	writeConfigValues(config, args.outDir)

