from functools import partial
import glob
from io import StringIO
import json
import numpy
import pandas
from pathlib import Path
import xarray

from .util import *
from util import Dataset

def intakeSingleCellData(dataDir):
	metadata = json.loads(readAndDecodeFile(dataDir / "metadata.json"))
	cellTypes = metadata["cellTypes"]

	measurableTypeMapPaths = list(dataDir.glob("measurable_type_map.*"))
	if len(measurableTypeMapPaths) == 0:
		raise Exception("No measurable type map found. Expected as [dataDir]/measurable_type_map.[json|csv].")
	measurableTypeMapPath = measurableTypeMapPaths[0]
	measurableTypeMapString = readAndDecodeFile(measurableTypeMapPath)
	if measurableTypeMapPath.suffix == ".json":
		measurableTypeMap = json.loads(measurableTypeMapString)
	elif measurableTypeMapPath.suffix == ".csv":
		measurableTypeMap = readClassificationCsv(measurableTypeMapString)
	else:
		raise Exception("Unrecognized file extension for measurable type map")
	inverseMeasurableTypeMap = {}
	for measurableType, measurables in measurableTypeMap.items():
		for measurable in measurables:
			inverseMeasurableTypeMap[measurable] = measurableType

	experimentNames = [experimentMetadata["name"] for experimentMetadata in metadata["experiments"]]
	measurables = list(inverseMeasurableTypeMap.keys())
	nans = numpy.empty((len(experimentNames), len(cellTypes), len(measurables)))
	nans[:] = numpy.nan
	pValues = xarray.DataArray(nans.copy(), dims=["experiment", "cellType", "measurable"], coords={"experiment": experimentNames, "cellType": cellTypes, "measurable": measurables})
	foldChanges = xarray.DataArray(nans.copy(), dims=["experiment", "cellType", "measurable"], coords={"experiment": experimentNames, "cellType": cellTypes, "measurable": measurables})

	allExperimentData = []
	for experimentMetadata in metadata["experiments"]:
		experimentName = experimentMetadata["name"]
		experimentDir = dataDir / experimentMetadata["dataDir"]
		if not experimentDir.exists():
			raise Exception("Experiment directory {} does not exist".format(experimentMetadata["dataDir"]))
		elif not experimentDir.is_dir():
			raise Exception("Experiment direectory {} is not a directory".format(experimentMetadata["dataDir"]))

		treatmentMapPaths = list(experimentDir.glob("treatment_map.*"))
		if len(treatmentMapPaths) == 0:
			raise Exception("No treatment map found in experiment \"{}\". Expected as [experimentDir]/treatment_map.[json|csv].".format(experimentName))
		treatmentMapPath = treatmentMapPaths[0]
		treatmentMapString = readAndDecodeFile(treatmentMapPath)
		if treatmentMapPath.suffix == ".json":
			treatmentMap = json.loads(treatmentMapString)
		elif treatmentMapPath.suffix == ".csv":
			treatmentMap = readClassificationCsv(treatmentMapString)
		else:
			raise Exception("Unrecognized file extension for experiment \"{}\" treatment map".format(experimentName))
		inverseTreatmentMap = {}
		for treatment, organisms in treatmentMap.items():
			for organism in organisms:
				inverseTreatmentMap[organism] = treatment

		experimentData = []
		for cellType in cellTypes:
			cellTypeDir = experimentDir / cellType
			if not cellTypeDir.exists():
				raise Exception("Directory for cell type \"{}\" in experiment \"{}\" does not exist".format(cellType, experimentName))
			elif not cellTypeDir.is_dir():
				raise Exception("Directory for cell type \"{}\" in experiment \"{}\" is not a directory".format(cellType, experimentName))

			organismFileMapPaths = list(cellTypeDir.glob("organism_file_map.*"))
			if len(organismFileMapPaths) == 0:
				raise Exception("No organism-file map found in {}. Expected as [cellTypeDir]/organism_file_map.[json|csv].".format(cellTypeDir))
			organismFileMapPath = organismFileMapPaths[0]
			organismFileMapString = readAndDecodeFile(organismFileMapPath)
			if organismFileMapPath.suffix == ".json":
				organismFileMap = json.loads(organismFileMapString)
			elif organismFileMapPath.suffix == ".csv":
				organismFileMap = readClassificationCsv(organismFileMapString)
			else:
				raise Exception("Unrecognized file extension for cell type \"{}\", experiment \"{}\" organism-file map".format(cellType, experimentName))

			for organism, organismFileName in organismFileMap.items():
				cellTypeDataString = readAndDecodeFile(cellTypeDir / organismFileName)
				cellTypeDataFrame = pandas.read_csv(StringIO(cellTypeDataString), index_col=0)
				cellTypeDataFrame = cellTypeDataFrame.replace(["na", "NA", "n/a", "N/A"], numpy.nan).apply(partial(pandas.to_numeric, errors="ignore"))
				cellTypeData = xarray.DataArray(cellTypeDataFrame)
				cellTypeData = cellTypeData.rename({cellTypeData.dims[0]: "measurable", cellTypeData.dims[1]: "cell"})

				cellCoords = ["{}_{}".format(cell, organism) for cell in cellTypeData.coords["cell"].data] # Ensure all cells have unique coords
				organismCoords = ["{}_{}".format(experimentName, organism)] * cellTypeData.sizes["cell"]
				cellTypeCoords = [cellType] * cellTypeData.sizes["cell"]
				treatmentCoords = [inverseTreatmentMap[organism]] * cellTypeData.sizes["cell"]

				cellTypeData = cellTypeData.assign_coords({"cell": cellCoords, "organism": ("cell", organismCoords), "cellType": ("cell", cellTypeCoords), "treatment": ("cell", treatmentCoords)})
				experimentData.append(cellTypeData)

			diffFilePath = cellTypeDir / "diff.csv"
			if not diffFilePath.exists():
				raise Exception("Diff file does not exist for cell type \"{}\" in experiment \"{}\". Expected as [cellTypeDir]/diff.csv.".format(cellType, experimentName))
			diffString = readAndDecodeFile(diffFilePath)
			differentials = pandas.read_csv(StringIO(diffString), index_col=0)
			differentials = differentials.replace(["na", "NA", "n/a", "N/A"], numpy.nan).apply(partial(pandas.to_numeric, errors="ignore"))
			missingMeasurables = set(differentials.axes[0]).difference(set(pValues.coords["measurable"].data))
			if len(missingMeasurables) != 0:
				raise Exception("Measurables in {} not found in organism data: {}".format(diffFilePath, missingMeasurables))
			differentialIndex = {"experiment": experimentName, "cellType": cellType, "measurable": list(differentials.axes[0])}
			pValues.loc[differentialIndex] = differentials["p_val_adj"]
			foldChanges.loc[differentialIndex] = differentials["avg_log2FC"]

		combinedExperimentData = xarray.concat(experimentData, dim="cell")
		experimentCoords = [experimentName] * combinedExperimentData.sizes["cell"]
		combinedExperimentData = combinedExperimentData.assign_coords({"experiment": ("cell", experimentCoords)})
		allExperimentData.append(combinedExperimentData)

	data = xarray.concat(allExperimentData, dim="cell")

	def addMeasurableTypeCoords(table):
		measurableCoords = table.coords["measurable"]
		typeCoords = [inverseMeasurableTypeMap[measurable.item()] for measurable in measurableCoords]
		return table.assign_coords({"measurableType": ("measurable", typeCoords)})
	data = addMeasurableTypeCoords(data)
	pValues = addMeasurableTypeCoords(pValues)
	foldChanges = addMeasurableTypeCoords(foldChanges)

	dataset = Dataset()
	dataset.add_table("cellData", data)
	dataset.add_table("correctedDifferencePValues", pValues)
	dataset.add_table("foldChanges", foldChanges)

	return dataset

