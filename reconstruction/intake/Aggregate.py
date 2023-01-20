from functools import partial
from io import StringIO
import json
import numpy
import pandas
import xarray

from .util import *

def intakeAggregateData(dataDir):
	metadata = json.loads(readAndDecodeFile(dataDir / "metadata.json"))
	measurableTypeMapString = readAndDecodeFile(dataDir / metadata["measurableTypeMapFile"])
	if metadata["measurableTypeMapFile"].endswith(".json"):
		measurableTypeMap = json.loads(measurableTypeMapString)
	elif metadata["measurableTypeMapFile"].endswith(".csv"):
		measurableTypeMap = readClassificationCsv(measurableTypeMapString)

	nextOrganismId = 0

	allExperimentData = []
	for experimentMetadata in metadata["experiments"]:
		experimentName = experimentMetadata["name"]

		experimentDataString = readAndDecodeFile(dataDir / experimentMetadata["dataFile"])
		experimentDataDataframe = pandas.read_csv(StringIO(experimentDataString), index_col=0)
		experimentDataDataframe = experimentDataDataframe.replace(["na", "NA", "n/a", "N/A"], numpy.nan).apply(partial(pandas.to_numeric, errors="ignore"))

		treatmentMapString = readAndDecodeFile(dataDir / experimentMetadata["treatmentMapFile"])
		if experimentMetadata["treatmentMapFile"].endswith(".json"):
			treatmentMap = json.loads(treatmentMapString)
		elif experimentMetadata["treatmentMapFile"].endswith(".csv"):
			treatmentMap = readClassificationCsv(treatmentMapString)

		experimentData = xarray.DataArray(experimentDataDataframe)
		experimentData = experimentData.rename({"dim_1": "organism", "ID": "measurable"})

		organismCoords = experimentData.coords["organism"]

		experimentCoords = [experimentName] * len(organismCoords)

		inverseTreatmentMap = {}
		for treatment, organisms in treatmentMap.items():
			for organism in organisms:
				inverseTreatmentMap[organism] = treatment
		treatmentCoords = [inverseTreatmentMap[organism.item()] for organism in organismCoords]

		organismCoords = [str(i) for i in range(nextOrganismId, nextOrganismId + len(organismCoords))]
		nextOrganismId += len(organismCoords)

		experimentData = experimentData.assign_coords({"organism": organismCoords, "treatment": ("organism", treatmentCoords), "experiment": ("organism", experimentCoords)})
		allExperimentData.append(experimentData)

	data = xarray.concat(allExperimentData, dim="organism")

	inverseMeasurableTypeMap = {}
	for measurableType, measurables in measurableTypeMap.items():
		for measurable in measurables:
			inverseMeasurableTypeMap[measurable] = measurableType

	measurableCoords = data.coords["measurable"]
	measurableTypeCoords = [inverseMeasurableTypeMap[measurable.item()] for measurable in measurableCoords]
	data = data.assign_coords({"measurableType": ("measurable", measurableTypeCoords)})

	return data
