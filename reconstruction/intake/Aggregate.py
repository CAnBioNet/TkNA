import json
import pandas
import xarray

from .util import readClassificationCsv

def intakeAggregateData(dataDir):
	with open(dataDir / "metadata.json") as metadataFile:
		metadata = json.load(metadataFile)
	with open(dataDir / metadata["measurableTypeMapFile"]) as measurableTypeMapFile:
		if metadata["measurableTypeMapFile"].endswith(".json"):
			measurableTypeMap = json.load(measurableTypeMapFile)
		elif metadata["measurableTypeMapFile"].endswith(".csv"):
			measurableTypeMap = readClassificationCsv(measurableTypeMapFile)

	nextOrganismId = 0

	allExperimentData = []
	for experimentMetadata in metadata["experiments"]:
		experimentName = experimentMetadata["name"]

		with open(dataDir / experimentMetadata["dataFile"]) as experimentDataFile:
			experimentDataDataframe = pandas.read_csv(experimentDataFile, index_col=0)
		with open(dataDir / experimentMetadata["treatmentMapFile"]) as treatmentMapFile:
			if experimentMetadata["treatmentMapFile"].endswith(".json"):
				treatmentMap = json.load(treatmentMapFile)
			elif experimentMetadata["treatmentMapFile"].endswith(".csv"):
				treatmentMap = readClassificationCsv(treatmentMapFile)

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
