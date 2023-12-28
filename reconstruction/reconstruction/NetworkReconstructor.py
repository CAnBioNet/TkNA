from collections import defaultdict
import networkx
import pandas
import tempfile
import xarray
from zipfile import ZipFile

# Checks that the treatments specified in the config file are present in the relevant experiments
def verifyTreatmentsPresent(data, config):
	metatreatments = config["metatreatments"]
	if metatreatments is not None:
		metatreatmentETMap = defaultdict(set)
		for metatreatmentMembers in metatreatments.values():
			for experimentAndTreatment in metatreatmentMembers:
				experiment = experimentAndTreatment[0]
				treatment = experimentAndTreatment[1]
				metatreatmentETMap[experiment].add(treatment)

	networkTreatment = config["networkTreatment"]
	groupedByExperiment = data.groupby("experiment")
	for experiment in groupedByExperiment:
		experimentName = experiment[0]
		experimentArray = experiment[1]

		if "comparisonTreatments" in config:
			for treatment in config["comparisonTreatments"]:
				if treatment not in experimentArray.treatment:
					raise Exception(f"Treatment {treatment} specified in \"comparisonTreatments\" is not present in experiment {experimentName}")

		if networkTreatment is not None and networkTreatment not in experimentArray.treatment:
			raise Exception(f"Treatment {networkTreatment} specified in \"networkTreatment\" is not present in experiment {experimentName}")
		if metatreatments is not None:
			for treatment in metatreatmentETMap[experimentName]:
				if treatment not in experimentArray.treatment:
					raise Exception(f"Treatment {treatment} specified in \"metatreatments\" is not present in experiment {experimentName}")

# Adapted from https://github.com/pydata/xarray/pull/7689
def resetDataarrayEncoding(da):
	ds = da._to_temp_dataset()
	reset_variables = {k: v._replace(encoding={}) for k, v, in ds.variables.items()}
	ds = ds._replace(variables=reset_variables, encoding={})
	return da._from_temp_dataset(ds)

class NetworkReconstructor:
	def runPipeline(self, stages, allData={}, start=None, stopStage=None, dataOutFilePath=None):
		started = True
		if start is not None:
			started = False
			dataInZip = ZipFile(start["startData"])
			for dataArrayName in dataInZip.namelist():
				dataArrayFile = dataInZip.open(dataArrayName)
				allData[dataArrayName] = xarray.open_dataarray(dataArrayFile)
			dataInZip.close()

		for stage in stages:
			if not started and stage.__name__ == start["startStage"]:
				started = True
			if started:
				stage(allData)
			if stopStage is not None and stage.__name__ == stopStage:
				break

		if dataOutFilePath is not None:
			dataOutFilePath.parent.mkdir(parents=True, exist_ok=True)
			dataOutZip = ZipFile(dataOutFilePath, mode="w")
			for dataName, dataArray in allData.items():
				tempFile = tempfile.NamedTemporaryFile()
				dims = list(dataArray.dims)
				for dim in dims:
					if isinstance(dataArray.get_index(dim), pandas.MultiIndex):
						dataArray = dataArray.reset_index(dim)
				dataArray = resetDataarrayEncoding(dataArray) # See https://github.com/pydata/xarray/issues/6352
				dataArray.to_netcdf(tempFile.name)
				dataOutZip.write(tempFile.name, arcname=dataName)
				tempFile.close()
			dataOutZip.close()

		return allData

	def reconstructNetwork(self, **kwargs):
		raise NotImplementedError()

