import networkx
import pandas
import tempfile
import xarray
from zipfile import ZipFile

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

