import networkx
import pandas
import tempfile
import xarray
from zipfile import ZipFile

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
				dataArray.to_netcdf(tempFile.name)
				dataOutZip.write(tempFile.name, arcname=dataName)
				tempFile.close()
			dataOutZip.close()

		return allData

	def reconstructNetwork(self, **kwargs):
		raise NotImplementedError()

