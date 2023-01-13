from zipfile import ZipFile
import xarray

def readDataZip(dataZipPath):
	data = {}
	dataZip = ZipFile(dataZipPath)
	for dataArrayName in dataZip.namelist():
		dataArrayFile = dataZip.open(dataArrayName)
		data[dataArrayName] = xarray.open_dataarray(dataArrayFile)
	dataZip.close()
	return data

