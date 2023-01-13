import csv
from itertools import chain

# TODO: Always use tuples of dimensions/coordinates

class DimTuple:
	def __init__(self, *dims):
		self.dims = dims

	def __getitem__(self, item):
		return self.dims[item]

class Column:
	def __init__(self, title, dataKey, rowCoordConversionDict=None, missingOk=False):
		self.title = title
		self.dataKey = dataKey
		self.missingOk = missingOk
		self.rowCoordConversionDict = rowCoordConversionDict

	def getHeaders(self, data):
		if self.dataKey not in data:
			if self.missingOk:
				return []
			else:
				raise Exception("Specified data array not present")
		return [self.title]

	def getValues(self, data, rowDim, rowCoord):
		if self.dataKey not in data:
			if self.missingOk:
				return []
			else:
				raise Exception("Specified data array not present")
		if self.rowCoordConversionDict is not None:
			rowCoord = self.rowCoordConversionDict[rowCoord]
		if isinstance(rowDim, DimTuple):
			return [data[self.dataKey].sel({dim: coord for dim, coord in zip(rowDim.dims, rowCoord)}).item()]
		else:
			return [data[self.dataKey].sel({rowDim: rowCoord}).item()]

	def getDataKeys(self):
		return [self.dataKey]

class Property:
	def __init__(self, title, dataKey, propertyDim):
		self.title = title
		self.dataKey = dataKey
		self.propertyDim = propertyDim

	def getHeaders(self, data):
		return [self.title]

	def getValues(self, data, rowDim, rowCoord):
		if isinstance(rowDim, DimTuple):
			return [data[self.dataKey].sel({dim: coord for dim, coord in zip(rowDim.dims, rowCoord)}).coords[self.propertyDim].item()]
		else:
			return [data[self.dataKey].sel({rowDim: rowCoord}).coords[self.propertyDim].item()]

	def getDataKeys(self):
		return [self.dataKey]

class PropertiesFormatted:
	def __init__(self, title, dataKey, formatStr, propertyDims):
		self.title = title
		self.dataKey = dataKey
		self.formatStr = formatStr
		self.propertyDims = propertyDims

	def getHeaders(self, data):
		return [self.title]

	def getValues(self, data, rowDim, rowCoord):
		if isinstance(rowDim, DimTuple):
			value = data[self.dataKey].sel({dim: coord for dim, coord in zip(rowDim.dims, rowCoord)})
		else:
			value = data[self.dataKey].sel({rowDim: rowCoord})
		return [self.formatStr.format(*[value.coords[propertyDim].item() for propertyDim in self.propertyDims])]

	def getDataKeys(self):
		return [self.dataKey]

class Per:
	def __init__(self, titleTemplate, dataKey, perDim, rowCoordConversionDict=None):
		self.titleTemplate = titleTemplate
		self.dataKey = dataKey
		self.perDim = perDim
		self.rowCoordConversionDict = rowCoordConversionDict

	def getHeaders(self, data):
		self.perCoords = data[self.dataKey].coords[self.perDim].data # Save coordinates to keep order consistent with values
		return [self.titleTemplate.format(coord) for coord in self.perCoords]

	def getValues(self, data, rowDim, rowCoord):
		if self.rowCoordConversionDict is not None:
			rowCoord = self.rowCoordConversionDict[rowCoord]
		if isinstance(rowDim, DimTuple):
			selDict = {dim: coord for dim, coord in zip(rowDim.dims, rowCoord)}
			return [data[self.dataKey].sel(dict(chain(selDict.items(), {self.perDim: perCoord}.items()))).item() for perCoord in self.perCoords]
		else:
			return [data[self.dataKey].sel({rowDim: rowCoord, self.perDim: perCoord}).item() for perCoord in self.perCoords]

	def getDataKeys(self):
		return [self.dataKey]

class Coordinate:
	def __init__(self, title):
		self.title = title

	def getHeaders(self, data):
		return [self.title]

	def getValues(self, data, rowDim, rowCoord):
		return [rowCoord]

	def getDataKeys(self):
		return []

class CoordinateFormatted:
	def __init__(self, title, valueTemplate):
		self.title = title
		self.valueTemplate = valueTemplate

	def getHeaders(self, data):
		return [self.title]

	def getValues(self, data, rowDim, rowCoord):
		if isinstance(rowDim, DimTuple):
			return [self.valueTemplate.format(*rowCoord)]
		else:
			return [self.valueTemplate.format(rowCoord)]

	def getDataKeys(self):
		return []

class CoordinateFunction:
	def __init__(self, title, func):
		self.title = title
		self.func = func

	def getHeaders(self, data):
		return [self.title]

	def getValues(self, data, rowDim, rowCoord):
		if isinstance(rowDim, DimTuple):
			return [self.func(*rowCoord)]
		else:
			return [self.func(rowCoord)]

	def getDataKeys(self):
		return []

class CoordComponentColumn:
	def __init__(self, title, dataKey, componentIndex, componentDim, rowCoordConversionDict=None):
		self.title = title
		self.dataKey = dataKey
		self.componentIndex = componentIndex
		self.componentDim = componentDim
		self.rowCoordConversionDict = rowCoordConversionDict

	def getHeaders(self, data):
		return [self.title]

	def getValues(self, data, rowDim, rowCoord):
		rowCoord = rowCoord[self.componentIndex]
		if self.rowCoordConversionDict is not None:
			rowCoord = self.rowCoordConversionDict[rowCoord]
		return [data[self.dataKey].sel({self.componentDim: rowCoord}).item()]

	def getDataKeys(self):
		return [self.dataKey]

class CoordComponentPer:
	def __init__(self, titleTemplate, dataKey, componentIndex, componentDim, perDim, rowCoordConversionDict=None):
		self.titleTemplate = titleTemplate
		self.dataKey = dataKey
		self.componentIndex = componentIndex
		self.componentDim = componentDim
		self.perDim = perDim
		self.rowCoordConversionDict = rowCoordConversionDict

	def getHeaders(self, data):
		self.perCoords = data[self.dataKey].coords[self.perDim].data # Save coordinates to keep order consistent with values
		return [self.titleTemplate.format(coord) for coord in self.perCoords]

	def getValues(self, data, rowDim, rowCoord):
		rowCoord = rowCoord[self.componentIndex]
		if self.rowCoordConversionDict is not None:
			rowCoord = self.rowCoordConversionDict[rowCoord]
		return [data[self.dataKey].sel({self.componentDim: rowCoord, self.perDim: perCoord}).item() for perCoord in self.perCoords]

	def getDataKeys(self):
		return [self.dataKey]

class CoordinateValueFunction:
	def __init__(self, title, dataKey, componentDim, func):
		self.title = title
		self.dataKey = dataKey
		self.componentDim = componentDim
		self.func = func

	def getHeaders(self, data):
		return [self.title]

	def getValues(self, data, rowDim, rowCoord):
		return [self.func(*[data[self.dataKey].sel({self.componentDim: component}).item() for component in rowCoord])]

	def getDataKeys(self):
		return [self.dataKey]

class Config:
	def __init__(self, rowDim, *columns):
		self.rowDim = rowDim
		self.columns = columns

	def getDataKeys(self):
		return list(chain(*[col.getDataKeys() for col in self.columns]))

class RowGenerator:
	def __init__(self, config):
		self.config = config

	def getHeader(self, data):
		return chain(*[col.getHeaders(data) for col in self.config.columns])

	def getRow(self, data, rowCoord):
		return chain(*[col.getValues(data, self.config.rowDim, rowCoord) for col in self.config.columns])

	# TODO: Ideally wouldn't need to pass rowCoords, as all datasets used must have the same set of rowCoords.
	# Should do some kind of verification of this and can get the coordinates through that.
	def getRows(self, data, rowCoords):
		for rowCoord in rowCoords:
			yield self.getRow(data, rowCoord)

def writeCsv(filePath, config, data, rowCoords):
	csvFile = open(filePath, "w")
	csvWriter = csv.writer(csvFile)

	rowGenerator = RowGenerator(config)
	csvWriter.writerow(rowGenerator.getHeader(data))
	for row in rowGenerator.getRows(data, rowCoords):
		csvWriter.writerow(row)

	csvFile.close()

