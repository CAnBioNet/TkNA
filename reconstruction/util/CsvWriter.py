import csv
from itertools import chain
import numpy

def getArr(data, dataKey, dim):
	dataArr = data[dataKey]
	if isinstance(dim, list):
		newDim = "_".join(dim)
		dataArr = dataArr.stack({newDim: dim})
		dim = newDim
	return dataArr, dim

class Column:
	def __init__(self, title, dataKey):
		self.title = title
		self.dataKey = dataKey

	def getHeaders(self, data):
		return [self.title]

	def getValues(self, data, dim, coords):
		dataArr, dim = getArr(data, self.dataKey, dim)
		return [dataArr.sel({dim: coords}).data]

	def getDataKeys(self):
		return [self.dataKey]

class Property:
	def __init__(self, title, dataKey, propertyDim):
		self.title = title
		self.dataKey = dataKey
		self.propertyDim = propertyDim

	def getHeaders(self, data):
		return [self.title]

	def getValues(self, data, dim, coords):
		dataArr, dim = getArr(data, self.dataKey, dim)
		return [dataArr.sel({dim: coords}).coords[self.propertyDim].data]

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

	def getValues(self, data, dim, coords):
		dataArr, dim = getArr(data, self.dataKey, dim)
		properties = numpy.array([dataArr.sel({dim: coords}).coords[propertyDim].data for propertyDim in self.propertyDims])
		return [numpy.apply_along_axis(lambda properties: self.formatStr.format(properties))]

	def getDataKeys(self):
		return [self.dataKey]

class Per:
	def __init__(self, titleTemplate, dataKey, perDim):
		self.titleTemplate = titleTemplate
		self.dataKey = dataKey
		self.perDim = perDim

	def getHeaders(self, data):
		self.perCoords = data[self.dataKey].coords[self.perDim].data # Save coordinates to keep order consistent with values
		return [self.titleTemplate.format(coord) for coord in self.perCoords]

	def getValues(self, data, dim, coords):
		dataArr, dim = getArr(data, self.dataKey, dim)
		return [dataArr.sel({self.perDim: perCoord, dim: coords}).data for perCoord in self.perCoords]

	def getDataKeys(self):
		return [self.dataKey]

class Coordinate:
	def __init__(self, title):
		self.title = title

	def getHeaders(self, data):
		return [self.title]

	def getValues(self, data, dim, coords):
		return [coords]

	def getDataKeys(self):
		return []

class CoordinateFormatted:
	def __init__(self, title, valueTemplate):
		self.title = title
		self.valueTemplate = valueTemplate

	def getHeaders(self, data):
		return [self.title]

	def getValues(self, data, dim, coords):
		if isinstance(dim, list):
			return [[self.valueTemplate.format(*coord) for coord in coords]]
		else:
			return [[self.valueTemplate.format(coord) for coord in coords]]

	def getDataKeys(self):
		return []

class CoordinateFunction:
	def __init__(self, title, func):
		self.title = title
		self.func = func

	def getHeaders(self, data):
		return [self.title]

	def getValues(self, data, dim, coords):
		if isinstance(dim, list):
			return [[self.func(*coord) for coord in coords]]
		else:
			return [[self.func(coord) for coord in coords]]

	def getDataKeys(self):
		return []

class CoordComponent:
	def __init__(self, title, componentIndex):
		self.title = title
		self.componentIndex = componentIndex

	def getHeaders(self, data):
		return [self.title]

	def getValues(self, data, dim, coords):
		return [[coord[self.componentIndex] for coord in coords]]

	def getDataKeys(self):
		return []

class CoordComponentColumn:
	def __init__(self, title, dataKey, componentIndex, componentDim):
		self.title = title
		self.dataKey = dataKey
		self.componentIndex = componentIndex
		self.componentDim = componentDim

	def getHeaders(self, data):
		return [self.title]

	def getValues(self, data, dim, coords):
		dataArr = data[self.dataKey]
		components = [coord[self.componentIndex] for coord in coords]
		return [dataArr.sel({self.componentDim: components}).data]

	def getDataKeys(self):
		return [self.dataKey]

class CoordComponentPer:
	def __init__(self, titleTemplate, dataKey, componentIndex, componentDim, perDim):
		self.titleTemplate = titleTemplate
		self.dataKey = dataKey
		self.componentIndex = componentIndex
		self.componentDim = componentDim
		self.perDim = perDim

	def getHeaders(self, data):
		self.perCoords = data[self.dataKey].coords[self.perDim].data # Save coordinates to keep order consistent with values
		return [self.titleTemplate.format(coord) for coord in self.perCoords]

	def getValues(self, data, dim, coords):
		dataArr = data[self.dataKey]
		components = [coord[self.componentIndex] for coord in coords]
		return [dataArr.sel({self.perDim: perCoord, self.componentDim: components}).data for perCoord in self.perCoords]

	def getDataKeys(self):
		return [self.dataKey]

class CoordComponentPropertyFormatted:
	def __init__(self, title, dataKey, formatStr, componentDim, propertyDim):
		self.title = title
		self.dataKey = dataKey
		self.formatStr = formatStr
		self.componentDim = componentDim
		self.propertyDim = propertyDim

	def getHeaders(self, data):
		return [self.title]

	def getValues(self, data, dim, coords):
		dataArr = data[self.dataKey]
		components = numpy.array(coords).T
		componentProperties = numpy.array([dataArr.sel({self.componentDim: componentList}).coords[self.propertyDim].data for componentList in components])
		return [numpy.apply_along_axis(lambda properties: self.formatStr.format(*properties), 0, componentProperties)]

	def getDataKeys(self):
		return [self.dataKey]

class Config:
	def __init__(self, rowDim, *columns):
		self.rowDim = rowDim
		self.columns = columns

	def getDataKeys(self):
		return list(chain(*[col.getDataKeys() for col in self.columns]))

def writeCsv(filePath, config, data, rowCoords):
	header = list(chain(*[col.getHeaders(data) for col in config.columns]))
	values = numpy.stack(list(chain(*[col.getValues(data, config.rowDim, rowCoords) for col in config.columns])))

	csvFile = open(filePath, "w")
	csvWriter = csv.writer(csvFile)

	csvWriter.writerow(header)
	numpy.apply_along_axis(lambda row: csvWriter.writerow(row), 0, values)

	csvFile.close()

