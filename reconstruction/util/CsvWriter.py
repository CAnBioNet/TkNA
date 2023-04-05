import csv
from itertools import chain
import numpy

class ColumnSpec:
	def __init__(self, coordMap=None):
		self.coordMap = coordMap
		self.dataKey = None

	def getHeaders(self, data):
		raise NotImplementedError()

	def _getValues(self, data, dim, coords, dataArr=None):
		raise NotImplementedError()

	def getValues(self, data, dim, coords):
		if self.dataKey is None:
			dataArr = None
		else:
			dataArr = data[self.dataKey]

			# If each row corresponds to multiple dimensions, stack them in order to select on them
			if isinstance(dim, list) and all([dimComponent in dataArr.dims for dimComponent in dim]):
				newDim = "_".join(dim)
				dataArr = dataArr.stack({newDim: dim})
				dim = newDim

			# Make it so that sel can be used on coordless dimensions by using indices as coordinates
			# A coordless dimension will not be present in the coords dict, but attempting to access its
			# coords anyway returns a list of indices, so just use that as the coords
			dimsWithoutCoords = [dim for dim in dataArr.dims if (dim not in dataArr.coords)]
			dummyCoords = {dim: dataArr.coords[dim] for dim in dimsWithoutCoords}
			dataArr = dataArr.assign_coords(dummyCoords)

		if self.coordMap is not None:
			if isinstance(coords[0], tuple):
				coords = [[self.coordMap[component] for component in coord] for coord in coords]
			else:
				coords = [self.coordMap[coord] for coord in coords]

		return self._getValues(data, dim, coords, dataArr=dataArr)

class Column(ColumnSpec):
	def __init__(self, title, dataKey, **kwargs):
		super().__init__(**kwargs)
		self.title = title
		self.dataKey = dataKey

	def getHeaders(self, data):
		return [self.title]

	def _getValues(self, data, dim, coords, dataArr):
		return [dataArr.sel({dim: coords}).data]

	def getDataKeys(self):
		return [self.dataKey]

class Property(ColumnSpec):
	def __init__(self, title, dataKey, propertyDim, **kwargs):
		super().__init__(**kwargs)
		self.title = title
		self.dataKey = dataKey
		self.propertyDim = propertyDim

	def getHeaders(self, data):
		return [self.title]

	def _getValues(self, data, dim, coords, dataArr):
		return [dataArr.sel({dim: coords}).coords[self.propertyDim].data]

	def getDataKeys(self):
		return [self.dataKey]

class PropertiesFormatted(ColumnSpec):
	def __init__(self, title, dataKey, formatStr, propertyDims, **kwargs):
		super().__init__(**kwargs)
		self.title = title
		self.dataKey = dataKey
		self.formatStr = formatStr
		self.propertyDims = propertyDims

	def getHeaders(self, data):
		return [self.title]

	def _getValues(self, data, dim, coords, dataArr):
		properties = numpy.array([dataArr.sel({dim: coords}).coords[propertyDim].data for propertyDim in self.propertyDims]).T
		return [[self.formatStr.format(*propertyList) for propertyList in properties]]

	def getDataKeys(self):
		return [self.dataKey]

class Per(ColumnSpec):
	def __init__(self, titleTemplate, dataKey, perDim, **kwargs):
		super().__init__(**kwargs)
		self.titleTemplate = titleTemplate
		self.dataKey = dataKey
		self.perDim = perDim

	def getHeaders(self, data):
		self.perCoords = data[self.dataKey].coords[self.perDim].data # Save coordinates to keep order consistent with values
		return [self.titleTemplate.format(coord) for coord in self.perCoords]

	def _getValues(self, data, dim, coords, dataArr):
		return [dataArr.sel({self.perDim: perCoord, dim: coords}).data for perCoord in self.perCoords]

	def getDataKeys(self):
		return [self.dataKey]

class Coordinate(ColumnSpec):
	def __init__(self, title, **kwargs):
		super().__init__(**kwargs)
		self.title = title

	def getHeaders(self, data):
		return [self.title]

	def _getValues(self, data, dim, coords, dataArr=None):
		return [coords]

	def getDataKeys(self):
		return []

class CoordinateFormatted(ColumnSpec):
	def __init__(self, title, valueTemplate, **kwargs):
		super().__init__(**kwargs)
		self.title = title
		self.valueTemplate = valueTemplate

	def getHeaders(self, data):
		return [self.title]

	def _getValues(self, data, dim, coords, dataArr=None):
		if isinstance(dim, list):
			return [[self.valueTemplate.format(*coord) for coord in coords]]
		else:
			return [[self.valueTemplate.format(coord) for coord in coords]]

	def getDataKeys(self):
		return []

class CoordinateFunction(ColumnSpec):
	def __init__(self, title, func, **kwargs):
		super().__init__(**kwargs)
		self.title = title
		self.func = func

	def getHeaders(self, data):
		return [self.title]

	def _getValues(self, data, dim, coords, dataArr=None):
		if isinstance(dim, list):
			return [[self.func(*coord) for coord in coords]]
		else:
			return [[self.func(coord) for coord in coords]]

	def getDataKeys(self):
		return []

class CoordComponent(ColumnSpec):
	def __init__(self, title, componentIndex, **kwargs):
		super().__init__(**kwargs)
		self.title = title
		self.componentIndex = componentIndex

	def getHeaders(self, data):
		return [self.title]

	def _getValues(self, data, dim, coords, dataArr=None):
		return [[coord[self.componentIndex] for coord in coords]]

	def getDataKeys(self):
		return []

class CoordComponentColumn(ColumnSpec):
	def __init__(self, title, dataKey, componentIndex, componentDim, **kwargs):
		super().__init__(**kwargs)
		self.title = title
		self.dataKey = dataKey
		self.componentIndex = componentIndex
		self.componentDim = componentDim

	def getHeaders(self, data):
		return [self.title]

	def _getValues(self, data, dim, coords, dataArr):
		components = [coord[self.componentIndex] for coord in coords]
		return [dataArr.sel({self.componentDim: components}).data]

	def getDataKeys(self):
		return [self.dataKey]

class CoordComponentPer(ColumnSpec):
	def __init__(self, titleTemplate, dataKey, componentIndex, componentDim, perDim, **kwargs):
		super().__init__(**kwargs)
		self.titleTemplate = titleTemplate
		self.dataKey = dataKey
		self.componentIndex = componentIndex
		self.componentDim = componentDim
		self.perDim = perDim

	def getHeaders(self, data):
		self.perCoords = data[self.dataKey].coords[self.perDim].data # Save coordinates to keep order consistent with values
		return [self.titleTemplate.format(coord) for coord in self.perCoords]

	def _getValues(self, data, dim, coords, dataArr):
		components = [coord[self.componentIndex] for coord in coords]
		return [dataArr.sel({self.perDim: perCoord, self.componentDim: components}).data for perCoord in self.perCoords]

	def getDataKeys(self):
		return [self.dataKey]

class CoordComponentPropertyFormatted(ColumnSpec):
	def __init__(self, title, dataKey, formatStr, componentDim, propertyDim, **kwargs):
		super().__init__(**kwargs)
		self.title = title
		self.dataKey = dataKey
		self.formatStr = formatStr
		self.componentDim = componentDim
		self.propertyDim = propertyDim

	def getHeaders(self, data):
		return [self.title]

	def _getValues(self, data, dim, coords, dataArr):
		components = numpy.array(coords).T
		componentProperties = numpy.array([dataArr.sel({self.componentDim: componentList}).coords[self.propertyDim].data for componentList in components]).T
		return [[self.formatStr.format(*propertyList) for propertyList in componentProperties]]

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
	numpy.apply_along_axis(csvWriter.writerow, 0, values)

	csvFile.close()

