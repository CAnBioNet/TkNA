def matchingCoords(dataArray, dimension):
	newCoords = {}
	for coordName, coord in dataArray.coords.items():
		if coord.dims[0] == dimension:
			newCoords[coordName] = coord
	return newCoords

def withoutDim(dataArray, dimension):
	newDims = list(dataArray.dims)
	newDims.remove(dimension)

	newCoords = {}
	for coordName, coord in dataArray.coords.items():
		if coord.dims[0] != dimension:
			newCoords[coordName] = coord
	
	return newDims, newCoords

class SharedArrayParams:
	def __init__(self, shape, dtype):
		self.shape = shape
		self.dtype = dtype

