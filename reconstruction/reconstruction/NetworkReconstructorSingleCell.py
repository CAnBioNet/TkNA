import bottleneck
import itertools
from multiprocessing import Pool, shared_memory
import numpy
from scipy import stats, special
from statsmodels.stats import multitest
import xarray

from .NetworkReconstructor import NetworkReconstructor
from .util import matchingCoords, withoutDim, SharedArrayParams

numpy.seterr(all="raise")

class CorrelationWorker:
	def __init__(self, dataParams, correlationsAndPValuesParams):
		self.dataParams = dataParams
		self.correlationsAndPValuesParams = correlationsAndPValuesParams

		self.initialized = False

	def initialize(self):
		self.dataSharedMemory = shared_memory.SharedMemory(name="expressionData")
		self.data = numpy.ndarray(shape=self.dataParams.shape, dtype=self.dataParams.dtype, buffer=self.dataSharedMemory.buf)

		self.correlationsSharedMemory = shared_memory.SharedMemory(name="correlations")
		self.correlations = numpy.ndarray(shape=self.correlationsAndPValuesParams.shape, dtype=self.correlationsAndPValuesParams.dtype, buffer=self.correlationsSharedMemory.buf)

		self.pValuesSharedMemory = shared_memory.SharedMemory(name="pValues")
		self.pValues = numpy.ndarray(shape=self.correlationsAndPValuesParams.shape, dtype=self.correlationsAndPValuesParams.dtype, buffer=self.pValuesSharedMemory.buf)

	def correlationMethod(self):
		raise NotImplementedError()

	def __call__(self, indices):
		# Initialize here so that initialization only occurs in the child processes
		if not self.initialized:
			self.initialize()
			self.initialized = True

		r, p = self.correlationMethod(indices)

		self.correlations[indices] = r
		self.pValues[indices] = p

	def __del__(self):
		# No shared memory to close in the parent process
		if self.initialized:
			self.dataSharedMemory.close()
			self.correlationsSharedMemory.close()
			self.pValuesSharedMemory.close()

class CorrelationWorkerSpearman(CorrelationWorker):
	def __init__(self, rankedParams, useCoefficient, coefficient=None, **kwargs):
		super().__init__(**kwargs)
		self.rankedParams = rankedParams
		self.useCoefficient = useCoefficient
		self.coefficient = coefficient

	def initialize(self):
		super().initialize()
		self.rankedSharedMemory = shared_memory.SharedMemory(name="spearmanRanked")
		self.ranked = numpy.ndarray(shape=self.rankedParams.shape, dtype=self.rankedParams.dtype, buffer=self.rankedSharedMemory.buf)

	def correlationMethod(self, indices):
		x_ranked = self.ranked[:, indices[0]]
		y_ranked = self.ranked[:, indices[1]]
		if self.useCoefficient:
			r = 1 - (self.coefficient * numpy.sum(numpy.square(x_ranked - y_ranked)))
		else:
			covariance = ((x_ranked - x_ranked.mean(axis=-1, keepdims=True)) * (y_ranked - y_ranked.mean(axis=-1, keepdims=True))).mean(axis=-1)
			xStd = x_ranked.std(axis=-1)
			yStd = y_ranked.std(axis=-1)
			if xStd == 0 or yStd == 0:
				r = numpy.sign(covariance)
			else:
				r = covariance / (xStd * yStd)

		if abs(r) == 1:
			p = 0
		else:
			dof = len(x_ranked) - 2
			t = r * numpy.sqrt(dof / (1 - numpy.square(r)))
			p = special.stdtr(dof, -numpy.abs(t)) * 2 # Two-sided

		return r, p

	def __del__(self):
		super().__del__()
		if self.initialized:
			self.rankedSharedMemory.close()

def combineCellsByType(config, data):
	combineMethod = config["combineCellsMethod"]

	def combineCellType(cellTypeData):
		if combineMethod == "mean":
			return cellTypeData.mean(dim="cell")
		else:
			raise Exception("Unknown cell combine method")

	combinedData = data.groupby("organism").map(lambda organismData: organismData.groupby("cellType").map(combineCellType))

	# Need to re-add treatment and experiment coords b/c xarray throws them away >:(
	treatmentMap = {}
	experimentMap = {}
	for i in range(len(data.coords["organism"])):
		treatmentMap[data.coords["organism"].item(i)] = data.coords["treatment"].item(i)
		experimentMap[data.coords["organism"].item(i)] = data.coords["experiment"].item(i)
	treatmentCoord = [treatmentMap[organism] for organism in combinedData.coords["organism"].data]
	experimentCoord = [experimentMap[organism] for organism in combinedData.coords["organism"].data]
	combinedData = combinedData.assign_coords({"treatment": ("organism", treatmentCoord), "experiment": ("organism", experimentCoord)})

	return combinedData

def stackCellTypeAndMeasurable(config, combinedCells):
	return combinedCells.stack({"measurableAndCellType": ("measurable", "cellType")})

def combineDifferencePValues(config, correctedPValues):
	combineMethod = config["differenceCombinePValuesMethod"]

	def fisher(pValues):
		statistic, pValue = stats.combine_pvalues(pValues)
		return pValue

	methodMap = {
		"fisher": fisher
	}
	if combineMethod in methodMap:
		combinedPValues = xarray.apply_ufunc(methodMap[combineMethod], correctedPValues, input_core_dims=[["differentialExperiment"]], vectorize=True)
	else:
		raise Exception("Unknown p-value combine method specified")

	return combinedPValues

def filterOnDifferencePValues(config, pValues, pValueType):
	# TODO: Thresholds for measurableType as well
	threshold = config["differencePValueThresholds"][pValueType]

	if isinstance(threshold, dict):
		def filterType(typePValues):
			cellType = str(typePValues.coords["differentialCellType"].item())
			if cellType not in threshold:
				raise Exception("No threshold for {} specified in difference p-value thresholds".format(cellType))
			return typePValues <= threshold[cellType]
		filterTable = pValues.groupby("differentialCellType").map(filterType)
	else:
		filterTable = pValues <= threshold

	return filterTable

def filterOnCorrectedDifferencePValues(config, correctedPValues):
	maxPValues = correctedPValues.max(dim="differentialExperiment")
	return filterOnDifferencePValues(config, maxPValues, "corrected")

def filterOnCombinedDifferencePValues(config, combinedPValues):
	return filterOnDifferencePValues(config, combinedPValues, "combined")

def getFoldChangeSigns(config, foldChanges):
	return numpy.sign(foldChanges)

def combineAndFilterFoldChanges(config, foldChanges, foldChangeSigns):
	filterMethod = config["foldChangeFilterMethod"]

	def allSameSign(experimentFoldChanges, experimentFoldChangeSigns):
		if numpy.all(experimentFoldChangeSigns == experimentFoldChangeSigns[0]):
			return experimentFoldChangeSigns[0]
		else:
			return 0

	percentThreshold = config["foldChangeFilterPercentAgreementThreshold"]
	def percentAgreement(experimentFoldChanges, experimentFoldChangeSigns):
		signFrac = numpy.sum(experimentFoldChangeSigns) / len(experimentFoldChangeSigns)
		if 0.5 + abs(signFrac) >= percentThreshold:
			return numpy.sign(signFrac)
		else:
			return 0

	methodMap = {
		"allsamesign": allSameSign,
		"percentagreement": percentAgreement
	}
	if filterMethod in methodMap:
		combinedSigns = xarray.apply_ufunc(methodMap[filterMethod], foldChanges, foldChangeSigns, input_core_dims=[["differentialExperiment"], ["differentialExperiment"]], vectorize=True)
	else:
		raise Exception("Unknown fold change filter method specified")
	filterTable = combinedSigns != 0

	return combinedSigns, filterTable

def stackDifferenceCellTypeAndMeasurable(data):
	return data.stack({"measurableAndCellType": ("measurable", "differentialCellType")})

def calculateCorrelations(config, filteredData):
	metatreatments = config["metatreatments"]
	if metatreatments is None:
		metatreatments = {}
		networkTreatment = config["networkTreatment"]
		for experiment in filteredData.coords["experiment"].data:
			metatreatments[experiment] = [(experiment, networkTreatment)]

	metatreatmentReverseMap = {}
	for metatreatmentName, experimentAndTreatments in metatreatments.items():
		for experimentAndTreatment in experimentAndTreatments:
			metatreatmentReverseMap[experimentAndTreatment] = metatreatmentName

	metatreatmentCoords = []
	for organismCoord in filteredData.coords["organism"]:
		experimentAndTreatment = (organismCoord.experiment.item(), organismCoord.treatment.item())
		if experimentAndTreatment in metatreatmentReverseMap:
			metatreatmentCoords.append(metatreatmentReverseMap[experimentAndTreatment])
		else:
			metatreatmentCoords.append(None)

	filteredData = filteredData.assign_coords(coords={"metatreatment": ("organism", metatreatmentCoords)})

	def prepareSpearman(treatmentData):
		spearmanKwargs = {}

		# See https://bottleneck.readthedocs.io/en/v1.3.4/reference.html#bottleneck.rankdata
		rankedParams = SharedArrayParams(treatmentData.shape, numpy.float64)
		rankedNBytes = numpy.prod(rankedParams.shape) * rankedParams.dtype().itemsize
		rankedSharedMemory = shared_memory.SharedMemory(create=True, size=rankedNBytes, name="spearmanRanked")
		ranked = numpy.ndarray(rankedParams.shape, rankedParams.dtype, buffer=rankedSharedMemory.buf)
		# Ranks begin at 1
		ranked[:] = bottleneck.rankdata(treatmentData, axis=0)

		spearmanKwargs["rankedParams"] = rankedParams

		n = treatmentData.sizes["organism"]

		# Check that all ranks are distinct before using formula
		modes = stats.mode(ranked, axis=1)
		spearmanKwargs["useCoefficient"] = numpy.all(modes.count == 1)
		if spearmanKwargs["useCoefficient"]:
			# From https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient#Definition_and_calculation
			spearmanKwargs["coefficient"] = 6.0 / (n * (numpy.square(n) - 1))

		spearmanEndArgs = [rankedSharedMemory]
		return spearmanKwargs, spearmanEndArgs

	def endSpearman(rankedSharedMemory):
		rankedSharedMemory.close()
		rankedSharedMemory.unlink()

	correlationMethod = config["correlationMethod"]
	methodMap = {
		"spearman": (CorrelationWorkerSpearman, prepareSpearman, endSpearman)
	}
	if correlationMethod not in methodMap:
		raise Exception("Unknown correlation method specified")

	treatmentData = filteredData

	def calculateForMetatreatment(treatmentData):
		# TODO: Create correlation/pvalue dtype?

		correlationsAndPValuesParams = SharedArrayParams((treatmentData.sizes["measurableAndCellType"], treatmentData.sizes["measurableAndCellType"]), numpy.float64)
		correlationsAndPValuesNBytes = numpy.prod(correlationsAndPValuesParams.shape) * correlationsAndPValuesParams.dtype().itemsize
		correlationsSharedMemory = shared_memory.SharedMemory(create=True, size=correlationsAndPValuesNBytes, name="correlations")
		correlations = numpy.ndarray(correlationsAndPValuesParams.shape, correlationsAndPValuesParams.dtype, buffer=correlationsSharedMemory.buf)
		pValuesSharedMemory = shared_memory.SharedMemory(create=True, size=correlationsAndPValuesNBytes, name="pValues")
		pValues = numpy.ndarray(correlationsAndPValuesParams.shape, correlationsAndPValuesParams.dtype, buffer=pValuesSharedMemory.buf)

		treatmentDataParams = SharedArrayParams(treatmentData.shape, treatmentData.dtype)
		treatmentDataSharedMemory = shared_memory.SharedMemory(create=True, size=treatmentData.nbytes, name="expressionData")
		dataCopy = numpy.ndarray(treatmentDataParams.shape, dtype=treatmentDataParams.dtype, buffer=treatmentDataSharedMemory.buf)
		dataCopy[:] = treatmentData.values[:]

		(methodWorker, methodPrepare, methodEnd) = methodMap[correlationMethod]
		if methodPrepare:
			methodKwargs, endArgs = methodPrepare(treatmentData)
		else:
			methodKwargs = {}
			endArgs = []

		worker = methodWorker(**methodKwargs, dataParams=treatmentDataParams, correlationsAndPValuesParams=correlationsAndPValuesParams)
		with Pool() as p:
			p.map(worker, itertools.combinations(range(treatmentData.sizes["measurableAndCellType"]), 2))

		treatmentDataSharedMemory.close()
		treatmentDataSharedMemory.unlink()

		diagonal = numpy.diag(numpy.ones(treatmentData.sizes["measurableAndCellType"]))
		correlations += correlations.transpose() + diagonal
		pValues += pValues.transpose() + diagonal

		# As far as I am able to discern, this is the easiest way to rename these levels
		renamed1 = treatmentData.unstack(dim="measurableAndCellType").rename({"measurable": "measurable1", "cellType": "cellType1"}).stack({"measurableAndCellType1": ("measurable1", "cellType1")})
		renamed2 = treatmentData.unstack(dim="measurableAndCellType").rename({"measurable": "measurable2", "cellType": "cellType2"}).stack({"measurableAndCellType2": ("measurable2", "cellType2")})

		renamed1 = renamed1.dropna("measurableAndCellType1")
		renamed2 = renamed2.dropna("measurableAndCellType2")

		results = xarray.Dataset(
			data_vars={"correlations": (("measurableAndCellType1", "measurableAndCellType2"), correlations.copy()), "pValues": (("measurableAndCellType1", "measurableAndCellType2"), pValues.copy())},
			coords={"measurableAndCellType1": renamed1.coords["measurableAndCellType1"], "measurableAndCellType2": renamed2.coords["measurableAndCellType2"], "measurableType1": ("measurableAndCellType1", treatmentData.coords["measurableType"].data), "measurableType2": ("measurableAndCellType2", treatmentData.coords["measurableType"].data)}
		)

		correlationsSharedMemory.close()
		correlationsSharedMemory.unlink()
		pValuesSharedMemory.close()
		pValuesSharedMemory.unlink()

		if methodEnd:
			methodEnd(*endArgs)

		return results

	correlationResults = treatmentData.groupby("metatreatment").map(calculateForMetatreatment)

	return correlationResults["correlations"], correlationResults["pValues"]

def combineCorrelationPValues(config, pValues):
	combineMethod = config["correlationCombinePValuesMethod"]

	metatreatmentsToCombine = config["metatreatmentsToCombine"]
	if metatreatmentsToCombine is None:
		metatreatmentsToCombine = [metatreatment for metatreatment in set(pValues.coords["metatreatment"].data) if metatreatment is not None]

	# Skip combining if there is only one metatreatment selected
	if len(metatreatmentsToCombine) == 1:
		return pValues.sel(metatreatment=metatreatmentsToCombine[0])

	selectedPValues = pValues.sel(metatreatment=metatreatmentsToCombine)

	def fisher():
		# Ignore errors that occur for p = 0
		with numpy.errstate(divide="ignore"):
			chi2 = -2 * (numpy.log(pValues)).sum(dim="metatreatment")
		dof = 2 * pValues.sizes["metatreatment"]
		combinedPValues = stats.chi2.sf(chi2, df=dof)

		newDims, newCoords = withoutDim(pValues, "metatreatment")
		combinedPValues = xarray.DataArray(combinedPValues, dims=newDims, coords=newCoords)
		return combinedPValues

	methodMap = {
		"fisher": fisher
	}
	if combineMethod in methodMap:
		combinedPValues = methodMap[combineMethod]()
	else:
		raise Exception("Unknown p-value combine method specified")

	return combinedPValues

def correctCorrelationPValues(config, combinedPValues):
	correctMethod = config["correlationCorrectPValuesMethod"]

	def fdr(pValues):
		nonNanIndices = numpy.argwhere(~numpy.isnan(pValues))
		nonNanPValues = pValues[nonNanIndices].flatten()
		rejected, correctedPValues = multitest.fdrcorrection(nonNanPValues)
		correctedPValuesWithNans = pValues
		numpy.put(correctedPValuesWithNans, nonNanIndices, correctedPValues)
		return correctedPValuesWithNans

	methodMap = {
		"fdr": fdr
	}
	if not correctMethod in methodMap:
		raise Exception("Unknown p-value correct method specified")

	edgeTypesCorrected = []
	def correctPValues(pValueMatrix):
		# If this edge type was already done, return zeros instead of recalculating
		firstValue = pValueMatrix[0, 0]
		edgeType = frozenset({(firstValue.measurableType1.item(), firstValue.measurableAndCellType1.item()[1]), (firstValue.measurableType2.item(), firstValue.measurableAndCellType2.item()[1])})
		if edgeType in edgeTypesCorrected:
			return numpy.zeros(pValueMatrix.shape)
		edgeTypesCorrected.append(edgeType)

		# If these correlations are within a type, they're on the diagonal,
		# so only correct the ones above the diagonal and return zeros everywhere else
		if len(edgeType) == 1:
			numMeasurables = pValueMatrix.sizes["measurableAndCellType1"]
			upperTriangleIndices = numpy.triu_indices(numMeasurables, k=1)
			upperTriangle = pValueMatrix.to_numpy()[upperTriangleIndices]
			correctedPValues = methodMap[correctMethod](upperTriangle)
			correctedPValueMatrix = numpy.zeros((numMeasurables, numMeasurables))
			correctedPValueMatrix[upperTriangleIndices] = correctedPValues
		# Otherwise, correct the entire set of edges
		else:
			correctedPValues = methodMap[correctMethod](pValueMatrix.to_numpy().flatten(order="C"))
			correctedPValueMatrix = numpy.reshape(correctedPValues, pValueMatrix.shape, order="C")
		return correctedPValueMatrix

	def correctTypeProduct(pValues):
		return pValues.groupby("cellType1").map(lambda type1Data: type1Data.groupby("cellType2").map(correctPValues))

	# Perform corrections
	correctedPValues = combinedPValues.groupby("measurableType1").map(lambda measurableType1Data: measurableType1Data.groupby("measurableType2").map(correctTypeProduct))

	# Copy data across the diagonal
	# (this works because the correction method makes all duplicate/diagonal entries 0)
	correctedPValues += correctedPValues.to_numpy().T

	return correctedPValues

def filterDiagonals(pValues):
	diagonalFilter = numpy.ones([pValues.sizes["measurableAndCellType1"], pValues.sizes["measurableAndCellType2"]], dtype=bool)
	numpy.fill_diagonal(diagonalFilter, 0)
	diagonalFilter = xarray.DataArray(diagonalFilter, dims=pValues.dims, coords=pValues.coords)
	return diagonalFilter

def filterOnCorrelationPValues(config, pValues, pValueType):
	# TODO: Thresholds for measurableType combinations as well
	threshold = config["correlationPValueThresholds"][pValueType]

	if isinstance(threshold, dict):
		def filterTypeCombo(typeComboPValues):
			cellTypes = [str(typeComboPValues.coords[typeCoord][0].values) for typeCoord in ["cellType1", "cellType2"]]
			cellTypeString = "({}, {})".format(*cellTypes)
			cellTypeStringReversed = "({1}, {0})".format(*cellTypes)
			if cellTypeString in threshold:
				typeComboThreshold = threshold[cellTypeString]
			elif cellTypeStringReversed in threshold:
				typeComboThreshold = threshold[cellTypeStringReversed]
			else:
				raise Exception("No threshold for {} specified in correlation p-value thresholds".format(cellTypeString))
			return typeComboPValues <= typeComboThreshold
		filterTable = pValues.groupby("cellType1").map(lambda cellType1Data: cellType1Data.groupby("cellType2").map(filterTypeCombo))
	else:
		filterTable = pValues <= threshold

	return filterTable

def filterOnIndividualCorrelationPValues(config, pValues):
	maxPValues = pValues.max(dim="metatreatment")
	return filterOnCorrelationPValues(config, maxPValues, "individual")

def filterOnCombinedCorrelationPValues(config, combinedPValues):
	return filterOnCorrelationPValues(config, combinedPValues, "combined")

def filterOnCorrectedCorrelationPValues(config, correctedPValues):
	return filterOnCorrelationPValues(config, correctedPValues, "corrected")

def combineAndFilterCorrelations(config, correlations, correlationSigns):
	filterMethod = config["correlationFilterMethod"]

	def allSameSign():
		filterTable = (correlationSigns.isel(metatreatment=0) == correlationSigns).all(dim="metatreatment")
		combinedSigns = correlationSigns.isel(metatreatment=0).where(filterTable, 0)
		return combinedSigns, filterTable

	# Should be greater than 0.5, as otherwise there can be ties between signs
	percentThreshold = config["correlationFilterPercentAgreementThreshold"]
	def percentAgreement():
		numMetatreatments = correlationSigns.sizes["metatreatment"]
		metatreatmentAxis = correlationSigns.get_axis_num("metatreatment")
		positiveFilterTable = (numpy.count_nonzero(correlationSigns == 1, axis=metatreatmentAxis) / numMetatreatments) > percentThreshold
		negativeFilterTable = (numpy.count_nonzero(correlationSigns == -1, axis=metatreatmentAxis) / numMetatreatments) > percentThreshold
		filterTable = positiveFilterTable | negativeFilterTable
		ones = xarray.ones_like(correlationSigns.isel(metatreatment=0))
		combinedSigns = ones.where(positiveFilterTable, 0) - ones.where(negativeFilterTable, 0)
		return combinedSigns, filterTable

	methodMap = {
		"allsamesign": allSameSign,
		"percentagreement": percentAgreement
	}
	if filterMethod in methodMap:
		combinedSigns, filterTable = methodMap[filterMethod]()
	else:
		raise Exception("Unknown correlation filter method specified")

	return combinedSigns, filterTable

def filterOnExpectedEdges(config, foldChangeSigns, correlationSigns, measurableFilter):
	stacked1 = foldChangeSigns.stack({"measurableAndCellType1": ("measurable", "differentialCellType")})
	stacked1 = stacked1.sel(measurableAndCellType1=stacked1.coords["measurableAndCellType1"][measurableFilter.rename(measurableAndCellType="measurableAndCellType1")])
	stacked2 = foldChangeSigns.stack({"measurableAndCellType2": ("measurable", "differentialCellType")})
	stacked2 = stacked2.sel(measurableAndCellType2=stacked2.coords["measurableAndCellType2"][measurableFilter.rename(measurableAndCellType="measurableAndCellType2")])
	foldChangeSignProducts = stacked1 * stacked2
	pucCompliant = foldChangeSignProducts == correlationSigns
	return pucCompliant, foldChangeSignProducts

class NetworkReconstructorSingleCell(NetworkReconstructor):
	def reconstructNetwork(self, config, data, **kwargs):
		def stageCombineCellsByType(allData):
			allData["cellsCombined"] = combineCellsByType(config, allData["cellData"])
			allData["stacked"] = stackCellTypeAndMeasurable(config, allData["cellsCombined"])

		def stageCombineDifferencePValues(allData):
			# NB: Combined p-values are generated from the corrected p-values, rather than vice versa
			allData["combinedDifferencePValues"] = combineDifferencePValues(config, allData["correctedDifferencePValues"])

		def stageFilterOnDifferences(allData):
			allData["correctedDifferencePValueFilter"] = filterOnCorrectedDifferencePValues(config, allData["correctedDifferencePValues"])
			allData["combinedDifferencePValueFilter"] = filterOnCombinedDifferencePValues(config, allData["combinedDifferencePValues"])
			allData["foldChangeSigns"] = getFoldChangeSigns(config, allData["foldChanges"])
			allData["combinedFoldChangeSigns"], allData["foldChangeFilter"] = combineAndFilterFoldChanges(config, allData["foldChanges"], allData["foldChangeSigns"])
			allData["measurableFilter"] = allData["correctedDifferencePValueFilter"] & allData["combinedDifferencePValueFilter"] & allData["foldChangeFilter"]

			allData["measurableFilterStacked"] = stackDifferenceCellTypeAndMeasurable(allData["measurableFilter"])
			allData["filteredData"] = allData["stacked"].sel(measurableAndCellType=allData["stacked"].coords["measurableAndCellType"][allData["measurableFilterStacked"]])

		def stageComputeCorrelations(allData):
			allData["correlationCoefficients"], allData["correlationPValues"] = calculateCorrelations(config, allData["filteredData"])
			allData["combinedCorrelationPValues"] = combineCorrelationPValues(config, allData["correlationPValues"])
			allData["correlationSigns"] = numpy.sign(allData["correlationCoefficients"])

			if not config["correctCorrelationPValuesAfterConsistencyFiltering"]:
				allData["correctedCorrelationPValues"] = correctCorrelationPValues(config, allData["combinedCorrelationPValues"])

		def stageFilterOnCorrelations(allData):
			allData["combinedCorrelationSigns"], allData["correlationFilter"] = combineAndFilterCorrelations(config, allData["correlationCoefficients"], allData["correlationSigns"])

			if config["correctCorrelationPValuesAfterConsistencyFiltering"]:
				allData["correctedCorrelationPValues"] = correctCorrelationPValues(config, allData["combinedCorrelationPValues"].where(allData["correlationFilter"], other=numpy.NaN))

			allData["diagonalFilter"] = filterDiagonals(allData["correctedCorrelationPValues"])
			allData["individualCorrelationPValueFilter"] = filterOnIndividualCorrelationPValues(config, allData["correlationPValues"])
			allData["combinedCorrelationPValueFilter"] = filterOnCombinedCorrelationPValues(config, allData["combinedCorrelationPValues"])
			allData["correctedCorrelationPValueFilter"] = filterOnCorrectedCorrelationPValues(config, allData["correctedCorrelationPValues"])
			allData["edgeFilter"] = allData["diagonalFilter"] & allData["individualCorrelationPValueFilter"] & allData["combinedCorrelationPValueFilter"] & allData["correctedCorrelationPValueFilter"] & allData["correlationFilter"]

		def stageFilterToExpectedEdges(allData):
			allData["expectedEdgeFilter"], allData["foldChangeSignProducts"] = filterOnExpectedEdges(config, allData["combinedFoldChangeSigns"], allData["combinedCorrelationSigns"], allData["measurableFilterStacked"])
			allData["edgeFilter"] &= allData["expectedEdgeFilter"]

		def stageCreateEdgeList(allData):
			allData["edges"] =  allData["combinedCorrelationSigns"].where(allData["edgeFilter"], other=0)

		allData = {}
		allData["cellData"] = data["cellData"]
		allData["correctedDifferencePValues"] = data["pValue"]
		allData["foldChanges"] = data["foldChange"]
		stages = [stageCombineCellsByType, stageCombineDifferencePValues, stageFilterOnDifferences, stageComputeCorrelations, stageFilterOnCorrelations, stageFilterToExpectedEdges, stageCreateEdgeList]
		return self.runPipeline(stages, allData, **kwargs)

