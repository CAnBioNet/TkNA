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
	def __init__(self, treatmentDataParams, correlationsAndPValuesParams):
		self.treatmentDataParams = treatmentDataParams
		self.correlationsAndPValuesParams = correlationsAndPValuesParams

		self.initialized = False

	def initialize(self):
		self.treatmentDataSharedMemory = shared_memory.SharedMemory(name="treatmentData")
		self.treatmentData = numpy.ndarray(shape=self.treatmentDataParams.shape, dtype=self.treatmentDataParams.dtype, buffer=self.treatmentDataSharedMemory.buf)

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
			self.treatmentDataSharedMemory.close()
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
		x_ranked = self.ranked[indices[0], :]
		y_ranked = self.ranked[indices[1], :]

		nans = numpy.isnan(x_ranked) | numpy.isnan(y_ranked)
		numNans = numpy.sum(nans)
		if numNans != 0:
			notNans = ~nans
			x = self.treatmentData[indices[0], :][notNans]
			y = self.treatmentData[indices[1], :][notNans]

			if x.size <= 1: # y has the same size as x, so only need to check x
				# TODO: Warning
				return numpy.nan, numpy.nan

			x_ranked = bottleneck.rankdata(x)
			y_ranked = bottleneck.rankdata(y)

			# TODO: Precompute size -> coefficient table instead?
			self.useCoefficient = False

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

		rEps = numpy.finfo(type(r)).eps
		if abs(abs(r) - 1) <= rEps:
			p = 0
			r = numpy.sign(r) # Round r to 1 or -1
		else:
			dof = len(x_ranked) - 2
			t = r * numpy.sqrt(dof / (1 - numpy.square(r)))
			p = special.stdtr(dof, -numpy.abs(t)) * 2 # Two-sided

		return r, p

	def __del__(self):
		super().__del__()
		if self.initialized:
			self.rankedSharedMemory.close()

class CorrelationWorkerPearson(CorrelationWorker):
	def correlationMethod(self, indices):
		x = self.treatmentData[indices[0], :]
		y = self.treatmentData[indices[1], :]
		r, p = stats.pearsonr(x, y)
		return r, p

def computeDifferencePValues(config, data):
	treatments = config["comparisonTreatments"]
	if not len(treatments) == 2:
		raise Exception("More than 2 treatments given to compare")

	differenceMethod = config["differenceMethod"]

	def mannWhitneyU(experimentData):
		grouped = experimentData.groupby("treatment")
		treatmentData = [grouped[treatment] for treatment in treatments]
		statistics, pValues = stats.mannwhitneyu(treatmentData[0], treatmentData[1], axis=1, nan_policy="omit")
		pValues = xarray.DataArray(pValues, dims=["measurable"], coords=matchingCoords(data, "measurable"))
		return pValues

	def independentTTest(experimentData):
		grouped = experimentData.groupby("treatment")
		treatmentData = [grouped[treatment] for treatment in treatments]

		# This `with` statement is a workaround for two separate bugs.
		with numpy.errstate(
			# Workaround for a numpy bug. Details follow:
			# ttest_ind, if using nan_policy="omit", uses masked-array division from numpy, which tests that the result isn't (effectively) infinity.
			# However, as part of this check, it multiplies the dividend (in this case, the sum of the squared deviations of the measurements in a treatment)
			# by the smallest possible float value, which causes an underflow error if it's <1 (i.e. for measurables with low variance).
			# Therefore we ignore underflow errors for this calculation.
			# See https://github.com/numpy/numpy/issues/4895 for this particular issue.
			# See https://github.com/numpy/numpy/issues/22347 for a related issue and more general information on this part of numpy.
			under="ignore",
			# Workaround for a scipy bug. Details follow:
			# In ttest_ind, for Student's t-test, scipy takes the reciprocals of the sample sizes to compute the t statistic.
			# However, these sample sizes exclude any nans in the data (since we're using nan_policy="omit"),
			# so if a treatment has no data for a particular measurable a division by 0 error will occur.
			# See https://github.com/scipy/scipy/issues/14651 for a mention of this issue.
			divide="ignore"
		):
			statistics, pValues = stats.ttest_ind(treatmentData[0], treatmentData[1], axis=1, nan_policy="omit")

		pValues = xarray.DataArray(pValues, dims=["measurable"], coords=matchingCoords(data, "measurable"))
		return pValues

	methodMap = {
		"mannwhitney": mannWhitneyU,
		"independentttest": independentTTest
	}
	if differenceMethod in methodMap:
		pValues = data.groupby("experiment").map(methodMap[differenceMethod])
	else:
		raise Exception("Unknown difference method specified")

	return pValues

def combineDifferencePValues(config, pValues):
	combineMethod = config["differenceCombinePValuesMethod"]

	def fisher(measurablePValues):
		statistic, pValue = stats.combine_pvalues(measurablePValues[~numpy.isnan(measurablePValues)])
		return pValue

	methodMap = {
		"fisher": fisher
	}
	if combineMethod in methodMap:
		combinedPValues = xarray.apply_ufunc(methodMap[combineMethod], pValues, input_core_dims=[["experiment"]], vectorize=True)
	else:
		raise Exception("Unknown p-value combine method specified")

	return combinedPValues

def correctDifferencePValues(config, combinedPValues):
	correctMethod = config["differenceCorrectPValuesMethod"]

	def fdr(typePValues):
		rejected, correctedPValues = multitest.fdrcorrection(typePValues)
		return correctedPValues

	methodMap = {
		"fdr": fdr
	}
	if correctMethod in methodMap:
		correctedPValues = combinedPValues.groupby("measurableType").map(methodMap[correctMethod])
	else:
		raise Exception("Unknown p-value combine method specified")

	return correctedPValues

def filterOnDifferencePValues(config, pValues, pValueType):
	threshold = config["differencePValueThresholds"][pValueType]

	if isinstance(threshold, dict):
		def filterType(typePValues):
			measurableType = str(typePValues.coords["measurableType"][0].values)
			if measurableType not in threshold:
				raise Exception("No threshold for {} specified in difference p-value thresholds".format(measurableType))
			return typePValues <= threshold[measurableType]
		filterTable = pValues.groupby("measurableType").map(filterType)
	else:
		filterTable = pValues <= threshold

	return filterTable

def filterOnIndividualDifferencePValues(config, pValues):
	maxPValues = pValues.max(dim="experiment")
	return filterOnDifferencePValues(config, maxPValues, "individual")

def filterOnCombinedDifferencePValues(config, combinedPValues):
	return filterOnDifferencePValues(config, combinedPValues, "combined")

def filterOnCorrectedDifferencePValues(config, correctedPValues):
	return filterOnDifferencePValues(config, correctedPValues, "corrected")

def computeFoldChanges(config, data):
	treatments = config["comparisonTreatments"]
	if not len(treatments) == 2:
		raise Exception("More than 2 treatments given to compare")

	foldChangeType = config["foldChangeType"]
	if not foldChangeType in ["median", "mean"]:
		raise Exception("Invalid fold change type")

	def computeFoldChanges_(experimentData):
		grouped = experimentData.groupby("treatment")
		treatmentData = [grouped[treatment] for treatment in treatments]
		combinedOrganismData = [getattr(treatmentDataset, foldChangeType)("organism") for treatmentDataset in treatmentData]
		return numpy.log2(combinedOrganismData[0] / combinedOrganismData[1])

	foldChanges = data.groupby("experiment").map(computeFoldChanges_)
	foldChangeSigns = numpy.sign(foldChanges)
	return foldChanges, foldChangeSigns

def combineAndFilterFoldChanges(config, foldChanges, foldChangeSigns):
	filterMethod = config["foldChangeFilterMethod"]

	def allSameSign(experimentFoldChanges, experimentFoldChangeSigns):
		if numpy.all((experimentFoldChangeSigns == experimentFoldChangeSigns[0]) & ~numpy.isnan(experimentFoldChangeSigns)):
			return experimentFoldChangeSigns[0]
		else:
			return 0

	percentThreshold = config["foldChangeFilterPercentAgreementThreshold"]
	def percentAgreement(experimentFoldChanges, experimentFoldChangeSigns):
		percentAgreementDecimals = 6

		numExperiments = len(experimentFoldChangeSigns)

		positiveOverThreshold = numpy.around(numpy.count_nonzero(correlationSigns == 1) / numExperiments, decimals=percentAgreementDecimals) >= percentThreshold
		if positiveOverThreshold:
			return 1

		negativeOverThreshold = numpy.around(numpy.count_nonzero(correlationSigns == -1) / numExperiments, decimals=percentAgreementDecimals) >= percentThreshold
		if negativeOverThreshold:
			return -1

		return 0

	methodMap = {
		"allsamesign": allSameSign,
		"percentagreement": percentAgreement
	}
	if filterMethod in methodMap:
		combinedSigns = xarray.apply_ufunc(methodMap[filterMethod], foldChanges, foldChangeSigns, input_core_dims=[["experiment"], ["experiment"]], vectorize=True)
	else:
		raise Exception("Unknown fold change filter method specified")
	filterTable = combinedSigns != 0

	return combinedSigns, filterTable

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
		ranked[:] = bottleneck.nanrankdata(treatmentData, axis=1)

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
		"spearman": (CorrelationWorkerSpearman, prepareSpearman, endSpearman),
		"pearson": (CorrelationWorkerPearson, None, None)
	}
	if correlationMethod not in methodMap:
		raise Exception("Unknown correlation method specified")

	treatmentData = filteredData

	def calculateForMetatreatment(treatmentData):
		# TODO: Create correlation/pvalue dtype?

		correlationsAndPValuesParams = SharedArrayParams((treatmentData.sizes["measurable"], treatmentData.sizes["measurable"]), numpy.float64)
		correlationsAndPValuesNBytes = numpy.prod(correlationsAndPValuesParams.shape) * correlationsAndPValuesParams.dtype().itemsize
		correlationsSharedMemory = shared_memory.SharedMemory(create=True, size=correlationsAndPValuesNBytes, name="correlations")
		correlations = numpy.ndarray(correlationsAndPValuesParams.shape, correlationsAndPValuesParams.dtype, buffer=correlationsSharedMemory.buf)
		pValuesSharedMemory = shared_memory.SharedMemory(create=True, size=correlationsAndPValuesNBytes, name="pValues")
		pValues = numpy.ndarray(correlationsAndPValuesParams.shape, correlationsAndPValuesParams.dtype, buffer=pValuesSharedMemory.buf)

		treatmentDataParams = SharedArrayParams(treatmentData.shape, treatmentData.dtype)
		treatmentDataSharedMemory = shared_memory.SharedMemory(create=True, size=treatmentData.nbytes, name="treatmentData")
		dataCopy = numpy.ndarray(treatmentDataParams.shape, dtype=treatmentDataParams.dtype, buffer=treatmentDataSharedMemory.buf)
		dataCopy[:] = treatmentData.values[:]

		(methodWorker, methodPrepare, methodEnd) = methodMap[correlationMethod]
		if methodPrepare:
			methodKwargs, endArgs = methodPrepare(treatmentData)
		else:
			methodKwargs = {}
			endArgs = []

		worker = methodWorker(**methodKwargs, treatmentDataParams=treatmentDataParams, correlationsAndPValuesParams=correlationsAndPValuesParams)
		with Pool() as p:
			p.map(worker, itertools.combinations(range(treatmentData.sizes["measurable"]), 2))

		treatmentDataSharedMemory.close()
		treatmentDataSharedMemory.unlink()

		diagonal = numpy.diag(numpy.ones(treatmentData.sizes["measurable"]))
		correlations += correlations.transpose() + diagonal
		pValues += pValues.transpose() + diagonal

		results = xarray.Dataset(
			data_vars={"correlations": (("measurable1", "measurable2"), correlations.copy()), "pValues": (("measurable1", "measurable2"), pValues.copy())},
			coords={"measurable1": treatmentData.coords["measurable"].data, "measurable2": treatmentData.coords["measurable"].data, "measurableType1": ("measurable1", treatmentData.coords["measurableType"].data), "measurableType2": ("measurable2", treatmentData.coords["measurableType"].data)}
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
		# Ignore errors that occur for p = 0 / NaN
		with numpy.errstate(divide="ignore", invalid="ignore"):
			chi2 = -2 * (numpy.log(selectedPValues)).sum(dim="metatreatment", skipna=True)
		dof = 2 * pValues.count(dim="metatreatment")
		combinedPValues = stats.chi2.sf(chi2, df=dof)
		combinedPValues = xarray.DataArray(combinedPValues, dims=["measurable1", "measurable2"], coords={"measurable1": pValues.coords["measurable1"].data, "measurable2": pValues.coords["measurable2"].data, "measurableType1": ("measurable1", pValues.coords["measurableType1"].data), "measurableType2": ("measurable2", pValues.coords["measurableType2"].data)})
		return combinedPValues

	methodMap = {
		"fisher": fisher
	}
	if combineMethod in methodMap:
		combinedPValues = methodMap[combineMethod]()
	else:
		raise Exception("Unknown p-value combine method specified")

	combinedPValues = combinedPValues.where(pValues.notnull().any(dim="metatreatment"))

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
		edgeType = frozenset({firstValue.measurableType1.item(), firstValue.measurableType2.item()})
		if edgeType in edgeTypesCorrected:
			return numpy.zeros(pValueMatrix.shape)
		edgeTypesCorrected.append(edgeType)

		# If these correlations are within a type, they're on the diagonal,
		# so only correct the ones above the diagonal and return zeros everywhere else
		if len(edgeType) == 1:
			numMeasurables = pValueMatrix.sizes["measurable1"]
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

	# Perform corrections
	correctedPValues = combinedPValues.groupby("measurableType1").map(lambda measurableType1Data: measurableType1Data.groupby("measurableType2").map(correctPValues))

	# Copy data across the diagonal
	# (this works because the correction method makes all duplicate/diagonal entries 0)
	correctedPValues += correctedPValues.to_numpy().T

	return correctedPValues

def filterDiagonals(pValues):
	diagonalFilter = numpy.ones([pValues.sizes["measurable1"], pValues.sizes["measurable2"]], dtype=bool)
	numpy.fill_diagonal(diagonalFilter, 0)
	diagonalFilter = xarray.DataArray(diagonalFilter, dims=pValues.dims, coords=pValues.coords)
	return diagonalFilter

def filterOnCorrelationPValues(config, pValues, pValueType):
	threshold = config["correlationPValueThresholds"][pValueType]

	if isinstance(threshold, dict):
		def filterTypeCombo(typeComboPValues):
			measurableTypes = [str(typeComboPValues.coords[typeCoord][0].values) for typeCoord in ["measurableType1", "measurableType2"]]
			measurableTypeString = "({}, {})".format(*measurableTypes)
			measurableTypeStringReversed = "({1}, {0})".format(*measurableTypes)
			if measurableTypeString in threshold:
				typeComboThreshold = threshold[measurableTypeString]
			elif measurableTypeStringReversed in threshold:
				typeComboThreshold = threshold[measurableTypeStringReversed]
			else:
				raise Exception("No threshold for {} specified in correlation p-value thresholds".format(measurableTypeString))
			return typeComboPValues <= typeComboThreshold
		filterTable = pValues.groupby("measurableType1").map(lambda measurableType1Data: measurableType1Data.groupby("measurableType2").map(filterTypeCombo))
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

	metatreatmentsToCombine = config["metatreatmentsForDirectionFiltering"]
	if metatreatmentsToCombine is None:
		metatreatmentsToCombine = [metatreatment for metatreatment in set(correlations.coords["metatreatment"].data) if metatreatment is not None]

	# Skip combining if there is only one metatreatment selected
	if len(metatreatmentsToCombine) == 1:
		return pValues.sel(metatreatment=metatreatmentsToCombine[0])

	selectedCorrelations = correlations.sel(metatreatment=metatreatmentsToCombine)

	def allSameSign():
		filterTable = ((correlationSigns == correlationSigns.isel(metatreatment=0)) & ~numpy.isnan(correlationSigns)).all(dim="metatreatment")
		passingMetatreatments = filterTable.expand_dims(dim={"metatreatment": correlationSigns.coords["metatreatment"]})
		combinedSigns = correlationSigns.isel(metatreatment=0).where(filterTable, 0)
		return combinedSigns, filterTable, passingMetatreatments

	# Should be greater than 0.5, as otherwise there can be ties between signs
	percentThreshold = config["correlationFilterPercentAgreementThreshold"]
	def percentAgreement():
		numMetatreatments = correlationSigns.sizes["metatreatment"]
		metatreatmentAxis = correlationSigns.get_axis_num("metatreatment")
		filterTableDims, filterTableCoords = withoutDim(correlationSigns, "metatreatment")
		percentAgreementDecimals = 6
		positiveFilterTable = xarray.DataArray(numpy.around(numpy.count_nonzero(correlationSigns == 1, axis=metatreatmentAxis) / numMetatreatments, decimals=percentAgreementDecimals) >= percentThreshold, dims=filterTableDims, coords=filterTableCoords)
		negativeFilterTable = xarray.DataArray(numpy.around(numpy.count_nonzero(correlationSigns == -1, axis=metatreatmentAxis) / numMetatreatments, decimals=percentAgreementDecimals) >= percentThreshold, dims=filterTableDims, coords=filterTableCoords)
		filterTable = positiveFilterTable | negativeFilterTable
		ones = xarray.ones_like(correlationSigns.isel(metatreatment=0))
		combinedSigns = ones.where(positiveFilterTable, 0) - ones.where(negativeFilterTable, 0)
		passingMetatreatments = xarray.apply_ufunc(lambda metatreatmentData: metatreatmentData == combinedSigns, correlationSigns, input_core_dims=[filterTableDims], output_core_dims=[filterTableDims], vectorize=True)
		return combinedSigns, filterTable, passingMetatreatments

	methodMap = {
		"allsamesign": allSameSign,
		"percentagreement": percentAgreement
	}
	if filterMethod in methodMap:
		combinedSigns, filterTable, passingMetatreatments = methodMap[filterMethod]()
	else:
		raise Exception("Unknown correlation filter method specified")

	return combinedSigns, filterTable, passingMetatreatments

def filterOnExpectedEdges(config, foldChangeSigns, correlationSigns):
	foldChangeSignProducts = foldChangeSigns.rename({"measurable": "measurable1"}) * foldChangeSigns.rename({"measurable": "measurable2"})
	pucCompliant = foldChangeSignProducts == correlationSigns
	return pucCompliant, foldChangeSignProducts

class NetworkReconstructorAggregate(NetworkReconstructor):
	def reconstructNetwork(self, config, data, **kwargs):
		skip = False

		def computeDifferences(allData):
			allData["differencePValues"] = computeDifferencePValues(config, allData["originalData"])
			allData["combinedDifferencePValues"] = combineDifferencePValues(config, allData["differencePValues"])
			allData["correctedDifferencePValues"] = correctDifferencePValues(config, allData["combinedDifferencePValues"])

			allData["foldChanges"], allData["foldChangeSigns"] = computeFoldChanges(config, allData["originalData"])

		def filterOnDifferences(allData):
			allData["individualDifferencePValueFilter"] = filterOnIndividualDifferencePValues(config, allData["differencePValues"])
			allData["combinedDifferencePValueFilter"] = filterOnCombinedDifferencePValues(config, allData["combinedDifferencePValues"])
			allData["correctedDifferencePValueFilter"] = filterOnCorrectedDifferencePValues(config, allData["correctedDifferencePValues"])
			allData["combinedFoldChangeSigns"], allData["foldChangeFilter"] = combineAndFilterFoldChanges(config, allData["foldChanges"], allData["foldChangeSigns"])
			allData["measurableFilter"] = allData["individualDifferencePValueFilter"] & allData["combinedDifferencePValueFilter"] & allData["correctedDifferencePValueFilter"] & allData["foldChangeFilter"]

			allData["filteredData"] = allData["originalData"].sel(measurable=allData["originalData"].coords["measurable"][allData["measurableFilter"]])

		def computeCorrelations(allData):
			nonlocal skip
			if allData["filteredData"].sizes["measurable"] == 0:
				print("WARNING: 0 measurables passed differential expression filters. No edge list generated.")
				skip = True
				return

			allData["correlationCoefficients"], allData["correlationPValues"] = calculateCorrelations(config, allData["filteredData"])

			allData["correlationSigns"] = numpy.sign(allData["correlationCoefficients"])
			allData["combinedCorrelationSigns"], allData["correlationFilter"], allData["metatreatmentsPassingCorrelationFilter"] = combineAndFilterCorrelations(config, allData["correlationCoefficients"], allData["correlationSigns"])

			allData["combinedCorrelationPValues"] = combineCorrelationPValues(config, allData["correlationPValues"].where(allData["metatreatmentsPassingCorrelationFilter"]))
			allData["correctedCorrelationPValues"] = correctCorrelationPValues(config, allData["combinedCorrelationPValues"].where(allData["correlationFilter"], other=numpy.NaN))

		def filterOnCorrelations(allData):
			nonlocal skip
			if skip:
				return

			allData["diagonalFilter"] = filterDiagonals(allData["correctedCorrelationPValues"])
			allData["individualCorrelationPValueFilter"] = filterOnIndividualCorrelationPValues(config, allData["correlationPValues"])
			allData["combinedCorrelationPValueFilter"] = filterOnCombinedCorrelationPValues(config, allData["combinedCorrelationPValues"])
			allData["correctedCorrelationPValueFilter"] = filterOnCorrectedCorrelationPValues(config, allData["correctedCorrelationPValues"])

			allData["edgeFilter"] = allData["diagonalFilter"] & allData["individualCorrelationPValueFilter"] & allData["combinedCorrelationPValueFilter"] & allData["correctedCorrelationPValueFilter"] & allData["correlationFilter"]

		def filterToExpectedEdges(allData):
			nonlocal skip
			if skip:
				return

			allData["expectedEdgeFilter"], allData["foldChangeSignProducts"] = filterOnExpectedEdges(config, allData["combinedFoldChangeSigns"], allData["combinedCorrelationSigns"])
			allData["edgeFilter"] &= allData["expectedEdgeFilter"]

		def createEdgeList(allData):
			nonlocal skip
			if skip:
				return

			allData["edges"] = allData["combinedCorrelationSigns"].where(allData["edgeFilter"], other=0)

		allData = {}
		allData["originalData"] = data
		stages = [computeDifferences, filterOnDifferences, computeCorrelations, filterOnCorrelations, filterToExpectedEdges, createEdgeList]
		return self.runPipeline(stages, allData, **kwargs)

