#!/usr/bin/env python3
import argparse
import itertools
import numpy
from pathlib import Path
from scipy import stats
import xarray
from zipfile import ZipFile

from util import CsvWriter

def getArgs():
	parser = argparse.ArgumentParser(description="Analyze the results of network reconstruction for data subsamples.", add_help=False)
	requiredArgGroup = parser.add_argument_group("required arguments")
	optionalArgGroup = parser.add_argument_group("optional arguments")
	requiredArgGroup.add_argument("--data-file", type=str, dest="dataFile", required=True, help="Input data ZIP containing subsample results")
	requiredArgGroup.add_argument("--out-file", type=str, dest="outFile", required=True, help="CSV file to write analysis results to")
	optionalArgGroup.add_argument("--singlecell", "-s", action="store_true", default=False, help="Work with single-cell rather than aggregate data")
	optionalArgGroup.add_argument("-h", "--help", action="help", help="Show this help message and exit")
	args = parser.parse_args()

	args.dataFile = Path(args.dataFile)
	if not args.dataFile.exists():
		raise Exception("Specified data file does not exist")

	args.outFile = Path(args.outFile)

	return args

if __name__ == "__main__":
	args = getArgs()

	dataZip = ZipFile(args.dataFile)
	nameList = dataZip.namelist()
	subsampleIndices = sorted({int(name.split("/")[0]) for name in nameList})
	neededArrays = {arrayName: [] for arrayName in ["correlationCoefficients", "combinedCorrelationPValues", "correctedCorrelationPValues", "edges"]}
	for index in subsampleIndices:
		for arrayName, arrayList in neededArrays.items():
			fileName = "{}/{}".format(index, arrayName)
			if fileName in nameList:
				arrayFile = dataZip.open(fileName)
				array = xarray.open_dataarray(arrayFile)
				arrayList.append(array)
				arrayFile.close()
			else:
				arrayList.append(None)
	dataZip.close()

	# This is a list of data array dicts. Each item corresponds to a subsample and is either None if one of the arrays in that subsample was None, or the full dict if they were all not None.
	# This is done so that a given subsample will only be included if all of the wanted data is present.
	filteredArrays = None
	for arrayName, arrayList in neededArrays.items():
		if filteredArrays is None:
			filteredArrays = [{arrayName: array} if array is not None else None for array in arrayList]
		else:
			for i in range(len(filteredArrays)):
				if filteredArrays[i] is not None:
					if arrayList[i] is None:
						filteredArrays[i] = None
					else:
						filteredArrays[i][arrayName] = arrayList[i]

	if sum(1 for arrayDict in filter(None.__ne__, filteredArrays)) == 0:
		print("WARNING: Missing arrays in all subsamples. No further processing performed.")
		quit()

	dataArrays = None
	for arrayDict in filteredArrays:
		if arrayDict is not None:
			if dataArrays is None:
				dataArrays = {arrayName: [array] for arrayName, array in arrayDict.items()}
			else:
				for arrayName, array in arrayDict.items():
					dataArrays[arrayName].append(array)

	rVals = xarray.concat(dataArrays["correlationCoefficients"], dim="subsample")
	combinedPVals = xarray.concat(dataArrays["combinedCorrelationPValues"], dim="subsample")
	correctedPVals = xarray.concat(dataArrays["correctedCorrelationPValues"], dim="subsample")
	edges = xarray.concat(dataArrays["edges"], dim="subsample")

	if args.singlecell:
		nodeDims = ["measurableAndCellType1", "measurableAndCellType2"]
	else:
		nodeDims = ["measurable1", "measurable2"]

	combinedPValMedians = combinedPVals.median(dim="subsample")
	correctedPValMedians = correctedPVals.median(dim="subsample")

	edgeIncluded = edges != 0
	subsampleAxis = edgeIncluded.get_axis_num("subsample")
	percentInclusions = xarray.apply_ufunc(numpy.count_nonzero, edgeIncluded, input_core_dims=[edgeIncluded.dims], output_core_dims=[nodeDims], kwargs=dict(axis=subsampleAxis)) / edgeIncluded.sizes["subsample"]

	combinedData = rVals.mean(dim="metatreatment")
	medians = combinedData.median(dim="subsample")
	maxs = combinedData.max(dim="subsample")
	mins = combinedData.min(dim="subsample")

	# See https://en.wikipedia.org/wiki/68%E2%80%9395%E2%80%9399.7_rule
	# 1, 2, 3, 4 standard deviations
	intervals = [0.6827, 0.9545, 0.9974, 0.9999]
	intervalStrings = [str(interval * 100) for interval in intervals]
	intervalLowerBoundData = [combinedData.quantile((1 - interval) / 2, dim="subsample") for interval in intervals]
	intervalUpperBoundData = [combinedData.quantile(1 - ((1 - interval) / 2), dim="subsample") for interval in intervals]
	intervalLowerBoundColumns = [CsvWriter.Column("{}% Confidence Interval Lower Bound".format(intervalString), intervalString + "_lower") for intervalString in reversed(intervalStrings)]
	intervalUpperBoundColumns = [CsvWriter.Column("{}% Confidence Interval Upper Bound".format(intervalString), intervalString + "_upper") for intervalString in intervalStrings]

	percentileOf0 = xarray.apply_ufunc(lambda values: stats.percentileofscore(values, 0) / 100, combinedData, input_core_dims=[["subsample"]], vectorize=True)

	data = {"percentInclusions": percentInclusions, "combinedPValMedians": combinedPValMedians, "correctedPValMedians": correctedPValMedians, "medians": medians, "maxs": maxs, "mins": mins, "percentileOf0": percentileOf0}
	for i in range(len(intervals)):
		data[intervalStrings[i] + "_lower"] = intervalLowerBoundData[i]
		data[intervalStrings[i] + "_upper"] = intervalUpperBoundData[i]

	dataColumns = [
		CsvWriter.Column("Percent Inclusion", "percentInclusions"),
		CsvWriter.Column("Median Combined p-Value", "combinedPValMedians"),
		CsvWriter.Column("Median Corrected p-Value", "correctedPValMedians"),
		CsvWriter.Column("Min", "mins"),
		*intervalLowerBoundColumns,
		CsvWriter.Column("Median", "medians"),
		*intervalUpperBoundColumns,
		CsvWriter.Column("Max", "maxs"),
		CsvWriter.Column("Percentile of 0", "percentileOf0")
	]

	coords = list(itertools.combinations(medians.coords[nodeDims[0]].data, 2))
	if args.singlecell:
		csvConfig = CsvWriter.Config(nodeDims,
			CsvWriter.PropertiesFormatted("Variable 1", "medians", "{}_{}", ["measurable1", "cellType1"]),
			CsvWriter.PropertiesFormatted("Variable 2", "medians", "{}_{}", ["measurable2", "cellType2"]),
			CsvWriter.PropertiesFormatted("Edge Name", "medians", "{}_{}<==>{}_{}", ["measurable1", "cellType1", "measurable2", "cellType2"]),
			CsvWriter.Property("Measurable 1", "medians", "measurable1"),
			CsvWriter.Property("Measurable 2", "medians", "measurable2"),
			CsvWriter.Property("Cell Type 1", "medians", "cellType1"),
			CsvWriter.Property("Cell Type 2", "medians", "cellType2"),
			CsvWriter.PropertiesFormatted("Edge Type", "medians", "{}<==>{}", ["cellType1", "cellType2"]),
			*dataColumns
		)
	else:
		csvConfig = CsvWriter.Config(nodeDims,
			CsvWriter.CoordinateFunction("Measurable 1", lambda m1, m2: m1),
			CsvWriter.CoordinateFunction("Measurable 2", lambda m1, m2: m2),
			CsvWriter.CoordinateFormatted("Edge Name", "{}<==>{}"),
			*dataColumns
		)
	CsvWriter.writeCsv(args.outFile, csvConfig, data, coords)

