#!/usr/bin/env python
import argparse
import json
import numpy
from pathlib import Path
import random
import xarray

def getArgs():
	parser = argparse.ArgumentParser(description="Create subsamples of a dataset.", add_help=False)
	requiredArgGroup = parser.add_argument_group("required arguments")
	optionalArgGroup = parser.add_argument_group("optional arguments")
	requiredArgGroup.add_argument("--dataFile", type=str, required=True, help="File in NetCDF format containing expression data")
	requiredArgGroup.add_argument("--subsampleFile", type=str, required=True, help="File to write the subsample descriptions to")
	requiredArgGroup.add_argument("-p", "--proportion", type=float, required=True, help="Proportion of samples to put in each subsample. Must satisfy 0 < p < 1.")
	requiredArgGroup.add_argument("-n", "--nsubsamples", type=int, required=True, help="Number of subsamples to generate.")
	optionalArgGroup.add_argument("-s", "--singlecell", action="store_true", default=False, help="Create subsamples for single-cell data instead of aggregate data")
	optionalArgGroup.add_argument("-o", "--organism", action="store_true", default=False, help="Subsample over organisms, even for single-cell data (this is the default behavior for aggregate data)")
	optionalArgGroup.add_argument("-h", "--help", action="help", help="Show this help message and exit")
	args = parser.parse_args()

	args.dataFile = Path(args.dataFile)
	if not args.dataFile.exists():
		raise Exception("Specified data file does not exist.")

	args.subsampleFile = Path(args.subsampleFile)

	if args.proportion <= 0 or args.proportion >= 1:
		raise Exception("Proportion of samples out of bounds. Must satisfy 0 < p < 1.")

	return args

if __name__ == "__main__":
	args = getArgs()

	subsamples = []
	if args.singlecell:
		data = xarray.open_dataset(args.dataFile)
		if args.organism:
			organisms = numpy.array(list(set(data["organism"].data)))
			nAllSamples = len(organisms)
			nSamples = round(args.proportion * nAllSamples)
			for i in range(args.nsubsamples):
				organismSubsample = random.sample(range(nAllSamples), nSamples)
				cellSubsample = []
				for organismIndex in organismSubsample:
					organism = organisms[organismIndex]
					organismCellIndices = numpy.argwhere((data["cellData"].organism == organism).data).flatten().tolist()
					cellSubsample.extend(organismCellIndices)
				cellSubsample.sort()
				subsamples.append(cellSubsample)
		else:
			for i in range(args.nsubsamples):
				coordSubsample = []
				def subsamplePopulation(data):
					nAllSamples = data.sizes["cell"]
					nSamples = round(args.proportion * nAllSamples)
					populationCoordSubsample = random.sample(list(data.coords["cell"].data), nSamples)
					coordSubsample.extend(populationCoordSubsample)
					return data
				data["cellData"].groupby("organism").map(lambda individualData: individualData.groupby("cellType").map(subsamplePopulation))
				indexSubsample = list(numpy.transpose(numpy.argwhere(numpy.isin(data.cell, coordSubsample))).astype(object)[0])
				subsamples.append(indexSubsample)
	else:
		data = xarray.open_dataarray(args.dataFile)
		nAllSamples = data.sizes["organism"]
		nSamples = round(args.proportion * nAllSamples)
		for i in range(args.nsubsamples):
			subsample = random.sample(range(nAllSamples), nSamples)
			subsample.sort()
			subsamples.append(subsample)

	subsampleFile = open(args.subsampleFile, "w")
	json.dump(subsamples, subsampleFile)
	subsampleFile.close()

