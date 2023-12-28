#!/usr/bin/env python
import argparse
import itertools
import json
import numpy
from pathlib import Path
import random

from util import Dataset

def getArgs():
	parser = argparse.ArgumentParser(description="Create subsamples of a dataset.", add_help=False)
	requiredArgGroup = parser.add_argument_group("required arguments")
	optionalArgGroup = parser.add_argument_group("optional arguments")
	requiredArgGroup.add_argument("--data-file", type=str, dest="dataFile", required=True, help="File in NetCDF format containing expression data")
	requiredArgGroup.add_argument("--subsample-file", type=str, dest="subsampleFile", required=True, help="File to write the subsample descriptions to")
	requiredArgGroup.add_argument("-p", "--proportion", type=float, required=True, help="Proportion of samples to put in each subsample. Must satisfy 0 < p < 1.")
	requiredArgGroup.add_argument("-n", "--nsubsamples", type=int, required=True, help="Number of subsamples to generate.")
	optionalArgGroup.add_argument("-s", "--singlecell", action="store_true", default=False, help="Create subsamples for single-cell data instead of aggregate data. By default, this will subsample cells rather than organisms. To override this behavior, use the -o/--organism option to subsample organisms or the -b/--both-organisms-and-cells option to subsample both.")
	optionalArgGroup.add_argument("-o", "--organism", action="store_true", default=False, help="Subsample over organisms, even for single-cell data (this is the default behavior for aggregate data)")
	optionalArgGroup.add_argument("-b", "--both-organisms-and-cells", dest="bothOrganismsAndCells", nargs=1, metavar="ORGANISM_PROPORTION", help="Subsample over organisms and cells (only applicable to single-cell, -s). The additional argument is the proportion of organisms in each treatment to use for each subsample and must satisfy 0 < p < 1.")
	optionalArgGroup.add_argument("-h", "--help", action="help", help="Show this help message and exit")
	args = parser.parse_args()

	args.dataFile = Path(args.dataFile)
	if not args.dataFile.exists():
		raise Exception("Specified data file does not exist.")

	args.subsampleFile = Path(args.subsampleFile)

	if args.proportion <= 0 or args.proportion >= 1:
		raise Exception("Proportion of samples out of bounds. Must satisfy 0 < p < 1.")

	bothOrganismsAndCells = args.bothOrganismsAndCells is not None

	if args.singlecell:
		if args.organism and bothOrganismsAndCells:
			raise Exception("Cannot use the -o/--organism and -b/--both-organisms-and-cells options simultaneously.")
	else:
		if args.organism:
			raise Exception("Can only use -o/--organism with single-cell data (-s/--singlecell)")
		if bothOrganismsAndCells:
			raise Exception("Can only use -b/--both-organisms-and-single-cells with single-cell data (-s/--singlecell).")

	if bothOrganismsAndCells:
		args.organismProportion = float(args.bothOrganismsAndCells[0])
		if args.organismProportion <= 0 or args.organismProportion >= 1:
			raise Exception("Proportion of organisms out of bounds. Must satisfy 0 < p < 1.")

	args.bothOrganismsAndCells = bothOrganismsAndCells

	return args

if __name__ == "__main__":
	args = getArgs()

	dataset = Dataset()
	dataset.load_from_file(args.dataFile)

	subsamples = []
	if args.singlecell:
		data = dataset.get_table("cellData")

		# Construct maps from dataset coordinates
		organismCellMap = {} # organism -> list of indices of cells belonging to the organism
		experimentMap = {} # experiment -> list of organisms in experiment
		treatmentMap = {} # treatment -> list of organisms in treatment group
		cellTypeMap = {} # cell type -> list of indices of cells of that type
		for index, cellData in enumerate(data.cell):
			organism = cellData.organism.item()
			if organism not in organismCellMap:
				organismCellMap[organism] = []

				experiment = cellData.experiment.item()
				if experiment not in experimentMap:
					experimentMap[experiment] = []
				experimentMap[experiment].append(organism)

				treatment = cellData.treatment.item()
				if treatment not in treatmentMap:
					treatmentMap[treatment] = []
				treatmentMap[treatment].append(organism)

			cellType = cellData.cellType.item()
			if cellType not in cellTypeMap:
				cellTypeMap[cellType] = []
			cellTypeMap[cellType].append(index)
			organismCellMap[organism].append(index)

		def subsampleCells(organisms):
			subsample = []
			for organism in organisms:
				organismCellIndices = set(organismCellMap[organism])
				for cellType in cellTypeMap:
					cellTypeCellIndices = set(cellTypeMap[cellType])
					matchingCellIndices = organismCellIndices.intersection(cellTypeCellIndices)
					subsetSubsampleSize = round(args.proportion * len(matchingCellIndices))
					subsetSubsample = random.sample(matchingCellIndices, subsetSubsampleSize)
					subsample.extend(subsetSubsample)
			subsample.sort()
			return subsample

		if args.organism or args.bothOrganismsAndCells:
			# Take a sample of the organisms from each treatment
			organismProportion = args.proportion if args.organism else args.organismProportion
			for i in range(args.nsubsamples):
				organismSubsample = []
				for experimentOrganisms in experimentMap.values():
					for treatmentOrganisms in treatmentMap.values():
						matchingOrganisms = set(experimentOrganisms).intersection(set(treatmentOrganisms))
						subsetSubsampleSize = round(organismProportion * len(matchingOrganisms))
						subsetSubsample = random.sample(matchingOrganisms, subsetSubsampleSize)
						organismSubsample.extend(subsetSubsample)
				if args.organism:
					subsample = list(itertools.chain(*[organismCellMap[organism] for organism in organismSubsample]))
					subsample.sort()
					subsamples.append(subsample)
				elif args.bothOrganismsAndCells:
					# Now that we have a sample of the organisms, take a sample of each selected organism's cells within each cell type
					subsamples.append(subsampleCells(organismSubsample))
				else:
					raise Exception("Argument state invalid.")
		else:
			for i in range(args.nsubsamples):
				subsamples.append(subsampleCells(organismCellMap.keys()))
	else:
		data = dataset.get_table("originalData")

		# Construct maps from dataset coordinates
		experimentMap = {} # experiment -> list of indices of organisms in experiment
		treatmentMap = {} # treatment -> list of indices of organisms in treatment group
		for index, organism in enumerate(data.organism):
			experiment = organism.experiment.item()
			if experiment not in experimentMap:
				experimentMap[experiment] = []
			experimentMap[experiment].append(organism)

			treatment = organism.treatment.item()
			if treatment not in treatmentMap:
				treatmentMap[treatment] = []
			treatmentMap[treatment].append(index)

		for i in range(args.nsubsamples):
			subsample = []
			for experimentOrganismIndices in experimentMap.values():
				for treatmentOrganismIndices in treatmentMap.values():
					matchingOrganismIndices = set(experimentOrganismIndices).intersection(set(treatmentOrganismIndices))
					subsetSubsampleSize = round(args.proportion * len(subsetOrganismIndices))
					subsetSubsample = random.sample(matchingOrganismIndices, subsetSubsampleSize)
					subsample.extend(subsetSubsample)
			subsample.sort()
			subsamples.append(subsample)

	subsampleFile = open(args.subsampleFile, "w")
	json.dump(subsamples, subsampleFile)
	subsampleFile.close()

