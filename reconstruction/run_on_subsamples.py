#!/usr/bin/env python3
import argparse
import json
from matplotlib import pyplot
import numpy
import pandas
from pathlib import Path
import tempfile
import xarray
from zipfile import ZipFile

from reconstruction import NetworkReconstructorAggregate, NetworkReconstructorSingleCell
from util import parseConfigFile, Dataset
from util.configs import aggregateConfigSpec, singleCellConfigSpec

def getArgs():
	parser = argparse.ArgumentParser(description="Generate causal networks for the specified subsamples of the data.", add_help=False)
	requiredArgGroup = parser.add_argument_group("required arguments")
	optionalArgGroup = parser.add_argument_group("optional arguments")
	requiredArgGroup.add_argument("--data-file", type=str, dest="dataFile", required=True, help="Input data file")
	requiredArgGroup.add_argument("--config-file", type=str, dest="configFile", required=True, help="JSON configuration file")
	requiredArgGroup.add_argument("--subsample-file", type=str, dest="subsampleFile", required=True, help="File specifying the subsamples to run on")
	requiredArgGroup.add_argument("--output-file", type=str, dest="outputFile", required=True, help="ZIP file to write out the output data to")
	optionalArgGroup.add_argument("--singlecell", "-s", action="store_true", default=False, help="Work with single-cell rather than aggregate data")
	optionalArgGroup.add_argument("-h", "--help", action="help", help="Show this help message and exit")
	args = parser.parse_args()

	args.dataFile = Path(args.dataFile)
	if not args.dataFile.exists():
		raise Exception("Specified data file does not exist")

	args.configFile = Path(args.configFile)
	if not args.configFile.exists():
		raise Exception("Specified config file does not exist")

	args.subsampleFile = Path(args.subsampleFile)
	if not args.subsampleFile.exists():
		raise Exception("Specified subsample file does not exist")

	args.outputFile = Path(args.outputFile)

	return args

def readJson(jsonFilePath):
	with open(jsonFilePath) as jsonFile:
		jsonContents = json.load(jsonFile)
	return jsonContents

if __name__ == "__main__":
	args = getArgs()

	configSpec = singleCellConfigSpec if args.singlecell else aggregateConfigSpec
	config = parseConfigFile(configSpec, args.configFile)

	subsamples = readJson(args.subsampleFile)

	dataset = Dataset()
	dataset.load_from_file(args.dataFile)

	if args.singlecell:
		network_reconstructor = NetworkReconstructorSingleCell()
	else:
		network_reconstructor = NetworkReconstructorAggregate()

	if args.singlecell:
		coordToSubsample = "cell"
	else:
		coordToSubsample = "organism"
	tablesToSubsampleNames = []
	otherTablesNames = []
	tableNames = dataset.get_table_names()
	for tableName in tableNames:
		if coordToSubsample in dataset.get_table(tableName).coords:
			tablesToSubsampleNames.append(tableName)
		else:
			otherTablesNames.append(tableName)

	allResults = []
	for subsample in subsamples:
		subsampledDataset = Dataset()
		for tableName in otherTablesNames:
			subsampledDataset.add_table(tableName, dataset.get_table(tableName))
		for tableName in tablesToSubsampleNames:
			tableToSubsample = dataset.get_table(tableName)
			subsampledTable = tableToSubsample.isel({coordToSubsample: subsample})
			subsampledDataset.add_table(tableName, subsampledTable)
		results = network_reconstructor.reconstructNetwork(config, subsampledDataset)
		allResults.append(results)

	args.outputFile.parent.mkdir(parents=True, exist_ok=True)
	dataOutZip = ZipFile(args.outputFile, mode="w")
	for i, results in enumerate(allResults):
		for dataName, dataArray in results.items():
			tempFile = tempfile.NamedTemporaryFile()
			dims = list(dataArray.dims)
			for dim in dims:
				if isinstance(dataArray.get_index(dim), pandas.MultiIndex):
					dataArray = dataArray.reset_index(dim)
			dataArray.to_netcdf(tempFile.name)
			dataOutZip.write(tempFile.name, arcname="{}/{}".format(i, dataName))
			tempFile.close()
	dataOutZip.close()

