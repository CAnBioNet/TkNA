#!/usr/bin/env python3
import argparse
import json
from pathlib import Path
import xarray

from intake import intakeAggregateData, intakeSingleCellData
from reconstruction import NetworkReconstructorAggregate, NetworkReconstructorSingleCell
from util import parseConfigFile
from util.configs import aggregateConfigSpec, singleCellConfigSpec

def getArgs():
	parser = argparse.ArgumentParser(description="Generate causal network.", add_help=False)
	requiredArgGroup = parser.add_argument_group("required arguments")
	optionalArgGroup = parser.add_argument_group("optional arguments")
	requiredArgGroup.add_argument("--data-source", type=str, dest="dataSource", required=True, help="Input data file/directory")
	requiredArgGroup.add_argument("--config-file", type=str, dest="configFile", required=True, help="JSON configuration file")
	requiredArgGroup.add_argument("--out-file", type=str, dest="outFile", required=True, help="Writes all generated data to a ZIP file with the specified path")
	optionalArgGroup.add_argument("--start", nargs=2, type=str, metavar=("startStage", "startData"), help="Start network generation from a particular stage")
	optionalArgGroup.add_argument("--stop", nargs=1, type=str, metavar="stopStage", help="Stop network generation at a particular stage")
	optionalArgGroup.add_argument("--cores", type=int, help="Number of cores to use for computation. If not provided, all available cores will be used.")
	optionalArgGroup.add_argument("--singlecell", "-s", action="store_true", default=False, help="Work with single-cell rather than aggregate data")
	optionalArgGroup.add_argument("-h", "--help", action="help", help="Show this help message and exit")
	args = parser.parse_args()

	args.dataSource = Path(args.dataSource)
	if not args.dataSource.exists():
		raise Exception("Specified data source does not exist")

	args.configFile = Path(args.configFile)
	if not args.configFile.exists():
		raise Exception("Specified config file does not exist")

	args.outFile = Path(args.outFile)

	if args.start:
		startFile = Path(args.start[1])
		if not startFile.exists():
			raise Exception("Specified start file does not exist")
		args.start = {"startStage": args.start[0], "startData": args.start[1]}

	if args.stop:
		args.stop = args.stop[0]

	return args

if __name__ == "__main__":
	args = getArgs()

	configSpec = singleCellConfigSpec if args.singlecell else aggregateConfigSpec
	config = parseConfigFile(configSpec, args.configFile)

	if args.dataSource.is_dir():
		if args.singlecell:
			data = intakeSingleCellData(args.dataSource)
		else:
			data = intakeAggregateData(args.dataSource)
	else:
		if args.singlecell:
			data = xarray.open_dataset(args.dataSource)
		else:
			data = xarray.open_dataarray(args.dataSource)

	if args.singlecell:
		network_reconstructor = NetworkReconstructorSingleCell()
	else:
		network_reconstructor = NetworkReconstructorAggregate()
	network_reconstructor.reconstructNetwork(config, data, start=args.start, stopStage=args.stop, dataOutFilePath=args.outFile, cores=args.cores)

