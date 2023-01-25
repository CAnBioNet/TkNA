#!/usr/bin/env python
import argparse
from pathlib import Path

from intake import intakeAggregateData, intakeSingleCellData

def getArgs():
	parser = argparse.ArgumentParser(description="Take in data in various formats and output it as a file in NetCDF.", add_help=False)
	requiredArgGroup = parser.add_argument_group("required arguments")
	optionalArgGroup = parser.add_argument_group("optional arguments")
	requiredArgGroup.add_argument("--data-dir", type=str, dest="dataDir", required=True, help="Directory to read data from")
	requiredArgGroup.add_argument("--out-file", type=str, dest="outFile", required=True, help="File to output data to in NetCDF")
	optionalArgGroup.add_argument("--singlecell", "-s", action="store_true", default=False, help="Intake single-cell data instead of aggregate data")
	optionalArgGroup.add_argument("-h", "--help", action="help", help="Show this help message and exit")
	args = parser.parse_args()

	args.dataDir = Path(args.dataDir)
	if not args.dataDir.exists():
		raise Exception("Specified data directory does not exist")
	elif not args.dataDir.is_dir():
		raise Exception("Specified data directory is not a directory")

	args.outFile = Path(args.outFile)

	return args

if __name__ == "__main__":
	args = getArgs()
	if args.singlecell:
		data = intakeSingleCellData(args.dataDir)
	else:
		data = intakeAggregateData(args.dataDir)

	args.outFile.parent.mkdir(parents=True, exist_ok=True)
	data.to_netcdf(args.outFile)

