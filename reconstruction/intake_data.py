#!/usr/bin/env python
import argparse
from pathlib import Path

from intake import intakeAggregateData, intakeSingleCellData

def getArgs():
	parser = argparse.ArgumentParser(description="Take in data in various formats and output it as a file in NetCDF.")
	parser.add_argument("dataDir", type=str, help="Directory to read data from")
	parser.add_argument("outFile", type=str, help="File to output data to in NetCDF")
	parser.add_argument("--singlecell", "-s", action="store_true", default=False, help="Intake single-cell data instead of aggregate data")
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

