#!/usr/bin/env python3
import argparse
from pathlib import Path
from zipfile import ZipFile

def getArgs():
	parser = argparse.ArgumentParser(description="Analyze the results of network reconstruction for data subsamples.")
	parser.add_argument("dataFile", type=str, help="Input data ZIP containing subsample results")
	parser.add_argument("keepArrays", type=str, nargs="+", help="Array to keep from each subsample")
	parser.add_argument("filteredDataFile", type=str, help="Output data ZIP for the filtered subsample results")
	args = parser.parse_args()

	if args.dataFile == args.filteredDataFile:
		print("Error: Cannot read from and write to the same ZIP file. Exiting...")
		quit(1)

	args.dataFile = Path(args.dataFile)
	if not args.dataFile.exists():
		raise Exception("Specified data file does not exist")

	args.filteredDataFile = Path(args.filteredDataFile)

	return args

if __name__ == "__main__":
	args = getArgs()

	args.filteredDataFile.parent.mkdir(parents=True, exist_ok=True)
	filteredDataZip = ZipFile(args.filteredDataFile, mode="w")
	dataZip = ZipFile(args.dataFile)

	nameList = dataZip.namelist()
	subsampleIndices = {name.split("/")[0] for name in nameList}
	for index in subsampleIndices:
		arrayFound = False
		for keepArray in args.keepArrays:
			keepFileName = "{}/{}".format(index, keepArray)
			if keepFileName in nameList:
				arrayFound = True
				fileBuf = dataZip.read(keepFileName)
				filteredDataZip.writestr(keepFileName, fileBuf)
		if not arrayFound:
			# Write a blank file to keep a record of the subsample's existence
			filteredDataZip.writestr("{}/_".format(index), "_")
		
	dataZip.close()
	filteredDataZip.close()

