#!/usr/bin/env python
import argparse
import csv
from functools import total_ordering
from pathlib import Path
import pickle

def getArgs():
	parser = argparse.ArgumentParser(description="Synthesize the stats for each network by determining the node with max BiBC and degree. By default, sorts first by BiBC and then by degree.")
	parser.add_argument("networkStatsFile", type=str, help="Pickle file containing the raw network stats")
	parser.add_argument("synthesizedStatsFile", type=str, help="CSV file to write the synthesized network stats to")
	parser.add_argument("--flip-priority", "-f", dest="flipPriority", action="store_true", help="Changes sorting method to use degree first and BiBC second")
	args = parser.parse_args()

	args.networkStatsFile = Path(args.networkStatsFile)
	if not args.networkStatsFile.exists():
		raise Exception("Specified network stats file does not exist")

	args.synthesizedStatsFile = Path(args.synthesizedStatsFile)

	return args

@total_ordering
class NoBibc:
	def __le__(self, other):
		return True

	def __eq__(self, other):
		return self is other

if __name__ == "__main__":
	args = getArgs()

	with open(args.networkStatsFile, "rb") as networkStatsFile:
		allStats = pickle.load(networkStatsFile)

	synthesizedStats = []
	for networkStats in allStats:
		# None can't be compared with ints, so replace it with instances of a special class that will always compare as lesser
		for node in networkStats:
			if networkStats[node]["bibc"] is None:
				networkStats[node]["bibc"] = NoBibc()

		if args.flipPriority:
			stat1, stat2 = "degree", "bibc"
		else:
			stat1, stat2 = "bibc", "degree"
		sortedStats = sorted(sorted(networkStats.items(), key=lambda nodeStats: nodeStats[1][stat1], reverse=True), key=lambda nodeStats: nodeStats[1][stat2], reverse=True)

		topNode = sortedStats[0]
		synthesizedStats.append({"node": topNode[0], "degree": topNode[1]["degree"], "bibc": topNode[1]["bibc"]})

	args.synthesizedStatsFile.parent.mkdir(parents=True, exist_ok=True)
	with open(args.synthesizedStatsFile, "w") as synthesizedStatsFile:
		synthesizedStatsWriter = csv.writer(synthesizedStatsFile)
		synthesizedStatsWriter.writerow(["Network Index", "Node", "Degree", "BiBC"])
		for index, networkStats in enumerate(synthesizedStats):
			synthesizedStatsWriter.writerow([index, networkStats["node"], networkStats["degree"], networkStats["bibc"]])

