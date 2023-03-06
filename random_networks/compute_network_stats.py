#!/usr/bin/env python
import argparse
from collections import defaultdict
from community import community_louvain
import csv
from functools import partial
import json
from multiprocessing import Pool
import networkx as nx
import os
from pathlib import Path
import pickle
import shutil
import tempfile

# TODO: Logging

def getArgs():
	parser = argparse.ArgumentParser(description="Computes degree and BiBC for each node in the given networks", add_help=False)
	requiredArgGroup = parser.add_argument_group("required arguments")
	optionalArgGroup = parser.add_argument_group("optional arguments")
	requiredArgGroup.add_argument("--networks-file", type=str, dest="networksFile", required=True, help="Pickle file containing the networks to calculate stats for")
	requiredArgGroup.add_argument("--bibc-groups", choices=["node_types", "modularity"], dest="bibcGroups", required=True, help="Group nodes for BiBC based on type or modularity")
	requiredArgGroup.add_argument("--bibc-calc-type", choices=["rbc", "bibc"], dest="bibcCalcType", required=True, help="Compute raw BiBC or normalize (rbc)")
	requiredArgGroup.add_argument("--stats-file", type=str, dest="statsFile", required=True, help="Pickle file to output the network stats to")
	optionalArgGroup.add_argument("--node-map", "-m", dest="nodeMap", help="CSV file mapping nodes to their types. Required if node_types is specified for --bibc-groups.")
	optionalArgGroup.add_argument("--node-groups", "-g", metavar="GROUP", nargs=2, dest="nodeGroups", help="Two types of nodes to use for BiBC grouping. Required if node_types is specified for --bibc-groups.")
	optionalArgGroup.add_argument("-h", "--help", action="help", help="Show this help message and exit")
	args = parser.parse_args()

	args.networksFile = Path(args.networksFile)
	if not args.networksFile.exists():
		raise Exception("Specified network file does not exist")

	args.statsFile = Path(args.statsFile)

	if args.bibcGroups == "node_types":
		if args.nodeMap is None:
			raise Exception("If using node_types, must supply a node map with --node-map")
		if args.nodeGroups is None:
			raise Exception("If using node_types, must supply node groups with --node-groups")

	if args.nodeMap is not None:
		args.nodeMap = Path(args.nodeMap)
		if not args.nodeMap.exists():
			raise Exception("Specified node map does not exist")

	return args

# Partitions a network and returns the nodes in the two largest clusters.
def largestClusters(network):
	# Note: As per https://python-louvain.readthedocs.io/en/latest/api.html,
	# node evaluation order is random, so this is nondeterministic
	bestPartition = community_louvain.best_partition(network)

	clusters = defaultdict(list)
	for node, clusterId in bestPartition.items():
		clusters[clusterId].append(node)

	sortedClusters = sorted(clusters.values(), key=len, reverse=True)

	return sortedClusters[0], sortedClusters[1]

# Determines the nodes in the two specified types using the map file.
# Returns two lists of nodes, one for each type.
def nodesByType(args, network):
	nodeMap = defaultdict(list)
	allTypes = []
	with open(args.nodeMap) as nodeMapFile:
		nodeMapReader = csv.reader(nodeMapFile)
		for row in nodeMapReader:
			node, nodeType = row[0], row[1]
			allTypes.append(nodeType)
			if node in network.nodes:
				nodeMap[nodeType].append(node)

	for nodeGroup in args.nodeGroups:
		if nodeGroup not in allTypes:
			raise Exception(f"Specified node group \"{nodeGroup}\" not in node map")
		if nodeGroup not in nodeMap:
			raise Exception(f"Specified node group \"{nodeGroup}\" not present in the largest component of a network. This is likely due to high sparcity in the generated networks.")

	return nodeMap[args.nodeGroups[0]], nodeMap[args.nodeGroups[1]]

# Calculates restricted betweeness centrality for two groups of nodes.
# Returns a dictionary of node -> BiBC for the nodes in the specified groups.
def restrictedBetweennessCentrality(network, nodes1, nodes2, normalize):
	rbcs = defaultdict(int)
	for node1 in nodes1:
		for node2 in nodes2:
			# Betweenness centrality does not count the endpoints, so only use paths with intervening nodes
			paths = [path[1:-1] for path in list(nx.all_shortest_paths(network, node1, node2)) if len(path) > 2]
			for path in paths:
				for node in path:
					rbcs[node] += 1 / len(paths)

	if normalize:
		for node in rbcs:
			if node in nodes1:
				normalizationFactor = (len(nodes1) - 1) * len(nodes2)
			elif node in nodes2:
				normalizationFactor = len(nodes1) * (len(nodes2) - 1)
			else:
				normalizationFactor = len(nodes1) * len(nodes2)
			rbcs[node] /= normalizationFactor

	return rbcs

def calculateNetworkStats(args, networksTempDirPath, statsTempDirPath, networkName):
	with open(networksTempDirPath / networkName) as networkFile:
		networkData = json.load(networkFile)
	network = nx.readwrite.json_graph.adjacency_graph(networkData)

	sortedSubgraphs = sorted([network.subgraph(component) for component in nx.connected_components(network)], key=len, reverse=True)

	# If the two largest components are the same size, then don't calculate BiBC
	if len(sortedSubgraphs) > 1 and len(sortedSubgraphs[0]) == len(sortedSubgraphs[1]):
		bibcs = {node: None for node in network.nodes}
	# Otherwise, calculate BiBC
	else:
		giantComponent = sortedSubgraphs[0]
		normalize = args.bibcCalcType == "rbc"
		if args.bibcGroups == "node_types":
			group1Nodes, group2Nodes = nodesByType(args, giantComponent)
			bibcs = restrictedBetweennessCentrality(giantComponent, group1Nodes, group2Nodes, normalize)
		else: # "modularity"
			cluster1, cluster2 = largestClusters(giantComponent)
			bibcs = restrictedBetweennessCentrality(giantComponent, cluster1, cluster2, normalize)

	networkStats = defaultdict(dict)
	for node in network.nodes:
		networkStats[node]["degree"] = network.degree(node)
		networkStats[node]["bibc"] = bibcs[node] if node in bibcs else None

	with open(statsTempDirPath / networkName, "wb") as statsFile:
		pickle.dump(networkStats, statsFile)

if __name__ == "__main__":
	args = getArgs()

	networksTempDir = tempfile.TemporaryDirectory()
	networksTempDirPath = Path(networksTempDir.name)
	shutil.unpack_archive(args.networksFile, networksTempDirPath, "zip")

	statsTempDir = tempfile.TemporaryDirectory()
	statsTempDirPath = Path(statsTempDir.name)

	with Pool() as pool:
		pool.map(partial(calculateNetworkStats, args, networksTempDirPath, statsTempDirPath), os.listdir(networksTempDirPath))

	networksTempDir.cleanup()

	args.statsFile.parent.mkdir(parents=True, exist_ok=True)
	shutil.make_archive(args.statsFile.parent / args.statsFile.stem, "zip", statsTempDirPath)

	statsTempDir.cleanup()

