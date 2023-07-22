#!/usr/bin/env python
import argparse
from functools import partial
import json
from math import factorial
from multiprocessing import Pool
import networkx as nx
import numpy as np
from pathlib import Path
import pickle
import shutil
import tempfile

def maxEdges(numNodes):
	return int(factorial(numNodes) / (2 * factorial(numNodes - 2)))

def generateNetwork(nodes, numEdges, tempDirPath, index, seed):
	network = nx.Graph()
	network.add_nodes_from(nodes)

	np.random.seed(seed)

	for _ in range(numEdges):
		while True:
			node1, node2 = np.random.choice(network.nodes, 2, replace=False)
			if not network.has_edge(node1, node2):
				network.add_edge(node1, node2)
				break

	networkData = nx.readwrite.json_graph.adjacency_data(network)
	with open(tempDirPath / str(index), "w") as networkFile:
		json.dump(networkData, networkFile)

def getArgs():
	parser = argparse.ArgumentParser(description="Generate random networks", add_help=False)
	inputArgGroup = parser.add_argument_group("input arguments")
	requiredArgGroup = parser.add_argument_group("required arguments")
	optionalArgGroup = parser.add_argument_group("optional arguments")

	inputArgGroup.add_argument("--template-network", type=str, dest="templateNetwork", help="Path to a file containing a pickled networkx network. Alternatively, use --node-list-file and --num-edges to specify the desired network properties.")
	inputArgGroup.add_argument("--node-list-file", type=str, dest="nodeListFile", help="Path to a file containing the names of the nodes in the generated networks. Must be one node per line. If used, --num-edges must also be provided. Alternatively, use --template-network to base the random networks on an existing one.")
	inputArgGroup.add_argument("--num-edges", type=int, dest="numEdges", help="Number of edges in each randomly-generated network. If used, --node-list-file must also be provided. Alternatively, use --template-network to base the random networks on an existing one.")

	requiredArgGroup.add_argument("--networks-file", type=str, dest="networksFile", required=True, help="ZIP file to output networks to")

	optionalArgGroup.add_argument("--num-networks", "-n", type=int, dest="numNetworks", default=10000, help="Number of networks to generate")
	optionalArgGroup.add_argument("--cores", type=int, help="Number of cores to use for computation. If not provided, all available cores will be used.")
	optionalArgGroup.add_argument("-h", "--help", action="help", help="Show this help message and exit")

	args = parser.parse_args()

	if args.templateNetwork is None and args.nodeListFile is None and args.numEdges is None:
		raise Exception("Must specify either --template-network or both --node-list-file and --num-edges")
	if args.templateNetwork is not None and (args.nodeListFile is not None or args.numEdges is not None):
		raise Exception("Cannot use both --template-network and --node-list-file/--num-edges")
	if args.nodeListFile is not None and args.numEdges is None:
		raise Exception("If using --node-list-file, must also specify --num-edges")
	if args.numEdges is not None and args.nodeListFile is None:
		raise Exception("If using --num-edges, must also specify --node-list-file")

	if args.templateNetwork is not None:
		args.templateNetwork = Path(args.templateNetwork)
		if not args.templateNetwork.exists():
			raise Exception("Specified template network does not exist")
	else:
		args.nodeListFile = Path(args.nodeListFile)
		if not args.nodeListFile.exists():
			raise Exception("Specified node list does not exist")

	args.networksFile = Path(args.networksFile)

	return args

if __name__ == "__main__":
	args = getArgs()

	if args.templateNetwork is not None:
		with open(args.templateNetwork, "rb") as templateNetworkFile:
			templateNetwork = pickle.load(templateNetworkFile)
		nodes = templateNetwork.nodes
		numEdges = len(templateNetwork.edges)
	else:
		with open(args.nodeListFile) as nodeListFile:
			nodes = nodeListFile.readlines()
		nodes = [node.strip() for node in nodes]
		nodes = [node for node in nodes if len(node) > 0]

		numMaxEdges = maxEdges(len(nodes))
		if args.numEdges > numMaxEdges:
			raise Exception(f"Too many edges desired (must be {numMaxEdges} or fewer for this node list)")
		numEdges = args.numEdges

	# Generate seeds ahead of time so that processes don't end up sharing seeds (and hence networks)
	seeds = np.random.randint(np.iinfo(np.int32).max, dtype=np.int32, size=args.numNetworks)

	tempDir = tempfile.TemporaryDirectory()
	tempDirPath = Path(tempDir.name)
	with Pool(args.cores) as pool:
		pool.starmap(partial(generateNetwork, nodes, numEdges, tempDirPath), enumerate(seeds))

	args.networksFile.parent.mkdir(parents=True, exist_ok=True)
	shutil.make_archive(args.networksFile.parent / args.networksFile.stem, "zip", tempDirPath)

	tempDir.cleanup()

