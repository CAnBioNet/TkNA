#!/usr/bin/env python
import argparse
from functools import partial
from math import factorial
from multiprocessing import Pool
import networkx as nx
import numpy as np
from pathlib import Path
import pickle

def maxEdges(numNodes):
	return int(factorial(numNodes) / (2 * factorial(numNodes - 2)))

def generateNetwork(nodes, numEdges, seed):
	network = nx.Graph()
	network.add_nodes_from(nodes)

	np.random.seed(seed)

	for _ in range(numEdges):
		while True:
			node1, node2 = np.random.choice(network.nodes, 2, replace=False)
			if not network.has_edge(node1, node2):
				network.add_edge(node1, node2)
				break

	return network

def getArgs():
	parser = argparse.ArgumentParser(description="Generate random networks")
	parser.add_argument("nodeListFile", type=str, help="Path to a file containing the names of the nodes in the generated networks. Must be one node per line.")
	parser.add_argument("numEdges", type=int, help="Number of edges in each randomly-generated network")
	parser.add_argument("networksFile", type=str, help="File to output pickled network list to")
	parser.add_argument("--numNetworks", "-n", type=int, default=10000, help="Number of networks to generate")
	args = parser.parse_args()

	args.nodeListFile = Path(args.nodeListFile)
	if not args.nodeListFile.exists():
		raise Exception("Specified node list does not exist")

	args.networksFile = Path(args.networksFile)

	return args

if __name__ == "__main__":
	args = getArgs()

	with open(args.nodeListFile) as nodeListFile:
		nodes = nodeListFile.readlines()
	nodes = [node.strip() for node in nodes]
	nodes = [node for node in nodes if len(node) > 0]

	numMaxEdges = maxEdges(len(nodes))
	if args.numEdges > numMaxEdges:
		raise Exception(f"Too many edges desired (must be {numMaxEdges} or fewer for this node list)")

	# Generate seeds ahead of time so that processes don't end up sharing seeds (and hence networks)
	seeds = np.random.randint(np.iinfo(np.int32).max, dtype=np.int32, size=args.numNetworks)

	with Pool() as pool:
		networks = pool.map(partial(generateNetwork, nodes, args.numEdges), seeds)

	args.networksFile.parent.mkdir(parents=True, exist_ok=True)
	with open(args.networksFile, "wb") as networksFile:
		pickle.dump(networks, networksFile)
