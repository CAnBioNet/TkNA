# -*- coding: utf-8 -*-
"""
Created on Thu Jan 5 23:13:06 2023

@author: Nolan Newman (infomap function by Kevin Brown)

Input: pickled network file created by assess_network.py

Output: CSV file containing the name of each node and which subnetwork (dictated by a number) the node was assigned via infomap
	
"""

import pickle
import argparse
import networkx as nx
from infomap import Infomap
    
if __name__ == '__main__':

    #def infomap_partition(G,n_mod=0):
    def infomap_partition(G):
        '''
        Wrapper to make networkx graph input compatible with the infomap
        package, which just calls wrapped C code and has an annoying API
        (i.e. does not talk to networkx directly, does not allow arbitrary
        hashable types as nodes, etc.)
        '''
        im = Infomap()
        # make node-to-int and int-to-node dictionaries
        j = 0
        node_to_int = {}
        int_to_node = {}
        for n in G.nodes():
            node_to_int[n] = j
            int_to_node[j] = n
            j += 1
        # copy the edges into InfoMap
        for e in G.edges():
            im.add_link(node_to_int[e[0]],node_to_int[e[1]])
        # now run in silent mode
        #options_string = '--silent --preferred-number-of-modules '+str(n_mod)
        options_string = '--silent'
        im.run(options_string)
        # set up the node->community id dictionary
        partition = {}
        for node in im.tree:
            if node.is_leaf:
                partition[int_to_node[node.node_id]] = node.module_id - 1
        return im.codelength,partition
    
    parser = argparse.ArgumentParser(description='Example: python infomap_assignment.py <pickled network file>')
    parser.add_argument("pickle", help = 'pickled network file from import_network_data.py')
    args = parser.parse_args()
    
    #Main driver of the code
    # Unpack the pickle
    p = open(args.pickle, "rb")
    p = pickle.load(p)
       
    G = p
    
    network_name = args.pickle[:-7] # remove the pickle extension from file name        

    im_output = infomap_partition(G)

    print(im_output)

    with open(network_name + "_infomap_partition.txt", "w") as file:
        file.write("Node,partition\n")
        [file.write(str(k) + "," + str(v) + "\n") for k,v in im_output[1].items()]
    file.close()
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
