# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 23:13:06 2020

@author: Nolan Newman <newmanno@oregonstate.edu>

Input: pickled network file created by assess_network.py

Output: CSV file containing the name of each node and which subnetwork (dictated by a number) the node was assigned via louvain

"""

from community import community_louvain
import pickle
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Example command: python louvain_partition.py <file.pickle>", add_help=False)

    requiredArgGroup = parser.add_argument_group('Required arguments')   
    requiredArgGroup.add_argument("--pickle", type=str, help="network.pickle file output by assess_network.py", required=True)
    
    
    # parser = argparse.ArgumentParser(description='Example: python louvain_partition.py <pickled network file>')
    # parser.add_argument("pickle", help = 'pickled network file from import_network_data.py')
    
    # parser.add_argument("--number", help = 'only a one-use instance where this is the number of the file being created in a loop')
    args = parser.parse_args()
    
    def louvain_assignment(G):
        p = community_louvain.best_partition(G, random_state=1, randomize = False)
        #[print(k,v) for k,v in p.items()]
        Q = community_louvain.modularity(p,G)
        pq = [p,Q]
        return(pq)        
        
    #Main driver of the code
    # Unpack the pickle
    p = open(args.pickle, "rb")
    p = pickle.load(p)
        
    G = p
    
    network_name = args.pickle[:-7] # remove the pickle extension from file name        

    lou = louvain_assignment(G)

    with open(network_name + "_louvain_partition.csv", "w") as file:
        file.write("Node,partition\n")
        [file.write(str(k) + "," + str(v) + "\n") for k,v in lou[0].items()]
    file.close()
    
    # with open(network_name + "_louvain_Q_" + str(args.number) + ".txt", "a") as qfile:
    #     qfile.write(str(lou[1]) + "\n")
    # qfile.close()
    
    print("\nModularity (Q) for this network:")
    print(str(lou[1]) + "\n")

    print("\nFile saved: " + network_name + "_louvain_partition.csv\n")