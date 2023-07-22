"""
Author: Nolan K Newman <newmanno@oregonstate.edu>
Last updated: 7/19/23

Written/tested in Python v3.8.10

Description:
Takes a network and uses Infomap to assign cluster numbers. Infomap clustering is performed for each data type in the network.

"""

import pickle
import argparse
import networkx as nx
from infomap import Infomap
import csv
import collections
import os
    
if __name__ == '__main__':

    def import_outside_nw(fname):
        '''
        Function to import the network into script if network was reconstructed via alternative methods other than TkNA
       
        Arguments:
            - fname: the name of the file
        '''
        row_count = 0 
        G = nx.Graph()
        
        with open(fname) as csvfile:     
            file = csv.reader(csvfile, delimiter = ',')
                    
            for row in file:
            
                #print(row)
                
                # Take the index of the source and target node in the header of the file
                if row_count == 0: 
                    p1 = int(row.index("partner1"))
                    
                    p2 = int(row.index("partner2"))
    
                parameter1 = row[p1]
                parameter2 = row[p2]  
                
                # Find each node that made it into the final network and the direction of the edge   
                if row_count != 0:
                    G.add_edge(parameter1, parameter2)
        
                row_count += 1
    
        csvfile.close()                
        return(G)
    
        # Function that takes as input the mapping file from the user and creates a dictionary of the type for each node. Then, for the nodes that are in the network input file, it assigns types to them based on the typing from the mapping file. 
    def assign_node_type(node_list_file):         
        '''
        Function that takes as input the mapping file from the user and creates a dictionary of the type for each node. Then, for the nodes that are in the network input file, it assigns types to them based on the typing from the mapping file. 
        
        Arguments:
            - node_list_file: the mapping file supplied by the user
        '''
        
        node_type_dict = collections.defaultdict(list)

        # Add all node-type pairs from the input file into the node_type_dict
        with open(node_list_file) as node_file:
            node_file = csv.reader(node_file, delimiter = ',')
            
            for row in node_file:
                node_type_dict[row[1]].append(row[0])
                
        return(node_type_dict)
 

    def infomap_partition(G, type_dict, subnw):
        '''
        Assigns nodes to subnetworks within user-defined subnetworks 
        
        Arguments:
            - G: the entire network
            - type_dict: dictionary that includes the node names as keys and the user-defined subnetwork as values
            - subnw: the user-defined subnetwork that is currently being analyzed
        '''
        im = Infomap()
        # make node-to-int and int-to-node dictionaries
        
        # Extract just the nodes from the subnetwork being analyzed
        subnw_nodes = type_dict[subnw]        
        subnw_graph = nx.subgraph(G, subnw_nodes)
        
        j = 0
        node_to_int = {}
        int_to_node = {}
        for n in subnw_graph.nodes():
            node_to_int[n] = j
            int_to_node[j] = n
            j += 1
            
        # copy the edges into InfoMap
        for e in subnw_graph.edges():
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
    
    
    parser = argparse.ArgumentParser(description="Example command: python ./analysis/infomap_assignment.py --network <file> --network-format <format> --map <file.csv> --out-dir <directory>", add_help=False)

    requiredArgGroup = parser.add_argument_group('Required arguments')  
    requiredArgGroup.add_argument("--network", type=str, help="The path to the network file, either in .pickle or .csv format; see --network-format", required=True)
    requiredArgGroup.add_argument("--network-format", type=str, dest="network_format", choices=['pickle', 'csv'], help="Format of the network file; Either use 'pickle' with the network.pickle file output made by assess_network.py (if network was reconstructed using the TkNA pipeline) or 'csv' if the network was reconstructed using an alternative pipeline (must be in .csv format and have 'partner1' and 'partner2' as the headers for the two node columns)", required=True)
    requiredArgGroup.add_argument("--map", help = 'CSV file with the name of the node in the first column and its data type in the second column (i.e. ENSMUSG00000030708, gene).', required=True)
    requiredArgGroup.add_argument("--out-dir", type=str, dest = "outdir", help="Path to output directory", required=True)

    optionalArgGroup  = parser.add_argument_group('Optional arguments')        
    optionalArgGroup.add_argument("-h", "--help", action="help", help="Show this help message and exit")
        
    args = parser.parse_args()
    outdir = args.outdir
    
    if not outdir[-1] == "/":
        outdir = outdir + "/"
        
    if not os.path.exists(outdir):
        os.makedirs(outdir)  
        
    #Main driver of the code
    
    # Load in the network
    if args.network_format == 'pickle':
        # Unpack the pickle
        p = open(args.network, "rb")
        p = pickle.load(p)
        G = p
        
        network_name = args.network.split("/")[-1]
        network_name = network_name[:-7] # remove the pickle extension from file name 
        
    elif args.network_format == 'csv':
        G = import_outside_nw(args.network)
       
        network_name = args.network.split("/")[-1]   
        network_name = network_name[:-4] # remove the csv extension from file name     

    # Assign each node a node type so infomap is run separately on each type
    node_type_dict = assign_node_type(args.map)

    # dict of infomap outputs, keyed on subnetwork name and values are list of infomap_partition results
    im_outputs = {}

    with open(outdir + network_name + "_infomap_partition.csv", "w") as file:
        file.write("Node,Subnetwork_partition\n")
        
        for subnet in node_type_dict.keys():
            print("Performing infomap in the " + subnet + " subnetwork...")
            im_outputs[subnet] = infomap_partition(G, node_type_dict, subnet)
            [file.write(str(k) + "," + subnet + "_" + str(v) + "\n") for k,v in im_outputs[subnet][1].items()]
    file.close()            
        
    print("\nFile saved: " + outdir + network_name + "_infomap_partition.csv\n")    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
