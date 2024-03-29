"""
Author: Nolan K Newman <newmanno@oregonstate.edu>
Last updated: 7/19/23

Written/tested in Python v3.8.10

Description:
Takes a network and calculates the number of shortest paths between two user-defined subnetworks. Output is a comma separated list of source 
and target nodes, along with their pathlength. For each pair of nodes, there can be multiple shortest paths. Each of these can be retained if 
specified by the user.

"""
import pickle 
import argparse
import networkx as nx
import numpy as np
import re
import csv
import os

if __name__ == '__main__':
    
    class dictionary(dict):
        def __init__(self):
            self = dict()

        def add(self, key, value):
            self[key] = value
            
    def connected_component_subgraphs(network):
    	return [network.subgraph(component) for component in nx.connected_components(network)]
    
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

    def assign_node_type(node_list_file, gc_nodes, type1, type2):
        '''
        Function that takes as input the node_type list from the user and creates a dictionary of the type for each node. Then, for the nodes
        that are in the netowrk input file, it assigns types to them based on the typing from the node_type file. It then outputs 2 dictionaries
        of nodes. 
        
         Arguments:
             - node_list_file: input file from user
             - gc_nodes: list of nodes in the giant component
             - type1: the type of nodes in group 1
             - type2: the type of nodes in group 2
        '''    
        
        subnet_dict = dictionary()
        
        # Empty lists to hold all the OTUs in the node_type input file from 
        # the command line
        type1_list = []
        type2_list = []
        
        node_type_dict = dictionary()
        
        # Add all node-type pairs from the input file into the node_type_dict
        with open(node_list_file) as node_file:
            node_file = csv.reader(node_file, delimiter = ',')
            
            for row in node_file:
                node_type_dict.add(row[0],row[1])
                
        # ensure the node groups supplied are both present in the mapping file
        for nodeGroup in args.node_groups:
            if nodeGroup not in list(node_type_dict.values()):
                raise Exception(f"Specified node group \"{nodeGroup}\" not in node map")
                
        # Generate a list of each specified node type in the giant component. Search the previously created dictionary and, for each 
        # value of type1 (i.e. 'otu') in the second column of the input file, assign the corresponding key to type1_list, then do the 
        # same thing for each type2 (i.e. 'pheno') value and its corresponding list
        for key,value in node_type_dict.items():
            try:
                if re.search(type1, value):
                    type1_list.append(key)
                elif re.match(type2, value):
                    type2_list.append(key)
            except:
                print("Unexpected value in the 'type' column of node_type input file.")

        # From the otu_list/pheno_list, only take the nodes that are present in the giant component. This is what the 
        # intersect1d function does. This prevents the need to generate  a new node type file for every network.                
        type1_for_dict = np.intersect1d(type1_list, gc_nodes)
        type2_for_dict = np.intersect1d(type2_list, gc_nodes)

        subnet_dict.add("Type1", type1_for_dict)
        subnet_dict.add("Type2", type2_for_dict)
        return(subnet_dict)

    def subnet_shortest_paths_length(G, spl_subnet1_node, spl_subnet2_node):
        sp_pair = spl_subnet1_node + "," + spl_subnet2_node
        sp_len = nx.shortest_path_length(G, spl_subnet1_node, spl_subnet2_node)
        sp_out = sp_pair + "," + str(sp_len)
        return sp_out
        
    def number_shortest_paths(G, nsp_subnet1_node, nsp_subnet2_node):
        all_sp = nx.all_shortest_paths(G, nsp_subnet1_node, nsp_subnet2_node)
        return all_sp
        
    # Main driver of the code
    parser = argparse.ArgumentParser(description="Example command: python ./analysis/find_all_shortest_paths_bw_subnets.py --network <file.pickle> --network-format <format> --map <map.csv> --node-groups <group1> <group2> --out-dir <directory>", add_help=False)

    requiredArgGroup = parser.add_argument_group('Required arguments')        
    requiredArgGroup.add_argument("--network", type=str, help="The path to the network file, either in .pickle or .csv format; see --network-format", required=True)
    requiredArgGroup.add_argument("--network-format", type=str, dest="network_format", choices=['pickle', 'csv'], help="Network file; Either use 'pickle' if you have already created the network.pickle file from calc_network_properties.py, or 'csv' if the network was reconstructed using an alternative pipeline (must be in .csv format and have 'partner1' and 'partner2' as the headers for the two node columns)", required=True)   
    requiredArgGroup.add_argument("--map", type=str, help="Mapping file (csv) of nodes and their subnetworks", required=True)
    requiredArgGroup.add_argument("--node-groups", nargs = 2, dest="node_groups", help = "The two groups in the mapping file you want to find the shortest paths between", required=True)
    requiredArgGroup.add_argument("--out-dir", type=str, dest = "outdir", help="Path to output directory", required=True)

    optionalArgGroup  = parser.add_argument_group('Optional arguments')        
    optionalArgGroup.add_argument("-h", "--help", action="help", help="Show this help message and exit")

    args = parser.parse_args()

    # set command arguments as variables
    map = args.map
    node_type1 = args.node_groups[0]
    node_type2 = args.node_groups[1]

    outdir = args.outdir
    if not outdir[-1] == "/":
        outdir = outdir + "/"
        
    if not os.path.exists(outdir):
        os.makedirs(outdir)  
        
    # Load in the network
    if args.network_format == 'pickle':
        # Unpack the pickle
        p = open(args.network, "rb")
        p = pickle.load(p)
        G = p
        
        network_name = args.network[:-7] # remove the pickle extension from file name    
        network_name = network_name.split('/')[-1]
        
    elif args.network_format == 'csv':
        G = import_outside_nw(args.network)
        
        network_name = args.network[:-4] # remove the csv extension from file name      
        network_name = network_name.split('/')[-1]

    filepath_nw_name = outdir + network_name
    print(filepath_nw_name)

    # Find only the giant component
    gc = max(connected_component_subgraphs(G), key=len)
    gc_nodes = gc.nodes()

    # find the nodes belonging to the two subnetworks
    two_subnets = assign_node_type(map, gc_nodes, node_type1, node_type2) 
    
    gc_pairs = len(two_subnets['Type1']) * len(two_subnets['Type2'])
    
    print("\nNumber of nodes in " + node_type1 + " group in the giant component: " + str(len(two_subnets['Type1'])))
    print("Number of nodes in " + node_type2 + " group in the giant component: " + str(len(two_subnets['Type2'])))
    print("Total number of pairs in the giant component: " + str(gc_pairs))

    with open(filepath_nw_name + "_shortest_path_bw_" + node_type1 + "_and_" + node_type2 + "_results.csv", "w") as out_file:
        
        # Find the shortest path length and number of shortest paths between each pair of nodes in the two subnetworks
        out_file.write(node_type1 + "_node," + node_type2 + "_node," + "Shortest_path_length,number_of_shortest_paths\n")
        for i in two_subnets['Type1']:
            for j in two_subnets['Type2']:
                spath = subnet_shortest_paths_length(G, i, j) # returns a string of pairs and their pathlength 
                sp_nodes = number_shortest_paths(G, i, j) # returns multiple generator objects of paths
                sp_len_and_num_sp = spath + "," + str(len(list(sp_nodes))) # Use len(list()) to find how many generators were created, then append
                out_file.write(sp_len_and_num_sp + "\n") # write to file

# =============================================================================
#   NOTE: This functionality should work, but will be implemented later.
#
#         # Do the same thing, but this time, find the nodes in the shortest path
#         out_file.write("\nNumber of shortest paths\n")
#         for x in two_subnets['Type1']:
#             for y in two_subnets['Type2']:
#                 sp_nodes = number_shortest_paths(G, x, y)
#                 num_sp = len(list(sp_nodes))
#                 out_file.write(x + "," + y + "," + str(num_sp) + "\n")
# =============================================================================
        
        # Overwrite any previous file with the same name instead of appending
        out_file.truncate()

print("\nFinished calculating all shortest paths. See results in " + filepath_nw_name + "_shortest_path_bw_" + node_type1 + "_and_" + node_type2 + "_results.csv\n")
