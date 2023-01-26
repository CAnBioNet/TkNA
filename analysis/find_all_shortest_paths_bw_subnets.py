# This script will take as input two subnetworks, then calculate all shortest paths between 
# each pair of nodes in the two subnetworks. Output is a comma separated list of source 
# and target nodes, along with their pathlength. For each pair of nodes, there can be 
# multiple shortest paths. Each of these can be retained if specified by the user.
import pickle 
import argparse
import networkx as nx
import numpy as np
import re
import csv
# from networkx import all_shortest_paths

parser = argparse.ArgumentParser(description="Example command: python find_all_shortest_paths_bw_subnets.py <network.pickle> --node_map <map.csv> --node_groups gene pheno", add_help=False)

requiredArgGroup  = parser.add_argument_group('Required arguments')        
requiredArgGroup.add_argument("--network", type=str, help="inputs_for_downstream_plots.pickle file output by dot_plots.py", required=True)
requiredArgGroup.add_argument("--node-map", type=str, dest="node_map", help="Mapping file (csv) of nodes and their subnetworks", required=True)
requiredArgGroup.add_argument("--node-groups", nargs = 2, dest="node_groups", help = "The two groups in the mapping file you want to find the shortest paths between", required=True)

optionalArgGroup  = parser.add_argument_group('Optional arguments')        
optionalArgGroup.add_argument("-h", "--help", action="help", help="Show this help message and exit")

# parser = argparse.ArgumentParser(description='Example: python find_all_shortest_paths_bw_subnets.py <network.pickle> --node_map <map.csv> --node_groups gene pheno')
# parser.add_argument("network", help = 'pickled network file from assess_network.py')
# parser.add_argument("--node_map", help = 'Mapping file (csv) of nodes and their subnetworks')
# parser.add_argument("--node_groups", nargs = 2, help = 'The two groups in the mapping file you want to find the shortest paths between')

args = parser.parse_args()

# set command arguments as variables
net = args.network
map = args.node_map
node_type1 = args.node_groups[0]
node_type2 = args.node_groups[1]

if __name__ == '__main__':
    
    class dictionary(dict):
        def __init__(self):
            self = dict()

        def add(self, key, value):
            self[key] = value

    # Read in the pickle file, which should already have the network saved into it
    p = open(net, "rb")
    p = pickle.load(p)
    
    G = p

    # Function that takes as input the node_type list from the user and creates a 
    # dictionary of the type for each node. Then, for the nodes that are in the 
    # network input file, it assigns types to them based on the typing from the 
    # node_type file. It then outputs 2 dictionaries of nodes, one for each of the
    # user-specifed node types 
    def assign_node_type(node_list_file, gc_nodes, type1, type2):
        
        subnet_dict = dictionary()
        
        # Empty lists to hold all the OTUs in the node_type input file from 
        # the command line
        type1_list = []
        type2_list = []
        
        node_list = list(gc_nodes)
        node_type_dict = dictionary()
        
        # Add all node-type pairs from the input file into the node_type_dict
        with open(node_list_file) as node_file:
            node_file = csv.reader(node_file, delimiter = ',')
            
            for row in node_file:
                node_type_dict.add(row[0],row[1])
        
        # Generate a list of each specified node type in the giant component.
        # Search the previously created dictionary and, for each value of type1 (i.e. 'otu') 
        # in the second column of the input file, assign the corresponding 
        # key to type1_list, then do the same thing for each type2 (i.e. 'pheno') value and 
        # its corresponding list
        for key,value in node_type_dict.items():
            try:
                if re.search(type1, value):
                    type1_list.append(key)
                elif re.match(type2, value):
                    type2_list.append(key)
            except:
                print("Unexpected value in the 'type' column of node_type input file.")

        # From the otu_list/pheno_list, only take the nodes that are present 
        # in the giant component. This is what the intersect1d function does. 
        # This is because I don't want to generate a new node type file for 
        # every network and this way I can keep using the same one.
        # node_pairs_whole_nw = len(type1_list) * len(type2_list)       
                
        # print("\nNumber of nodes in " + node_type1 + " group: " + str(len(type1_list)))
        # print("Number of nodes in " + node_type2 + " group: " + str(len(type2_list)))
        # print("Total number of pairs between the two groups in the entire network: " + str(node_pairs_whole_nw))

        type1_for_dict = np.intersect1d(type1_list, gc_nodes)
        #print("\nCommon nodes between node_type1 input and giant component:")
        #print(type1_for_dict)

        type2_for_dict = np.intersect1d(type2_list, gc_nodes)
        #print("\nCommon nodes between node_type2 input and giant component:")
        #print(type2_for_dict)

        # Add the nodes from node_type1 and node_type2 that are exclusive to 
        # the network to their respective dictionaries, then return the 
        # dictionaries and use them to call the 
        # restricted_betweenness_centrality function
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
        

    # Find only the giant component
    gc = max(nx.connected_component_subgraphs(G), key=len)
    gc_nodes = gc.nodes()

    # Call assign_node_types to find the nodes belonging to the two subnetworks
    # takes mapping file, the nodes in the gc, and the type of each node 
    # (from the mapping file) as input, then assigns 
    two_subnets = assign_node_type(map, gc_nodes, node_type1, node_type2) 
    
    gc_pairs = len(two_subnets['Type1']) * len(two_subnets['Type2'])
    
    print("\nNumber of nodes in " + node_type1 + " group in the giant component: " + str(len(two_subnets['Type1'])))
    print("Number of nodes in " + node_type2 + " group in the giant component: " + str(len(two_subnets['Type2'])))
    print("Total number of pairs in the giant component: " + str(gc_pairs))

    with open("shortest_path_bw_" + node_type1 + "_and_" + node_type2 + "_results.txt", "w") as out_file:
        
        # Find the shortest path length and number of shortest paths between each pair of nodes in the two subnetworks
        out_file.write(node_type1 + "_node," + node_type2 + "_node," + "Shortest_path_length,number_of_shortest_paths\n")
        for i in two_subnets['Type1']:
            for j in two_subnets['Type2']:
                spath = subnet_shortest_paths_length(G, i, j) # returns a string of pairs and their pathlength 
                sp_nodes = number_shortest_paths(G, i, j) # returns multiple generator objects of paths
                sp_len_and_num_sp = spath + "," + str(len(list(sp_nodes))) # Use len(list()) to find how many generators were created, then append
                out_file.write(sp_len_and_num_sp + "\n") # write to file


        # Do the same thing, but this time, find the nodes in the shortest path
        #out_file.write("\nNumber of shortest paths\n")
        #for x in two_subnets['Type1']:
        #    for y in two_subnets['Type2']:
        #        sp_nodes = number_shortest_paths(G, x, y)
        #        num_sp = len(list(sp_nodes))
        #        out_file.write(x + "," + y + "," + str(num_sp) + "\n")

        
        # Overwrite any previous file with the same name instead of appending
        out_file.truncate()

print("\nFinished calculating all shortest paths. See shortest_path_bw_" + node_type1 + "_and_" + node_type2 + "_results.txt\n")
