"""
Author: Nolan K Newman <newmanno@oregonstate.edu>
Last updated: 1/5/23

Written in Python v3.5.3

Description:
Takes as input the pickled network file created from import_network_data.py and calculates various node and network properties for it

Developer note: As of 9/29/19, the code no longer relies strictly on positional arguments and instead uses an argparser
Developer note: As of 2/19/20, the script now has flags to specify whether PUC and BiBC should be calculated. Also added
                eigenvalue centrality as another node property. On top of that, this version now incorporates positive    
                and negative ratios of edges in the giant component now.
Developer note: As of 6/29/20, I removed the dictionary class since I don't think it served much of a purpose. Now I just
                add to and modify dictionaries like normal. Additionally, I added an argument where the user can specify whether 
                they would like to calculate node fragmentation. This was needed due to how time intensive that property is to 
                calculate.
Developer note: As of 7/17/20, infomap was added as an option to maximize modularity
Developer note: As of 7/28/20, the code takes a pickled input file of both a networkx graph object and a dictionary of deviation from expected values, not calculated
                with import_network_data.py
"""

import pickle
import networkx as nx
import numpy as np 
import re
import argparse
#import community
from community import community_louvain
from statistics import mean, median
from sortedcontainers import SortedDict
import csv
from collections import defaultdict
from networkx import all_shortest_paths
from collections import OrderedDict
import copy
import os
import matplotlib.pyplot as plt

#from infomap import Infomap

####### Get user input ########

parser = argparse.ArgumentParser(description="Example command: python calc_network_properties.py <pickled network file> --bibc --bibc-groups node_types --bibc-calc-type bibc --node-map <node map csv> --node_groups micro pheno", add_help=False)

requiredArgGroup = parser.add_argument_group('Required arguments')        
requiredArgGroup.add_argument("--pickle", type=str, help="The pickle file created with import_network_data.py", required=True)

optionalArgGroup = parser.add_argument_group('Optional arguments') 
optionalArgGroup.add_argument("-h", "--help", action="help", help="Show this help message and exit")
optionalArgGroup.add_argument("--frag", help = 'Flag; Do you want to compute node fragmentation centrality? (Significantly increases run-time)', action = 'store_true')
optionalArgGroup.add_argument("--bibc", help= 'Flag; Do you want to compute BiBC? (Significantly increases run-time)', action = 'store_true')
optionalArgGroup.add_argument("--bibc-groups", dest = "bibc_groups", choices = ['node_types', 'modularity'], help= "What to compute BiBC on, either distinct groups (node_types) or on the two most modular regions (modularity) of the network (found using the Louvain method). Required if --bibc is set")
optionalArgGroup.add_argument("--bibc-calc-type", dest="bibc_calc_type", choices = ['rbc', 'bibc'], help= "Would you like to normalize the bibc value based on amount of nodes in each group (rbc) or not (bibc)? Required if --bibc is set.")
optionalArgGroup.add_argument("--node-map", type=str, dest="node_map", help= "Required if node_types is specified for --bibc-groups. CSV of nodes and their types (i.e. otu, pheno, gene, etc.)")
optionalArgGroup.add_argument("--node-groups", nargs = 2, type=str, dest="node_groups", help= "2 args; Required if node_types is specified for --bibc-groups. The two groups of nodes to calculate BiBC/RBC on")



# parser = argparse.ArgumentParser(description='Example: python calc_network_properties.py <pickled network file> --bibc --bibc_groups node_types --bibc_calc_type bibc --node_map <node map csv> --node_groups micro pheno')

# # Required args
# # pickle output from import_network_data2.py
# parser.add_argument('pickle', help = 'The pickle file created with import_network_data.py')

# # Flags and optional arguments
# parser.add_argument("--frag", help = 'Flag; Do you want to compute node fragmentation centrality? (Significantly increases run-time)', action = 'store_true')
# parser.add_argument("--bibc", help = 'Flag; Do you want to compute BiBC? (Significantly increases run-time)', action = 'store_true')
# parser.add_argument("--bibc_groups", choices = ['node_types', 'modularity'], help = 'What to compute BiBC on, either distinct groups or on the two most modular regions of the network (found using the Louvain method)')
# parser.add_argument("--bibc_calc_type", choices = ['rbc', 'bibc'], help = 'Would you like to normalize based on amount of nodes in each group (rbc) or not (bibc)?')
# parser.add_argument("--node_map", help = 'Required if node_types is specified for --bibc_groups. CSV of nodes and their types (i.e. otu, pheno, gene, etc.)')
# parser.add_argument("--node_groups", nargs = 2, help = '2 args; Required if node_types is specified for --bibc_groups. Its the two groups of nodes to calculate BiBC/RBC on')
# parser.add_argument("--subnw_mean_deg", help = 'Required if you wish to calculate subnetwork mean degree. Also requires a mapping file to be supplied with --node_map')

args = parser.parse_args()

if args.bibc:
    if args.bibc_groups == "node_types":
        bibc_choice = "node_types"
        node_input_file = args.node_map
        node_type1 = args.node_groups[0]
        node_type2 = args.node_groups[1]
    elif args.bibc_groups == "modularity":
        bibc_choice = "modularity"
        
    bibc_calc_type = args.bibc_calc_type

def connected_component_subgraphs(network):
	return [network.subgraph(component) for component in nx.connected_components(network)]

if __name__ == '__main__':

    # Unpack the pickle
    p = open(args.pickle, "rb")
    p = pickle.load(p)
    G = p
    
    class dictionary(dict):
        def __init__(self):
            self = dict()
        def add(self, key, value):
            self[key] = value 
    
    # Function that calculates the number of second neighbors for each node
    def second_neighbors(G,n):
        nn = []
        for x in list([y for y in G.neighbors(n)]):
            nn.append([y for y in G.neighbors(x)])
        flat_nn = [item for sublist in nn for item in sublist]
        return [i for i in flat_nn if i not in G.neighbors(n) and i !=  n]

    # Function that takes as input the node_type list from the user and creates a dictionary of the type for each node. Then, for the nodes
    # that are in the netowrk input file, it assigns types to them based on the typing from the node_type file. It then outputs 2 dictionaries
    # of nodes, one that includes phenotypes and one that contains OTUs. These then get passed to the restricted_betweenness_centrality function,
    # where rbc is calculated for each node.
    def assign_node_type(node_list_file, gc_nodes, type1, type2):

        otu_and_pheno_dict = {}

        # Empty lists to hold all the OTUs in the node_type input file from the command line
        type1_list = []
        type2_list = []            

        node_type_dict = {}

        # Add all node-type pairs from the input file into the node_type_dict
        with open(node_list_file) as node_file:
            node_file = csv.reader(node_file, delimiter = ',')
            
            for row in node_file:
                node_type_dict[row[0]] = row[1]

        # Search the previously created dictionary and, for each 'otu' value in the second column of the input file, assign 
        # the corresponding key to otu_list, then do the same thing for each 'pheno' value and its corresponding list 
        for key,value in node_type_dict.items():
            try:
                if re.search(type1, value):
                    type1_list.append(key)
                elif re.match(type2, value):
                    type2_list.append(key)
            except:
                print("Unexpected value in the 'type' column of node_type input file.")

        # From the otu_list/pheno_list, only take the nodes that are present in the giant component. This is what the intersect1d function does. 
        # This is because I don't want to generate a new node type file for every network and this way I can keep using the same one. 

        type1_for_dict = np.intersect1d(type1_list, gc_nodes)
        print("\nCommon nodes between node_type1 input and giant component:")
        print(type1_for_dict)        

        type2_for_dict = np.intersect1d(type2_list, gc_nodes)
        print("\nCommon nodes between node_type2 input and giant component:")
        print(type2_for_dict)        

        # Add the nodes from node_type1 and node_type2 that are exclusive to the network to their respective dictionaries, 
        # then return the dictionaries and use them to call the restricted_betweenness_centrality function
        otu_and_pheno_dict['Type1'] = type1_for_dict
        otu_and_pheno_dict['Type2'] = type2_for_dict

        return(otu_and_pheno_dict)
    
    # Take only the nodes in the giant component and assign each to a subnetwork based on the node mapping file supplied by
    # the user. Function returns a dictionary of nodes keyed by their subnetwork name
    def create_node_dict(node_list_file, gc_nodes):
            
            subnet_dict = dictionary()
            node_list = list(gc_nodes)        
            node_type_dict = dictionary()
            uniq_subnw = []
            
            # Add all node-type pairs from the input file into the subnet_dict
            with open(node_list_file) as node_file:
                node_file = csv.reader(node_file, delimiter = ',')
    
                for row in node_file:
                    parameter = row[0]
                    sub = row[1]
                    
                    node_type_dict.add(parameter,sub)
    
                    # Generate a unique list of user-specified subnetworks
                    if sub not in uniq_subnw: 
                        uniq_subnw.append(sub)
            
                # Iterate through all the subnetworks and generate a dictionary keyed on subnw name
                for subnet in uniq_subnw:
                    templist = [] # Temporary list to hold unique nodes for each subnw 
                    
                    for k,v in node_type_dict.items():
                       
                        # If the node is part of the GC...
                        if k in node_list:
                        
                            # and if the current subnw in the list is the same as the the current subnetwork in
                            # the dictionary then add the list key (the node) to the list of nodes in that subnw
                            if subnet == v: 
                                #print(subnet + " and " + v + " are the same.")
                                templist.append(k)
                                #print(templist)
                                subnet_dict.add(subnet, templist) # subnet is current subnw in list, k is node name
                
            return(subnet_dict)


    # Function that finds which nodes belong to the two most modular portions of the giant component, then returns
    # those nodes as a dictionary. These then get passed to the restricted_betweenness_centrality function below.
    def bibc_mod(nodes_from_gc):
        # Split the gc into 'clusters', clustering by the modularity
        part = community_louvain.best_partition(nodes_from_gc)
        
        # The previous method returns a dictionary, where keys are nodes and
        # values are the cluster they belong in. We want to find the two largest 
        # clusters, so we first find which unique clusters there are, then 
        # parse through those and assign nodes
        unique_clusters = set(part.values())

        # Add the clusters to a dictionary, where keys are clusters and values are nodes in that cluster
        mod_dict = {}

        for val in unique_clusters:
            nodes_in_cluster = []
            for key,value in part.items():
                if value == val:
                    nodes_in_cluster.append(key)
    
            mod_dict[val] = nodes_in_cluster

        # Sort the mod_dict by the length (size) of the clusters
        sorted_clust = sorted(mod_dict, key = lambda k: len(mod_dict[k]), reverse = True)
    
        # make lists of nodes in the largest clusters. sorted_clust[0] is the key of the largest cluster in mod_dict 
        large_mod = mod_dict[sorted_clust[0]]
        second_large_mod = mod_dict[sorted_clust[1]]

        # Add the lists to a dictionary, which then get returned to be passed to restricted_betweenness_centrality()
        return_mod = {}
        return_mod['mod1'] = large_mod
        return_mod['mod2'] = second_large_mod
        
        print("\n Network partitions")
        print("Largest GC partition (", str(len(return_mod['mod1'])), " nodes):")
        print(return_mod['mod1'])
        print("\n")
        print("Second largest GC partition (", str(len(return_mod['mod2'])), " nodes):")    
        print(return_mod['mod2'])
        print("\n")

        return(return_mod)        
    
    # Calculates BiBC (more correctly called restricted betweenness centrality) for each node
    def restricted_betweenness_centrality(G,nodes_0,nodes_1,type):
        '''
        Restricted betweenness centrality that only computes centrality using paths with sources
        in nodes_0 and targets in nodes_1 (or vice versa, which double counts).

        Returns three dictionaries of centralities, one for nodes_0, one for nodes_1,
        and one for the other nodes in G not in either set.

        If one thinks carefully about normalization:
        -centrality values for nodes in group 0 should be divided by (N0-1)N1
        -group 1: N0(N1-1)
        -others: N0*N1

        (Should just be able to normalize by len(nodes_0)*len(nodes_1))
        '''
        
        flatten = lambda l: [item for sublist in l for item in sublist]
        rbc = defaultdict(int)
        for s in nodes_0:
            for t in nodes_1:
                # might want to avoid putting the whole thing in memory
                # betweenness centrality does not count the endpoints (v not in s,t)
                paths_st = [x for x in list(all_shortest_paths(G,s,t)) if len(x) > 2]
                n_paths = len(paths_st)
                nodes_to_update = flatten([p[1:-1] for p in paths_st])
                for n in nodes_to_update:
                    rbc[n] += 1/n_paths
        
        # split the dictionary in three
        rbc_n0 = {}.fromkeys(nodes_0)
        rbc_n1 = {}.fromkeys(nodes_1)
        rbc_other = {}
        
        for n in G.nodes():
            if n in nodes_0:
                rbc_n0[n] = rbc[n]
            elif n in nodes_1:
                rbc_n1[n] = rbc[n]
            else:
                rbc_other[n] = rbc[n]
    
        # If the user specifies rbc as the bibc calculation type, then normalize each node
        if type == "rbc":
            # Normalize each node - code suggested by Kevin and added by Nolan
            for i in rbc_n0:
                rbc_n0[i] = rbc_n0[i]/((len(nodes_0)-1) * len(nodes_1))
    
            for i in rbc_n1:
                rbc_n1[i] = rbc_n1[i]/((len(nodes_1)-1) * len(nodes_0))
    
            for i in rbc_other:
                rbc_other[i] = rbc_other[i]/(len(nodes_0) * len(nodes_1))

        elif type == "bibc":
            print("Normalization was not conducted on the set of nodes in the rbc function. BiBC was calculated instead.")

        # Return the list of BiBC of each node from each of the two groups
        return rbc_n0,rbc_n1,rbc_other

    def distance_fragmentation(G):
        '''
        Another fragmentation measure of Borgatti's, this one based on distance:
    
                dF = 1 - sum_ij(1/d_ij)/N(N-1)
        '''
        N = len(G.nodes())
        sum_inv_dij = 0.0
        for n in nx.all_pairs_shortest_path_length(G):
            for k in n[1]:
                if k != n[0]:
                    sum_inv_dij += 1.0/n[1][k]
        if N != 1:
            dF = 1 - sum_inv_dij/(N*(N-1))
        else:
            dF = 1
        return dF
    
    def plot_deg_dist(G):
        '''
        Creates a png of the degree distribution of the nodes in the reconstructed network
    
        This code was obtained from yatu's post here: https://stackoverflow.com/questions/53958700/plotting-the-degree-distribution-of-a-graph-using-nx-degree-histogram'
    
        '''        
        
        degree_freq = nx.degree_histogram(G)
        degrees = range(len(degree_freq))
        plt.figure(figsize=(12, 8)) 
        plt.loglog(degrees, degree_freq,'go-', linestyle='None')  
        plt.xlabel('Degree')
        plt.ylabel('Frequency')
     
    # Deprecated, now plotted in dot_plots.py
    # def plot_deg_bibc(G):
    #     '''
    #     Creates a png of the degree-BiBC distribution of the nodes in the reconstructed network
        
    #     '''        
        
    #     degree_freq = nx.degree_histogram(G)
    #     degrees = range(len(degree_freq))
    #     plt.figure(figsize=(12, 8)) 
    #     plt.loglog(degrees, degree_freq,'go-', linestyle='None')  
    #     plt.xlabel('Degree')
    #     plt.ylabel('Frequency')

    # Deprecated, now calculated in infomap_assignment.py
    # def infomap_partition(G,n_mod=0):
    #     '''
    #     Wrapper to make networkx graph input compatible with the infomap
    #     package, which just calls wrapped C code and has an annoying API
    #     (i.e. does not talk to networkx directly, does not allow arbitrary
    #     hashable types as nodes, etc.)
    #     '''
    #     im = Infomap()
    #     # make node-to-int and int-to-node dictionaries
    #     j = 0
    #     node_to_int = {}
    #     int_to_node = {}
    #     for n in G.nodes():
    #         node_to_int[n] = j
    #         int_to_node[j] = n
    #         j += 1
    #     # copy the edges into InfoMap
    #     for e in G.edges():
    #         im.add_link(node_to_int[e[0]],node_to_int[e[1]])
    #     # now run in silent mode
    #     options_string = '--silent --preferred-number-of-modules '+str(n_mod)
    #     im.run(options_string)
    #     # set up the node->community id dictionary
    #     partition = {}
    #     for node in im.tree:
    #         if node.is_leaf:
    #             partition[int_to_node[node.node_id]] = node.module_id - 1
    #     return im.codelength,partition

    ################################################################################
    ######################## Calculate network properties ##########################
    ################################################################################

    network_name = args.pickle[:-7]

    filedir = os.path.dirname(os.path.abspath(args.pickle))

    with open(filedir + "/network_properties.txt", "w") as file:

        #------------------------------------------------#
        ###             Network properties             ###
        #------------------------------------------------#
        file.write("### Network properties ###\n")
    
        ### Number of nodes ###
        print("Finding number of nodes...")
        nnodes = nx.number_of_nodes(G)
        nnodes_str = str(nx.number_of_nodes(G))
        file.write("Number_of_nodes\t" + nnodes_str + "\n")
    
        ### Number of edges ###
        print("Finding number of edges...")
        nedges = nx.number_of_edges(G)
        nedges_str = str(nx.number_of_edges(G))
        file.write("Number_of_edges\t" + nedges_str + "\n")

        ### Mean degree ###
        print("Finding mean degree...")
        mean_degree = str((2*nedges)/nnodes)
        file.write("Mean_degree\t" + mean_degree + "\n")

        ### Average clustering coefficient ###
        # The following code works exactly like the nx method, but the nx version was not outputting any values at first
        print("Finding average clustering coefficient...")
        clust_coef = nx.clustering(G)
        cc_holder = []
        for key,value in clust_coef.items():
            cc_holder.append(value)
        
        file.write("Average_clustering_coefficient\t" + str(round(mean(cc_holder), 5)) + "\n")
        #file.write("Average_clustering_coefficient\t" + str(nx.average_clustering(G)) + "\n")
        

        ### Mean geodesic path length ###
        # Finds ASPL for each node in the giant component, then find the average of those nodes
        #print("Finding average shortest path...")
        gc = max(connected_component_subgraphs(G), key=len)
        #ASPL = nx.average_shortest_path_length(gc)
        #ASPL = ASPL + str(a) + "\t"  
        #file.write("Mean_geodesic_path_length\t" + str(round(ASPL, 5)) + "\n")    

        ### Giant component ###
        print("Finding size of giant component...")
        giant = len(max(connected_component_subgraphs(G), key=len))/len(G)
        file.write("Giant_component\t" + str(round(giant, 5)) + "\n")

        ### Number of components ###
        print("Finding number of components...")
        ncc = str(nx.number_connected_components(G))
        file.write("Number_of_connected_components\t" + ncc + "\n")
    
        ### Freeman centrality ###
        print("Finding freeman centrality...")
        max_degree = 0
        for i in G:
            if G.degree[i] > max_degree:
                max_degree = G.degree[i]
        
        FC_list = []
        if nnodes > 2:
            for i in G:
                FC_node = (max_degree - G.degree[i])/((nnodes - 1) * (nnodes - 2))
                FC_list.append(FC_node)
            freeman_cent = sum(FC_list)
            file.write("Freeman_centralization\t" + str(round(freeman_cent, 5)) + "\n")    
        
        else:
            file.write("Freeman_centralization\tnan \n")    
        
        ### Mean closeness centrality ###
        print("Finding mean closeness centrality...")
        closeness_cent = nx.closeness_centrality(G)
        num_holder = []
        for key,value in closeness_cent.items():
            num_holder.append(value)
        mean_closeness_cent = mean(num_holder)
        file.write("Mean_closeness_centrality\t" + str(round(mean_closeness_cent, 5)) + "\n")
        
        ### Modularity ###
        print("Finding modularity...")
        p = community_louvain.best_partition(G, random_state=1, randomize = False)
        Q = community_louvain.modularity(p,G)
        print("Modularity of best partition of graph: ", str(Q))
        file.write("Modularity\t" + str(round(Q, 5)) + "\n")

        ### Median comp size over total number of nodes ###
        print("Finding median component size over number of nodes...")
        med_list = []
        for g in connected_component_subgraphs(G):
            med_list.append(nx.number_of_nodes(g))
        med_over_nnodes = median(med_list)/nnodes
        file.write("Median_comp_size_over_#_nodes\t" + str(round(med_over_nnodes, 5)) + "\n")

        ### Degree assortativity ###
        print("Finding degree assortativity...")
        degree_assortativity = nx.degree_assortativity_coefficient(G)
        file.write("Degree_assortativity\t" + str(round(degree_assortativity, 5)) + "\n") 

        ### Network fragmentation ###
        print("Finding network fragmentation...")
        '''
            Obtained from Borgatti's key player problem paper, equation 3
            This equation specifically takes into account the size of the components
            rather than just the number of components. But it does not take into 
            account distances, which is what dF (sdistance based fragmentation) does.
            We calculate dF later on as a node property
        '''
        sum_numerator = 0
        for comp in nx.connected_components(G):
            comp_size = len(comp)
            sum_numerator += comp_size*(comp_size-1)
        F = 1 - (sum_numerator/(nnodes * (nnodes - 1)))         
        file.write("Network_fragmentation\t " + str(F))
        
        # Overwrite any previous file with same name instead of appending    
        file.truncate()
                
    file.close()
    

        ################################################################################
        ########################## Calculate subnetwork properties #####################
        ################################################################################       
    with open(filedir + "/subnetwork_properties.txt", "w") as file:

        # mean_deg_subnw = mean_degree_each_subnet_move_argparse_1_2_2023.subnw_mean_degree(G, node_input_file)
        
        # print("\n\n\nMEAN DEGREE\n")
        # print(mean_deg_subnw)
        # print("\n\n\n")
       
        file.write("\n### Subnetwork properties ###\n")
        print("Finding subnetwork mean degree...")
        
        # Find only the giant component
        gc = max(connected_component_subgraphs(G), key=len)
        gc_nodes = gc.nodes()
       
        # Assign nodes to subnetworks from the mapping file
        subnets = create_node_dict(args.node_map, gc_nodes)
        
        meandeg_dict = {}
        for k,v in subnets.items():
            total_degree = 0
            sg = G.subgraph(v)
            for i in sg.nodes():
                total_degree += sg.degree[i]
            mean_degree = total_degree / nx.number_of_nodes(sg)
            meandeg_dict[k] = mean_degree # Make a new key with the name of k and the value of mean_degree
            file.write(str(k) + "_mean_degree\t" + str(mean_degree) + "\n")
    
        # Overwrite any previous file with same name instead of appending    
        file.truncate()
                
    file.close()
        
        ################################################################################
        ########################## Calculate node properties ###########################
        ################################################################################
    with open(filedir + "/node_properties.txt", "w") as file:
    
        ###    Empty cell   ###
        file.write("name\t")    

        ###    Node list    ###
        node_names = list(G.nodes)
        node_names_sort = sorted(node_names)    
        [file.write(i + "\t") for i in node_names_sort]
        file.write("\n")    

        ###   Node degree   ###
        print("Finding each node's degree...")
        node_degrees = ""
        for i in node_names_sort:
            node_degrees = node_degrees + str(G.degree(i)) + "\t"  
        file.write("Node_degrees\t" + node_degrees + "\n")    

        ### Node strength ###
        print("Finding each node's strength...")
        strength = ""
        for i in node_names_sort:
            strength = strength + str(round(G.degree(i,weight='weight'), 5)) + "\t"
        file.write("Node_strength\t" + strength + "\n")
    
        ###   Node clustering   ###
        # checked using https://www.e-education.psu.edu/geog597i_02/node/832
        print("Finding each node's clustering coefficient...")
        node_clust = ""
        for i in node_names_sort:
            node_clust = node_clust + str(round(nx.clustering(G,i), 5)) + "\t"  
        file.write("Node_clustering\t" + node_clust + "\n")    

        ###   Node closeness   ###
        print("Finding each node's closeness...")
        node_close = ""
        for i in node_names_sort:
            node_close = node_close + str(round(nx.closeness_centrality(G,i), 5)) + "\t"
        file.write("Node_closeness\t" + node_close + "\n")    
        
        ###   Eigenvector centrality   ###
        # Note: un-commented this on 2/19/2020. Had issues with it before. Trying again.
        print("Finding each node's eigenvector centrality...")
        ecen_dict = nx.eigenvector_centrality(G)
        ecen_dict_sorted = SortedDict(ecen_dict)
        ecen = ""
        for key,value in ecen_dict_sorted.items():
            ecen = ecen + str(round(value, 5)) + "\t"    
        file.write("Eigenvalue_centrality\t" + ecen + "\n")

        ###   Betweenness centrality   ###
        '''
            The numeber of times a single node appears on the shortest path between all other pairs of 
            nodes in a network
        '''
        print("Finding each node's betweeness centrality...")
        bc_dict = nx.betweenness_centrality(G)
        bc_dict_sorted = SortedDict(bc_dict)
        bc = ""
        for key,value in bc_dict_sorted.items():
            bc = bc + str(round(value, 5)) + "\t"
        file.write("Betweenness_centrality\t" + bc + "\n")


        ### Fragmentation centrality ###
        '''
            This is similar to the network fragmentation calculated previously, but in this case, 
            this property takes into account the distances between the nodes and what happens after
            each node of the network is removed. It is technically a network property that gets 
            calculated for each iteration of removing nodes and then that value gets assigned to the
            node that was removed from the network that iteration
        '''

        if args.frag:
            print("\nCalculating node fragmentation centrality. This might be a while...")
            node_frag_dict = {}
    
            for i in node_names_sort: # For each node...
                cumul_div = 0
                H = copy.deepcopy(G)
                
                #sum_term = 0.0
                #frag_denominator = nnodes * (nnodes - 1)
                
                #print("All nodes in new graph H after removing node", i)
                H.remove_node(i) # Remove current node
    
                node_frag = distance_fragmentation(H)         
                node_frag_dict[i] = node_frag
    
                #print("The node fragmentation value for", i, "is", node_frag)
    
            frag = ""
            for key,value in SortedDict(node_frag_dict).items():
                frag = frag + str(value) + "\t"
                
            file.write("Node_fragmentation\t" + frag + "\n")
    
            # For testing purposes
            dist_frag = distance_fragmentation(G)
            print("Distance fragmentation for the entire network: " + str(dist_frag))


        ### Number of second neighbors ###
        print("Finding each node's number of second neighbors...")
        nsn_list = ""
        for i in node_names_sort:
            nsn = nx.single_source_shortest_path_length(G, i, cutoff=2)
            num = 0
            for key,value in nsn.items():
                if value == 2:
                    num = num + 1    
            nsn_list = nsn_list + str(num) + "\t" 
        file.write("Number_second_neighbors\t" + nsn_list + "\n")

        ### Katz centrality ###
        #kc = nx.katz_centrality(G,0.5)
        #print("Katz centrality")
        #for key,value in kc.items():
        #    print(key, ";", value)

        ### Bi-BC ###
        # Get the nodes from /just/ the giant component because the rbc function won't work on networks with multiple subgraphs        
        gc_nodes = gc.nodes() # gc was already made earlier in the mean geodesic pathlength function

        # If the giant component is the same size as the second largest component, then return no values for BiBC.
        # Find the size of the second largest component to see if it is the same size as the largest comp
        subg = sorted(connected_component_subgraphs(G), key = len, reverse = True)

        # Make an empty list and string to add the output to. These will be updated right away if there are multiple giant components, 
        # or will be updated at the end if there is only one gc 
        otu_pheno_value_list = []        
        otu_pheno_value_str = ""
         
        def BiBC(choice, node_file, nodes_in_gc, type1, type2, giantcomp, calc_type, otu_pv_list, otu_pv_str):
            '''
            Calculates RBC/BiBC on the specified nodes 
        
            This code was obtained from yatu's post here: https://stackoverflow.com/questions/53958700/plotting-the-degree-distribution-of-a-graph-using-nx-degree-histogram'
        
            '''   
            
            # If the user wants to calculate based on pre-defined node types
            if choice == "node_types":
                # Pass the gc nodes to the function that will assign the correct node types to each of the nodes    
                otu_pheno_types = assign_node_type(node_file, nodes_in_gc, type1, type2)

                # Calculate rbc using the above function that Kevin wrote, which takes the outputs of assign_node_type
                rbc = restricted_betweenness_centrality(giantcomp, otu_pheno_types['Type1'], otu_pheno_types['Type2'], bibc_calc_type)
                #print(rbc)                

            # Otherwise, if they wish to use modularity as the BiBC parameter...
            elif choice == "modularity":
                nodes_in_gc_for_bibc_mod = max(connected_component_subgraphs(G), key=len)
                nodes_from_bibc_mod = bibc_mod(nodes_in_gc_for_bibc_mod)
                rbc = restricted_betweenness_centrality(gc, nodes_from_bibc_mod['mod1'], nodes_from_bibc_mod['mod2'], bibc_calc_type)
    
            #print(rbc)

            # Combine the rbc function output into one single dictionary
            merged_rbc = {**rbc[0], **rbc[1], **rbc[2]}

            # Loop through the list of sorted node names and for each one create a new listing in bibc_dict_w_NAs that describes
            # 1) whether the node is present in the network or not and 2) what the BiBC is of that node
            bibc_dict_w_NAs = {}
            for i in node_names_sort:
                # Check if node is in giant comp by checking to see if it was output from the rbc function (which only accepts the giant comp as input)
                # If it is present, then just create a key of that node, using the calculated BiBC value
                if i in merged_rbc:
                    bibc_dict_w_NAs[i] = merged_rbc[i]
                # If it is not present, then still create a key of the current node, but assign it NA since it was not in the giant comp. 
                else:
                    bibc_dict_w_NAs[i] = "NA"

            # Order the previously made dictionary by key name so it can be input into the properties file
            ordered_bibc = OrderedDict(sorted(bibc_dict_w_NAs.items()))
                
            for key,value in ordered_bibc.items():
                otu_pv_list.append(value) # otu_pheno_value_list created prior to checking if there are multiple giant components
        
            # Loop through the list containing just the values, inj order of node names and add each to the otu_pheno_value_str
            for i in otu_pv_list:
                otu_pv_str = otu_pv_str + str(i) + "\t"
    
            return(otu_pv_str)        
            
        # The following try-except statement is intended to be used ONLY if there is an issue in calculating rbc, such as there not being 
        # enough nodes (too restrictive of cutoffs), there only being one modular region of the giant component (if the suer specifies 
        # modularity for bibc parameter), etc. It should (ideally) not be used to catch other errors.
        if args.bibc:
            print("Finding each node's BiBC...")

            # First check if there are multiple connected components and if so then...
            if len(subg) > 1:
                print("Multiple components were found in the graph")

                # Check if the two giant comps are the same size. If so then do not bother calculating rbc
                if len(subg[0]) == len(subg[1]):
                    print("There are multiple giant components. BiBC of each node = NA.")
    
                    # Make a list that is the length of the number of nodes in the whole network and write "NA" for the BiBC of each 
                    # node, then add those to a string and write the string to the output file
                    for i in range(len(node_names_sort)):
                        otu_pheno_value_list.append("NA")
                
                    for i in otu_pheno_value_list:
                        otu_pheno_value_str = otu_pheno_value_str + i + "\t"
    
                    file.write("BiBC\t" + otu_pheno_value_str + "\n")
                
                # Otherwise, calculate BiBC      
                elif (len(subg[0]) != len(subg[1])):
                    print("BiBC being calculated for giant component.")
                    output_str = BiBC(bibc_choice, node_input_file, gc_nodes, node_type1, node_type2, gc, bibc_calc_type, otu_pheno_value_list, otu_pheno_value_str)
                    file.write("BiBC\t" + output_str + "\n")


            # If there is only one giant comp, then compute BiBC on whole graph.        
            else:
                print("There is only one component. BiBC being calculated on the entire graph.")
                output_str = BiBC(bibc_choice, node_input_file, gc_nodes, node_type1, node_type2, gc, bibc_calc_type, otu_pheno_value_list, otu_pheno_value_str)
                file.write("BiBC\t" + output_str + "\n")
    
        # Overwrite any previous file with same name instead of appending    
        file.truncate()
                
    file.close()

print("\nNetwork and node properties have been calculated. Check the *properties.txt files in " + filedir + "/.\n")
