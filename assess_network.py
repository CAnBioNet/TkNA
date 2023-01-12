#-*- coding: utf-8 -*-
"""
Created on Mon Apr 29 12:01:46 2019

Written in Python v3.5.3

@author: Nolan K Newman <newmanno@oregonstate.edu>

Input: the correlations_bw_signif_measurables.csv file output by to_csv.py. This is the a CSV file consisting of correlations 
       between significant nodes (as dictated by the config file)

Output: a pickled file of a networkx graph, as well as a .txt file containing basic properties of the network

"""

###################################################################################
 ##################### import network and make networkx object ###################
###################################################################################

import argparse
import csv
import numpy as np 
import pickle
import networkx as nx
from collections import Counter
from math import factorial

fc_parameters = {}
# nodes_in_final_nw = []


parser = argparse.ArgumentParser(description='Example: python assess_network.py <network file> --source partner1 --target parnter2\n\n')
parser.add_argument("input", help = 'Correlation file (correlations_bw_signif_measurables.csv)')
parser.add_argument("--source", default = 'partner1', help = 'Source node for each edge')
parser.add_argument("--target", default = 'partner2', help = 'Target node for each edge')

args = parser.parse_args()

net_file = args.input
net_file_trimmed = net_file[:-4] # trim the ".csv" or ".txt" from the input file string   


# Counters 
puc_compliant = 0 # number of puc-compliant correlations
puc_noncompliant = 0 # number of puc-compliant correlations
positive_corr_all = 0 # number of positive correlations, regardless of if they are in the final graph
positive_corr_nw = 0 # number of positive correlations, ONLY if they are in the final graph
negative_corr_all = 0 # number of negative correlations, regardless of if they are in the final graph
negative_corr_nw = 0 # number of negative correlations, ONLY if they are in the final graph
pos_nw_edge = 0 # number of positive edges (in the final graph)
neg_nw_edge = 0 # number of positive edges (in the final graph) 
row_count = 0 

G = nx.Graph() 
    
# import specified file into python
with open(net_file) as csvfile:
    file = csv.reader(csvfile, delimiter = ',')
        
    for row in file:
    
        
        # Take the index of the source and target node in the header of the file
        if row_count == 0: 
            p1 = int(row.index(args.source))
            p2 = int(row.index(args.target)) + 1
        
        parameters = row[p1:p2]
        list_to_tuple = tuple(parameters)

        fc_node1_column = len(row) - 5 # column that holds FC direction of source node
        fc_node2_column = len(row) - 4 # column that holds FC direction of target node
        final_nw_edge_column = len(row) - 1 # column that identifies whether edge made it into the final network and its edge value
        
        # Find FC direction of each parameter
        fc_parameters[parameters[0]] = row[fc_node1_column].strip()
        fc_parameters[parameters[1]] = row[fc_node2_column].strip()
        
        # Find the number of total positive and negative corrs, regardless if they made it into the network
        corr_dir_column = len(row) - 6
        
        if row[corr_dir_column] == str(1):
            positive_corr_all += 1
        elif row[corr_dir_column] == str(-1):
            negative_corr_all += 1
            
        # Find each node that made it into the final network and the direction of the edge   
        if row_count != 0:
            
            #print(row[final_nw_edge_column])
            
            if row[final_nw_edge_column] == str(1.0):
                pos_nw_edge += 1
                positive_corr_nw += 1
                G.add_edge(parameters[0], parameters[1])
            elif row[final_nw_edge_column] == str(-1.0):
                neg_nw_edge += 1 
                negative_corr_nw += 1
                G.add_edge(parameters[0], parameters[1])
            
      
        # Is each edge PUC-compliant?
        puc_col = len(row) - 2
        
        #print(row[puc_col].strip())
        
        if row[puc_col].strip() == str(1):
            puc_compliant += 1
        elif row[puc_col].strip() == str(0):
            puc_noncompliant += 1
            
            
        row_count += 1

    csvfile.close()

del fc_parameters[args.source]
del fc_parameters[args.target]

# function to get unique values
def unique(list1):
    x = np.array(list1)
    return np.unique(x)

#print(fc_parameters)

# Subset dictionary holding FC of all params for just nodes in nw

#print(G.nodes())
    
fc_parameters_final_nw = dict((k, fc_parameters[k]) for k in G.nodes())   

# for k in G.nodes():
#     print(k)
#     print(fc_parameters[k])
    
# Find the ratio of positive:negative edges (and vice versa) in the observed graph
if int(negative_corr_nw) != 0:
    obs_posneg_ratio = int(positive_corr_nw)/int(negative_corr_nw)
else:
    obs_posneg_ratio = 1.0
    
if int(positive_corr_nw) != 0:
    obs_negpos_ratio = int(negative_corr_nw)/int(positive_corr_nw)
else:
    obs_negpos_ratio = 1.0


# Count the number of positive and negative nodes in the network  
nodedir = Counter(fc_parameters_final_nw.values())
pos_nodes = nodedir['1.0']
neg_nodes = nodedir['-1.0']
        
#print(pos_nodes)
#print(neg_nodes)
total_nodes = int(pos_nodes) + int(neg_nodes)

obs_edge_node_ratio = G.number_of_edges() / total_nodes

obs_posneg_node_ratio = int(pos_nodes) / int(neg_nodes)
obs_negpos_node_ratio = int(neg_nodes) / int(pos_nodes)



# Find the number of edges in a full graph  
expec_pos = int(factorial(pos_nodes)/(2 * factorial(pos_nodes - 2)) + factorial(neg_nodes)/(2 * factorial(neg_nodes - 2)))               
expec_neg = pos_nodes * neg_nodes
expec_total = expec_pos + expec_neg          
expec_edge_node_ratio = expec_total / total_nodes

# Find the ratio of positive:negative edges (and vice versa) in a full graph
ideal_ratio_posneg = round(expec_pos/expec_neg, 2)
ideal_ratio_negpos = round(expec_neg/expec_pos, 2)

# Calculate the non-normalized deviation from the expected (full) graph
dev_posneg = round(obs_posneg_ratio/ideal_ratio_posneg, 2)
dev_negpos = round(obs_negpos_ratio/ideal_ratio_negpos, 2)

#Calculate the normalized deviation from the expected (full) graph
dev_norm_posneg = round((obs_posneg_ratio - ideal_ratio_posneg) / ideal_ratio_posneg, 2)
dev_norm_negpos = round((obs_negpos_ratio - ideal_ratio_negpos) / ideal_ratio_negpos, 2)

# calculate the normalized deviation of the edge:node (density) from the full graph
dens_dev = (abs(obs_edge_node_ratio - expec_edge_node_ratio)) / expec_edge_node_ratio



# Calculate PUC (the proportion of edges that do not follow the expected direction). Remember row_count is number of rows in input file
puc = round((100 * (puc_noncompliant / (row_count - 1))), 2)

# mean degree
mdeg = 2 * G.number_of_edges() / G.number_of_nodes()

# Edges in network over the edges in a full graph
nw_edge_full_graph_ratio = round(G.number_of_edges()/expec_total, 2)



print("\nPUC: " + "{}%".format(str(puc)))
print("Mean degree: " + str(mdeg))
print("Edges in network over the edges in a full graph: " + str(nw_edge_full_graph_ratio))
print("Positive:negative edge ratio: " + str(obs_posneg_ratio) + "\n")

print("Expected positive:negative edge ratio in full graph: " + str(ideal_ratio_posneg))
print("Deviation from expected positive:negative edge ratio in full graph: " + str(dev_norm_posneg) + "\n")

print("Number of positive nodes: " + str(pos_nodes))
print("Number of negative nodes: " + str(neg_nodes))
print("Number of total nodes: " + str(G.number_of_nodes()) + "\n")

print("Number of positive edges: " + str(pos_nw_edge))
print("Number of negative edges: " + str(neg_nw_edge))
print("Number of total edges: " + str(G.number_of_edges()) + "\n")


with open("network_quality_assessment.txt", "w") as file:
    file.write("PUC: " + "{}%".format(str(puc)) + "\n")
    file.write("Mean degree: " + str(mdeg) + "\n")
    file.write("Edges in network over the edges in a full graph: " + str(nw_edge_full_graph_ratio) + "\n") 
    file.write("Positive:negative edge ratio: " + str(obs_posneg_ratio) + "\n\n")
    
    file.write("Expected positive:negative edge ratio in full graph: " + str(ideal_ratio_posneg) + "\n")
    file.write("Deviation from expected positive:negative edge ratio in full graph: " + str(dev_norm_posneg) + "\n\n")
    
    file.write("Number of positive nodes: " + str(pos_nodes) + "\n")
    file.write("Number of negative nodes: " + str(neg_nodes) + "\n")
    file.write("Number of total nodes: " + str(G.number_of_nodes()) + "\n\n")
    
    file.write("Number of positive edges: " + str(pos_nw_edge) + "\n")
    file.write("Number of negative edges: " + str(neg_nw_edge) + "\n")
    file.write("Number of total edges: " + str(G.number_of_edges()) + "\n")


    file.truncate()

file.close()

pickle.dump(G, open("network.pickle", "wb"))

# pickle_list = [G, dev_dict] # output both tht network and deviation dictionary
# pickle.dump(pickle_list, pickle_out)
# pickle_out.close()

print("Assessment calculations have been completed. Results have been saved in network_quality_assessment.txt.\n")











