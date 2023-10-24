#-*- coding: utf-8 -*-
"""
Created on Mon Apr 29 12:01:46 2019

Last updated: 7/19/23

Written/tested in Python v3.8.10

@author: Nolan K Newman <newmanno@oregonstate.edu>

Usage
	python ./analysis/assess_network.py --file <network file> --out-dir <directory>
Example command
	python ./ analysis/assess_network.py --file ./project_folder/output/network_output/correlations_bw_signif_measurables.csv --out-dir ./project_folder/output/network_output/
Inputs
	--file: correlations_bw_signif_measurables.csv file created with to_csv.py
	--out-dir: Path to the directory to output results to
Outputs
    network_quality_assessment.csv: Contains the network quality statistics of the reconstructed network, calculated per edge type. 

"""

###################################################################################
 ##################### import network and make networkx object ###################
###################################################################################

import argparse
import pandas as pd
import csv
import numpy as np 
import networkx as nx
from collections import Counter
from math import factorial
import os
from itertools import combinations_with_replacement

fc_parameters = {}

testing = False

if testing:
    net_file = "correlations_bw_signif_measurables.csv"
    outdir = '.'
    
else:

    parser = argparse.ArgumentParser(description="Example command: python assess_network.py --file <correlation file> --out-dir <directory>", add_help=False)
    requiredArgGroup = parser.add_argument_group('Required arguments')        
    requiredArgGroup.add_argument("--file", type=str, help="correlations_bw_signif_measurables.csv file output by to_csv.py", required=True)
    requiredArgGroup.add_argument("--out-dir", type=str, dest = "outdir", help="Path to output directory", required=True)
    
    optionalArgGroup = parser.add_argument_group('Optional arguments')   
    optionalArgGroup.add_argument("-h", "--help", action="help", help="Show this help message and exit")
    
    args = parser.parse_args()
    
    net_file = args.file
    outdir = args.outdir
    
def find_unique_networks(df):
    '''
    Function to find redundant pairs of edge types (e.g. microbe<==>gene and gene<==>microbe) and extract subnetworks based on those unique types

     Arguments:
         - df: data frame with all pairs of data types in the "Edge_Type column"

    '''    
    # Find all unique pairs of data types in input file
    all_params = []
    uniq_pairs = df.Edge_Type.unique()
    
    # Get a list of all parameters
    for i in uniq_pairs:  
        all_params.append(i.split('<==>')[0])
        all_params.append(i.split('<==>')[1])
        
    # Find all unique pairs of parameters
    uniq_params = [*set(all_params)]
    set_possible_combinations = set(combinations_with_replacement(uniq_params, r=2))
    
    # Convert the set to a list
    all_possible_combinations = []
    for val in set_possible_combinations:
        all_possible_combinations.append(str(val[0]) + "<==>" + str(val[1]))
    
    # Rename the edge in the data frame by flipping parameter 1 and 2 if the edge is not found in the all_possible_combinations list. 
    for index, row in df.iterrows():
        if row['Edge_Type'] not in all_possible_combinations:
            new_edge_type = row['Edge_Type'].split('<==>')
            new_edge_type = new_edge_type[1] + '<==>' + new_edge_type[0]
            df.at[index, 'Edge_Type'] = new_edge_type

    #Take all unique network types and save in list
    nws_by_type = []
    for i in all_possible_combinations:
        subnw = df.loc[df['Edge_Type'] == i]
        nws_by_type.append(subnw)
        
    return(nws_by_type)
        
        
net_file_trimmed = net_file[:-4] # trim the ".csv" or ".txt" from the input file string   
    
corr_file = pd.read_csv(net_file)
corr_file.columns = corr_file.columns.str.replace(' ', '_')
corr_file.columns = [*corr_file.columns[:-1], 'Final_Network_Value']

corr_file['combined_Coefficient_correlation_Direction'] = corr_file['combined_Coefficient_correlation_Direction'].astype('int')
corr_file['partner1_FC_direction'] = corr_file['partner1_FC_direction'].astype('int')
corr_file['partner2_FC_direction'] = corr_file['partner2_FC_direction'].astype('int')
corr_file['IfFoldChangeDirectionMatch'] = corr_file['IfFoldChangeDirectionMatch'].astype('int')
corr_file['PUC'] = corr_file['PUC'].astype('int')
corr_file['Final_Network_Value'] = corr_file['Final_Network_Value'].astype('int')


nws_by_type = find_unique_networks(corr_file)

output_properties_df = pd.DataFrame()
    
#print(nws_by_type[0])

for nw in nws_by_type:
    
    nw_type = str(nw.Edge_Type.unique()[0])

    fc_parameters = {}
    G = nx.Graph() 
    output_list = []
    
    # Counters 
    puc_noncompliant = 0 # number of puc-noncompliant correlations
    positive_corr_all = 0 # number of positive correlations, regardless of if they are in the final graph
    positive_corr_nw = 0 # number of positive correlations, ONLY if they are in the final graph
    negative_corr_all = 0 # number of negative correlations, regardless of if they are in the final graph
    negative_corr_nw = 0 # number of negative correlations, ONLY if they are in the final graph
    pos_nw_edge = 0 # number of positive edges (in the final graph)
    neg_nw_edge = 0 # number of positive edges (in the final graph) 
    row_count = 0 
    signif_meta_edge = 0 # number of edges that are significant and pass the meta-analysis thresholds
    
    for index, row in nw.iterrows():

        parameters = row['partner1'], row['partner2']
    
        list_to_tuple = tuple(parameters)

        # Find FC direction of each parameter
        fc_parameters[parameters[0]] = row['partner1_FC_direction']
        fc_parameters[parameters[1]] = row['partner2_FC_direction']
    
        # Find the number of total positive and negative corrs, regardless if they made it into the network
        if row['combined_Coefficient_correlation_Direction'] == 1:
            positive_corr_all += 1
        elif row['combined_Coefficient_correlation_Direction'] == -1:
            negative_corr_all += 1
        
        if row['Final_Network_Value'] == 1:
            pos_nw_edge += 1
            positive_corr_nw += 1
            G.add_edge(parameters[0], parameters[1])
        elif row['Final_Network_Value'] == -1:
            neg_nw_edge += 1 
            negative_corr_nw += 1
            G.add_edge(parameters[0], parameters[1])
    
        # Count the number of significant edges that pass meta-analysis thresholds for the
        # denominator of the PUC calculation
        if str(row['All_Non-PUC_Filters_Passed']) == "True":
            signif_meta_edge += 1
     
            # For just those edges that are significant and pass those thresholds, also 
            # check if they are an expected edge for the numerator of the PUC calculation           
            rowpuc = row['PUC'] 
            if row['PUC'] == 0:
                puc_noncompliant += 1             


    # Subset dictionary holding FC of all params for just nodes in current nw    
    fc_parameters_current_nw = dict((k, fc_parameters[k]) for k in G.nodes())   

    # Find the ratio of positive:negative edges (and vice versa) in the observed graph
    if int(negative_corr_nw) != 0:
        obs_posneg_ratio = int(positive_corr_nw)/int(negative_corr_nw)
    else:
        obs_posneg_ratio = "Inf"
        
    if int(positive_corr_nw) != 0:
        obs_negpos_ratio = int(negative_corr_nw)/int(positive_corr_nw)
    else:
        obs_negpos_ratio = "Inf"


    # Count the number of positive and negative nodes in the current network  
    nodedir = Counter(fc_parameters_current_nw.values())
    pos_nodes = nodedir[1]
    neg_nodes = nodedir[-1]
    
    total_nodes = int(pos_nodes) + int(neg_nodes)

    if total_nodes != 0:
        obs_edge_node_ratio = G.number_of_edges() / total_nodes
    else: 
        obs_edge_node_ratio = "NaN"

    if neg_nodes != 0:
        obs_posneg_node_ratio = int(pos_nodes) / int(neg_nodes)
    else:
        obs_posneg_node_ratio = "NaN"
        
    if pos_nodes != 0:            
        obs_negpos_node_ratio = int(neg_nodes) / int(pos_nodes)
    else:
        obs_negpos_node_ratio = "NaN"

    print("\n\n\nSubnetwork type: " + nw_type)

    # Find the number of edges in a full graph
    if pos_nodes > 2 and neg_nodes > 2:
        
        # Calculate PUC (the proportion of edges that do not follow the expected direction). 
        puc = round((100 * (puc_noncompliant / signif_meta_edge)), 2)
        
        # mean degree
        mdeg = 2 * G.number_of_edges() / G.number_of_nodes()
                
        expec_pos = int(factorial(pos_nodes)/(2 * factorial(pos_nodes - 2)) + factorial(neg_nodes)/(2 * factorial(neg_nodes - 2)))               
        expec_neg = pos_nodes * neg_nodes
        expec_total = expec_pos + expec_neg          
        expec_edge_node_ratio = expec_total / total_nodes
            
        # Find the ratio of positive:negative edges (and vice versa) in a full graph
        if expec_neg != 0:
            ideal_ratio_posneg = round(expec_pos/expec_neg, 2)
        else: 
            ideal_ratio_posneg = "Inf"
        
        if expec_pos != 0:
            ideal_ratio_negpos = round(expec_neg/expec_pos, 2)
        else:
            ideal_ratio_negpos = "Inf"
            
        # Calculate the non-normalized deviation from the expected (full) graph
        if obs_posneg_ratio != "Inf" and ideal_ratio_posneg != "Inf":
            dev_posneg = round(obs_posneg_ratio/ideal_ratio_posneg, 2)
        else:
            dev_posneg = "NA"
            
        if obs_negpos_ratio != "Inf" and ideal_ratio_negpos != "Inf":
            dev_negpos = round(obs_negpos_ratio/ideal_ratio_negpos, 2)
        else:
            dev_negpos = "NA"
        
        #Calculate the normalized deviation from the expected (full) graph
        if obs_posneg_ratio != "Inf" and ideal_ratio_posneg != "Inf":
            dev_norm_posneg = round((obs_posneg_ratio - ideal_ratio_posneg) / ideal_ratio_posneg, 2)
        else:
            dev_norm_posneg = "NA"
            
        if obs_negpos_ratio != "Inf" and ideal_ratio_negpos != "Inf":
            dev_norm_negpos = round((obs_negpos_ratio - ideal_ratio_negpos) / ideal_ratio_negpos, 2)
        else:
            dev_norm_negpos = "NA"
        
        # calculate the normalized deviation of the edge:node (density) from the full graph
        dens_dev = (abs(obs_edge_node_ratio - expec_edge_node_ratio)) / expec_edge_node_ratio
    
        # Edges in network over the edges in a full graph
        nw_edge_full_graph_ratio = round((G.number_of_edges()/expec_total) * 100, 2)
    
    else:
        print("\n\n!!!\nWarning: Not enough positive or negative nodes to calculate some or all basic properties of the " + nw_type + " network. \nThese calculations require the input network to have at least two positive log2 foldchange nodes and two negative log2 foldchange nodes.\n!!!")
        
        if signif_meta_edge != 0:
            # Calculate PUC (the proportion of edges that do not follow the expected direction). 
            puc = round((100 * (puc_noncompliant / signif_meta_edge)), 2)
        else:
            puc = "NA"  
            
        # mean degree
        if G.number_of_nodes() != 0:
            mdeg = 2 * G.number_of_edges() / G.number_of_nodes()      
        else:
            mdeg = "NaA"
            
        expec_pos = "NA"
        expec_neg = "NA"
        expec_total = "NA"          
        expec_edge_node_ratio = "NA"
    
        # Find the ratio of positive:negative edges (and vice versa) in a full graph
        ideal_ratio_posneg = "NA"
        ideal_ratio_negpos = "NA"
            
        # Calculate the non-normalized deviation from the expected (full) graph
        dev_posneg = "NA"
        dev_negpos = "NA"
        
        #Calculate the normalized deviation from the expected (full) graph
        dev_norm_posneg = "NA"
        dev_norm_negpos = "NA"
        
        # calculate the normalized deviation of the edge:node (density) from the full graph
        dens_dev = "NA"
    
        # Edges in network over the edges in a full graph
        nw_edge_full_graph_ratio = "NA"
        
  
    output_list.append(puc)
    output_list.append(mdeg)     
    output_list.append(nw_edge_full_graph_ratio)
    output_list.append(obs_posneg_ratio)
    output_list.append(ideal_ratio_posneg)
    output_list.append(dev_norm_posneg)
    output_list.append(pos_nodes)
    output_list.append(neg_nodes)
    output_list.append(total_nodes)
    output_list.append(positive_corr_nw)
    output_list.append(negative_corr_nw)
    output_list.append(G.number_of_edges())

    output_properties_df[nw_type] = output_list
        

    
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

#Set column names for output file
output_properties_df.index = ['PUC (%)', 
                              'Mean degree', 
                              'Edges in network over edges in full graph (%)', 
                              'Observed positive:negative edge ratio', 
                              'Expected positive:negative edge ratio', 
                              'Deviation from expected positive:negative edge ratio in full graph', 
                              'Positive nodes', 
                              'Negative nodes', 
                              'Total nodes', 
                              'Positive edges', 
                              'Negative edges', 
                              'Total edges']
output_properties_df.index.name = 'Property'

# Correct the path if needed to the output dir
if not outdir[-1] == "/":
    outdir = outdir + "/"
    
if not os.path.exists(outdir):
    os.makedirs(outdir)  

output_properties_df.to_csv(f"{outdir}network_quality_assessment.csv", na_rep = "NA")

print("Assessment calculations have been completed. Results have been saved in " + outdir + "network_quality_assessment.csv\n")



