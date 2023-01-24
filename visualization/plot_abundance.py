# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 14:35:54 2023

@author: Nolan
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 09:23:46 2023

@author: Nolan
"""

import argparse
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import pickle
import seaborn as sns
import math
from matplotlib.legend_handler import HandlerTuple
import numpy as np
import os


if __name__ == '__main__':

    def assign_type(nw):
        '''
        Function that generates a dictionary of node types from the user-supplied network_file
        
        Arguments:
            - nw: panda data frame
        '''
        
        type_dict_p1 = {}
        type_dict_p2 = {}
        all_nodes_dict = {}
        
        nw.columns = [c.replace(' ', '_') for c in nw.columns] # rename columns to not contain spaces
        
        nw[['type1', 'type2']] = nw['Edge_Type'].str.split("<==>", 1, expand = True) # split Edge_Type column to two columns
        
        # Make dictionaries for the type of each partner, then combine into one dict
        type_dict_p1 = nw.set_index('partner1').to_dict()['type1']
        type_dict_p2 = nw.set_index('partner2').to_dict()['type2']
        all_nodes_dict = {**type_dict_p1,**type_dict_p2}
        
        # print(all_nodes_dict)                
        
        return(all_nodes_dict)
    
    def plot_abund(dat, x, y, group_by, pal):
        '''        
        Function that creates abundance barplots from u
        
        Arguments:
            - dat: a single data frame containinig all abundance values, treatment conditions, and experiment numbers for all nodes across all experiments
            - x: Name of the property to plot on the x-axis, obtained from args.x_axis_abund
            - y: Name of the parameter to plot on the y-axis
            - group_by: Variable to group the data by (i.e. Treatment or Experiment)
        '''
        sns.set(rc={'figure.figsize':(11.7,8.27)})
        sns.set_style(style='white') 
        
        sns.stripplot(data=concat_df, x=x, y=y, hue=group_by, dodge=True, color='Black')

        
        ax = sns.boxplot(data=concat_df, x=x, y=y, hue=group_by, palette = pal)    
        handles, labels = ax.get_legend_handles_labels()
        # ax.legend(handles=handles,
        #     labels=treatments,
        #     loc='upper right', handlelength=4,
        #     colorhandler_map={tuple: HandlerTuple(ndivide=None)})   
            
        # Find the number of unique groups
        num_groups = len(dat[group_by].unique())
                
        plt.legend(handles[0:num_groups], labels[0:num_groups], bbox_to_anchor=(0.02, 0.98), loc='upper left', borderaxespad=0.)

        #ax.legend_.remove()
        # handles, labels = ax.get_legend_handles_labels()
        # ax.legend(handles=handles,
        #           labels=treatments,
        #           loc='upper right', handlelength=4,
        #           handler_map={tuple: HandlerTuple(ndivide=None)})   
           
        return(ax)


    testing = False

    if testing:
        # node_props = "node_properties_modified.txt"
        pick = "inputs_for_downstream_plots.pickle"
        # network_file = "network_output_comp.csv"
        # top_abund_prop = "BiBC"
        # node_type = "pheno"
        abund_data = ["Expt1_new.csv", "Expt2_new.csv", "Expt3_new.csv", "Expt4_new.csv"]
        metadata = ["group_map_1.csv", "group_map_2.csv", "group_map_3.csv", "group_map_4.csv"]
        color_group = "Treatment"
        x_axis_abund = "Experiment"
        user_node_list = []
    
    else:
        print("")

        parser = argparse.ArgumentParser(description="Example command: python plot_abundance.py --pickle inputs_for_downstream_plots.pickle --abund_data Expt1.csv Expt2.csv --metadata Expt1_meta.csv Expt2_meta.csv --color_group Treatment --x_axis_abund Experiment --nodes_to_plot geneABC geneDEF", add_help=False)

        requiredArgGroup  = parser.add_argument_group('Required arguments')        
        requiredArgGroup.add_argument("--pickle", type=str, help="inputs_for_downstream_plots.pickle file output by dot_plots.py", required=True)
        # requiredArgGroup.add_argument("--network_file", type=str, help="network_output_comp.csv file output by to_csv.py", required=True)
        requiredArgGroup.add_argument("--abund_data", type=str, nargs = '+', help="List of data files containing expressions/abundances, these are the original data files used for intake_data.py", required=True)
        requiredArgGroup.add_argument("--metadata", type=str, nargs = '+', help="List of metadata files containing Experiment/Treatment columns, these are the original metadata files used for intake_data.py", required=True)
        requiredArgGroup.add_argument("--color_group", type=str, default='Treatment', choices=['Treatment', 'Experiment'], help="Variable to color the plot by", required=True)
        requiredArgGroup.add_argument("--x_axis", type=str, default="Experiment", choices=['Treatment', 'Experiment'], help="Variable you wish to group the data by on the x-axis", required=True)  
        
        optionalArgGroup  = parser.add_argument_group('Optional arguments')        
        optionalArgGroup.add_argument("-h", "--help", action="help", help="Show this help message and exit")
        optionalArgGroup.add_argument("--nodes_to_plot", type=str, nargs='+', help="Specify any specific node (or a list of node) to plot the abundance of")


        # parser = argparse.ArgumentParser(description="Example command: python plot_abundance.py inputs_for_downstream_plots.pickle --abund_data Expt1.csv Expt2.csv --metadata Expt1_meta.csv Expt2_meta.csv --color_group Treatment --x_axis_abund Experiment --nodes_to_plot geneABC geneDEF")
        # # Args for abundance/expression plots
        # parser.add_argument("pickle", type=str, help="inputs_for_downstream_plots.pickle file output by dot_plots.py")
        # parser.add_argument("network_file", type=str, help="network_output_comp.csv file output by to_csv.py")
        # parser.add_argument("--abund_data", type=str, nargs = '+', help="List of data files containing expresions/abundances, these are the original data files used for intake_data.py", required=True)
        # parser.add_argument("--metadata", type=str, nargs = '+', help="List of metadata files containing Experiment/Treatment columns, these are the original metadata files used for intake_data.py", required=True)
        # # parser.add_argument("--top_abund_prop", type=str, help="Property from node_properties.txt to select top nodes for abundance plots.  Name must match property name in node_properties.txt")  
        # # parser.add_argument("--node_type", type=str, help="The node type (from network_output_comp.csv) to subset for before finding top nodes")
        # parser.add_argument("--color_group", type=str, default='Treatment', choices=['Treatment', 'Experiment'], help="Variable to color the plot by", required=True)        
        # parser.add_argument("--x_axis_abund", type=str, default="Experiment", choices=['Treatment', 'Experiment'], help="Variable you wish to group the data by", required=True)
        # parser.add_argument("--nodes_to_plot", type=str, nargs='+', help="Optional; User can specify any specific node (or a list of node) they wish to plot the abundance of")

    
        args = parser.parse_args()    
        pick = args.pickle
        # network_file = args.network_file
        abund_data = args.abund_data
        metadata = args.metadata
        # node_type = args.node_type
        color_group = args.color_group
        x_axis_abund = args.x_axis
        user_node_list = args.nodes_to_plot

    # nw_df = pd.read_csv(network_file)
    # nw_w_types = assign_type(nw_df)




    ###################################
    # Plot abundances of top nodes
    ###################################
    # Load the pickle file that already contains the top nodes to plot and their properties
    p = open(pick, "rb")
    pload = pickle.load(p)
    node_props = pload[0]
    top_nodes = pload[1]
    top_nodes_per_type = pload[2] # not used in this script


    # # Loop through each unique node type and extract just the top nodes for each type
    # # subset df for just the current data type
    # sub = node_props[node_props['Data_type']  == node_type]

    # # Find the top nodes    
    # try: 
    #     sub = sub.sort_values(by=top_abund_prop, ascending=False)
    #     sub_top = sub.head(top_abund_num)
    #     top_nodes = list(sub_top['index'])
        
    # except:
    #     raise Exception("Could not create plot for " + node_type + " data type based on top " + str(top_abund_num) + " " + top_abund_prop + " nodes. Please ensure ")
        
    count = 0

    # Final list of dfs that are subset by treatment and experiment
    subset_dfs = []
        
    # For each abundance data csv and metadata...
    while count < len(abund_data):
        #print(abund_data[count])
        
        # Subset abundance data for just top nodes
        ab = pd.read_csv(abund_data[count])
        all_node_names = list(ab['ID'])

        # print("BEFORE")
        # print(ab_sub)
        # count += 1 
        
        # ab['ID'] =  pd.Categorical(ab.ID, categories = list(ab['ID']), ordered = True)
        # print('AFTER')
        # ab = ab.sort_values(by='ID')
        # print(ab_sub)
        # count += 1 

        met = pd.read_csv(metadata[count], header=None, names = ['sample_name', 'Treatment'])
        treatments = list(met['Treatment'].unique())

        # Transform metadata to dictionary of treatment:sample 
        grouped_dict = met.groupby('Treatment').sample_name.apply(list).to_dict()
            
        # Create a data frame that contains the abundance/expression of each top parameter per treatment per experiment
        for k,v in grouped_dict.items():
            # print(k,v)
            col_subset = ab[ab.columns & v]
        
            col_subset = col_subset.T

            col_subset.insert(len(col_subset.columns), 'Treatment', k)
            col_subset.insert(len(col_subset.columns), 'Experiment', abund_data[count])
            head_col_subset = []
            head_col_subset.extend(all_node_names)
            head_col_subset.append('Treatment')
            head_col_subset.append('Experiment')
            # print(head_col_subset)

            # head_col_subset = top_nodes
            col_subset.set_axis(head_col_subset, axis=1,inplace=True)
            # print(col_subset)

            # append this new df to a list. We will then combine all dfs in that list
            subset_dfs.append(col_subset)           
        
        count += 1 

    # Combine all elements of the subset_dfs list into one data frame so it can be plotted
    concat_df = pd.concat(subset_dfs)
    concat_df.to_csv('./plots/table_for_manual_abundance_plotting.csv')
    
    
    # if the user has entered any specific nodes to plot, plot them
    if user_node_list:
        for i in user_node_list:
            # print(i)
            plot =  plot_abund(concat_df, x_axis_abund, i, color_group, 'Greys')
            plot.figure.savefig("./plots/abundance_of_" + i + "_grouped_by_" + color_group + ".png", bbox_inches='tight') 
            print("Plot saved: abundance_of_" + i + "_grouped_by_" + color_group + ".png")
            plot.figure.clf()
    
    # Plot all the nodes that were in the zoomed in 'per-type' plots of dot_plots.py. Note that by default, 10 nodes per group will be plotted, but that changes if there are not 10 nodes in that group in the network, or if the user specified a smaller number via the --top_num_per_type argument in dot_plots.py
    for i in top_nodes:
        # print(i)
        #sns.set_style(style='white') 
        #sns.set_context("notebook")
        
        plot =  plot_abund(concat_df, x_axis_abund, i, color_group, 'Greys')
        plot.figure.savefig("./plots/abundance_of_" + i + "_grouped_by_" + color_group + ".png", bbox_inches='tight') 
        print("Plot saved: abundance_of_" + i + "_grouped_by_" + color_group + ".png")
        plot.figure.clf()


    
