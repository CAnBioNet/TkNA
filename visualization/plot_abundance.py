# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 14:35:54 2023

@author: Nolan K Newman <newmanno@oregonstate.edu>
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

        return(all_nodes_dict)

    def plot_abund(dat, grouping, nodename, group_by, colordict):
        '''
        Function that creates abundance barplots from u

        Arguments:
            - dat: a single data frame containinig all abundance values, treatment conditions, and experiment numbers for all nodes across all experiments
            - grouping: Name of the property to plot on the x-axis, obtained from args.x_axis
            - nodename: Name of the parameter to plot on the y-axis
            - group_by: Variable to group the data by (i.e. Treatment or Experiment)
            - colordict: A dictionary with keys of treatment group names and values of the colors
        '''
        sns.set(rc={'figure.figsize':(11.7,8.27)})
        sns.set_style(style='white')

        sns.stripplot(data=concat_df, x=grouping, y=nodename, hue=group_by, dodge=True, color='Black')

        ax = sns.boxplot(data=concat_df, x=grouping, y=nodename, hue=group_by, palette = colordict, fliersize=0)
        handles, labels = ax.get_legend_handles_labels()

        # Find the number of unique groups
        num_groups = len(dat[group_by].unique())

        plt.legend(handles[num_groups:], labels[num_groups:], bbox_to_anchor=(0.02, 0.98), loc='upper left', borderaxespad=0.)

        return(ax)

    def rename_data(name_list):
        '''
        Function that checks a list of file names and if they include paths, it keeps just the file name (not the path) in the list

        Arguments:
            - name_list: the list of names the user supplies in --abund-data
        '''
        new_names = []

        for i in name_list:
            if '/' in i:
                newname = i.split('/')[-1]
                newname = newname[:-4]
                new_names.append(newname)
            else:
                newname = i[:-4]
                new_names.append(newname)

        return(new_names)

    testing = False

    if testing:
        # node_props = "node_properties_modified.txt"
        pick = "inputs_for_downstream_plots.pickle"
        abund_data = ["./fold/expt1_new.csv", "./fold/expt2_new.csv", "./fold/expt3_new.csv", "./fold/expt4_new.csv"]
        metadata = ["./fold/group_map_1.csv", "./fold/group_map_2.csv", "./fold/group_map_3.csv", "./fold/group_map_4.csv"]
        color_group = "Treatment"
        x_axis_abund = "Experiment"
        user_node_list = []

    else:
        print("")

        parser = argparse.ArgumentParser(description="Example command: python plot_abundance.py --pickle inputs_for_downstream_plots.pickle --abund_data Expt1.csv Expt2.csv --metadata Expt1_meta.csv Expt2_meta.csv --color_group Treatment --x_axis_abund Experiment --nodes_to_plot geneABC geneDEF", add_help=False)

        requiredArgGroup  = parser.add_argument_group('Required arguments')
        requiredArgGroup.add_argument("--pickle", type=str, help="inputs_for_downstream_plots.pickle file output by dot_plots.py", required=True)
        requiredArgGroup.add_argument("--abund-data", type=str, nargs = '+', dest="abund_data", help="List of data files containing expressions/abundances, these are the original data files used for intake_data.py", required=True)
        requiredArgGroup.add_argument("--metadata", type=str, nargs = '+', help="List of metadata files containing Experiment/Treatment columns, these are the original metadata files used for intake_data.py", required=True)
        requiredArgGroup.add_argument("--x-axis", type=str, default="Experiment", choices=['Treatment', 'Experiment'], dest="x_axis", help="Variable you wish to group the data by on the x-axis.", required=True)
        requiredArgGroup.add_argument("--group-names", type=str, nargs = '+', dest="group_names",  help="A list of the names of the treatment groups; must match the order of names in --group-colors. For example, if 'Experiment' is chosen for --x-axis then list the names of the experimental groups in --metadata here. And if 'Treatment' is chosen for --x-axis then list the names of the experiment files, without filename extensions, that were used in --abund-data.", required=True)
        requiredArgGroup.add_argument("--group-colors", type=str, nargs = '+', dest="group_colors",  help="A list of the names of the colors to use for each specified group; must match the order of colors in --group-names. Accepted colors can be found at https://matplotlib.org/stable/gallery/color/named_colors.html", required=True)


        optionalArgGroup  = parser.add_argument_group('Optional arguments')
        optionalArgGroup.add_argument("-h", "--help", action="help", help="Show this help message and exit")
        optionalArgGroup.add_argument("--nodes-to-plot", type=str, nargs='+', dest="nodes_to_plot", help="Specify any specific node (or a list of node) to plot the abundance of")


        args = parser.parse_args()
        pick = args.pickle
        abund_data = args.abund_data
        metadata = args.metadata
        x_axis_abund = args.x_axis
        color_group = ''
        if x_axis_abund == 'Experiment':
            color_group = 'Treatment'
        elif x_axis_abund == 'Treatment':
            color_group = 'Experiment'
        else:
            raise Exception("Error: Please select either 'Experiment' or 'Treatment' for how to group your data on the resulting abundance plot.")
        user_node_list = args.nodes_to_plot

    # Rename the input abundance data if it contains a path to a file
    abund_data_renamed = rename_data(abund_data)



    ###################################
    # Plot abundances of top nodes
    ###################################
    # Load the pickle file that already contains the top nodes to plot and their properties
    p = open(pick, "rb")
    pload = pickle.load(p)
    node_props = pload[0]
    top_nodes = pload[1]
    top_nodes_per_type = pload[2] # not used in this script
    plotdir = pload[3] # location of previously created plots directory

    count = 0

    # Final list of dfs that are subset by treatment and experiment
    subset_dfs = []

    # For each abundance data csv and metadata...
    while count < len(abund_data):

        # Subset abundance data for just top nodes
        ab = pd.read_csv(abund_data[count])


        all_node_names = list(ab['ID'])
        met = pd.read_csv(metadata[count], header=None, names = ['sample_name', 'Treatment'])
        treatments = list(met['Treatment'].unique())

        # Transform metadata to dictionary of treatment:sample
        grouped_dict = met.groupby('Treatment').sample_name.apply(list).to_dict()
        #count += 1

        # Create a data frame that contains the abundance/expression of each top parameter per treatment per experiment
        for k,v in grouped_dict.items():
            col_subset = ab[ab.columns & v]

            col_subset = col_subset.T

            col_subset.insert(len(col_subset.columns), 'Treatment', k)
            col_subset.insert(len(col_subset.columns), 'Experiment', abund_data_renamed[count])
            head_col_subset = []
            head_col_subset.extend(all_node_names)
            head_col_subset.append('Treatment')
            head_col_subset.append('Experiment')

            col_subset.set_axis(head_col_subset, axis=1,inplace=True)

            # append this new df to a list. We will then combine all dfs in that list
            subset_dfs.append(col_subset)

        count += 1

    # Combine all elements of the subset_dfs list into one data frame so it can be plotted
    concat_df = pd.concat(subset_dfs)
    concat_df.to_csv(plotdir + 'table_for_manual_abundance_plotting.csv')

    # Generate the color dictionary for each treatment group for plotting
    coldict = {}
    color_index = 0
    while color_index < len(args.group_colors):
        for treatment in args.group_names:
            coldict[treatment] = args.group_colors[color_index]
            color_index += 1


    # if the user has entered any specific nodes to plot, plot them
    if user_node_list:
        for i in user_node_list:

            plot =  plot_abund(concat_df, x_axis_abund, i, color_group, coldict)
            plot.figure.savefig(plotdir + "abundance_of_" + i + "_grouped_by_" + x_axis_abund + ".png", bbox_inches='tight')
            print("Plot saved: " + plotdir + "abundance_of_" + i + "_grouped_by_" + x_axis_abund + ".png")
            plot.figure.clf()

    # Plot all the nodes that were in the zoomed in 'per-type' plots of dot_plots.py. Note that by default, 10 nodes per group 
    # will be plotted, but that changes if there are not 10 nodes in that group in the network, or if the user specified a smaller 
    # number via the --top_num_per_type argument in dot_plots.py
    for i in top_nodes:


        plot =  plot_abund(concat_df, x_axis_abund, i, color_group, coldict)
        plot.figure.savefig(plotdir + "abundance_of_" + i + "_grouped_by_" + x_axis_abund + ".png", bbox_inches='tight')
        print("Plot saved: " + plotdir + "abundance_of_" + i + "_grouped_by_" + x_axis_abund + ".png")
        plot.figure.clf()

