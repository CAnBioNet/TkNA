# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 11:52:11 2022


1. Reads in the original data, network created, and the properties file that was created
2. Creates the following graphs
    a. Degree distribution
    b. Degree vs. BiBC
    c. Abundance of top degree/BiBC nodes
    d. Degree-BiBC density distribution

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


pd.options.mode.chained_assignment = None  # default='warn'


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

    # subset a df for just the top X%
    def subset_for_top(df, nrow):
        '''
        Function that takes a data frame and extracts just the top n rows

        Arguments:
            - df: panda data frame
            - nrow: number of rows
        '''

        df = df.head(nrow)
        return(df)

    def plot_scatter(choicex, choicey, df, h, col):
        '''
        Function that creates a seaborn scatterplot

        Arguments:
            - choicex: Name of the property to plot on the x-axis, obtained from args.propx
            - choicey: Name of the property to plot on the y-axis, obtained from args.propy
            - df: data frame to use for plotting
            - h: hue to color nodes by
            - col: color/color palette to use for colors of dots
        '''

        degbibc = sns.scatterplot(data=df, x=df[choicex], y=df[choicey], hue = h, s=200, color = col)

        degbibc.set(xlabel = choicex, ylabel = choicey)

        degbibc.set(ylim=((min(df[choicey]) - 0.1 * min(df[choicey])), max(df[choicey])+(0.15 * (max(df[choicey])-min(df[choicey]))))) # add 0.1 for buffer at top
        degbibc.figure.set_size_inches(8, 6)

        return(degbibc)

    # Add dots on top of density plot
    def add_text(df, xchoice, ychoice, fig):
        '''
        Function that creates a dictionary keyed on node name, with a value of a list of the values calculated for the two user-selected properties, default x is degree, y is BiBC

        Arguments:
            - df: data frame pre-filtered to contain desired nodes in 'index' column containing node names and other columns containing the x and y properties
            - xchoice: user-defined property to plot on the X-axis
            - ychoice: user-defined property to plot on the Y-axis
            - fig: seaborn scatterplot figure to add text labels to
        '''
        name_prop_prop_dict = {}

        for index, row in df.iterrows():
            name_prop_prop_dict[row[0]] = [row[xchoice], row[ychoice]]

        for k,v in name_prop_prop_dict.items():
            fig.annotate(str(k), (v[0],v[1]), fontsize=18)

        return(fig)

    testing = False

    if testing:
        pickle_file = "network.pickle"
        node_props = "node_properties_modified.txt"
        network_file = "correlations_bw_signif_measurables.csv"
        propx = "BiBC"
        propy = "Node_closeness"
        top_num = 5
        top_pct = False
        top_pct_num = 40
        top_num_per_type = 4

    else:
        print("")

        parser = argparse.ArgumentParser(description="Example command: python dot_plots.py --pickle network.pickle --node_props node_properties_modified.txt --network_file network_output_comp.csv --propx BiBC --propy Node_degrees --top_num 5 --top_num_per_type 3", add_help=False)

        requiredArgGroup  = parser.add_argument_group('Required arguments')
        requiredArgGroup.add_argument("--pickle", type=str, help="Pickled network file output by calc_network_properties.py", required=True)
        requiredArgGroup.add_argument("--node-props", dest="node_props", type=str, help="node_properties.txt file output by calc_network_properties.py", required=True)
        requiredArgGroup.add_argument("--network-file", type=str, dest="network_file", help="network_output_comp.csv file output by to_csv.py", required=True)
        requiredArgGroup.add_argument("--propx", default="Node_Degrees", type=str, help="Node property to plot on X-axis. Name must match property name in node_properties.txt", required=True)
        requiredArgGroup.add_argument("--propy", default="BiBC", type=str, help="Node property to plot on Y-axis. Name must match property name in node_properties.txt", required=True)
        requiredArgGroup.add_argument("--top-num", default=10, type=int, dest="top_num", help="Number of nodes you want to zoom in to on the property v property plot", required=True)
        requiredArgGroup.add_argument("--top-num-per-type", type=int, dest="top_num_per_type", help="The number of nodes to plot for each data type when zoomed in on the dot plot", required=True)
        requiredArgGroup.add_argument("--out-dir", dest="outdir", help="The number of nodes to plot for each data type when zoomed in on the dot plot", required=True)

        optionalArgGroup  = parser.add_argument_group('Optional arguments')
        optionalArgGroup.add_argument("-h", "--help", action="help", help="Show this help message and exit")
        optionalArgGroup.add_argument("--top-pct", action = 'store_true', dest="top_pct", help="Flag; Do you want to plot the top X percent of nodes instead of specifying a number with --top_num? Will calculate which nodes to zoom in on based on propx argument. Requires --top_pct_num")
        optionalArgGroup.add_argument("--top-pct-num", type=int, dest="top_pct_num", help="Percent (in integer format) of top nodes to plot; requires --top_pct")


        args = parser.parse_args()

        pickle_file = args.pickle
        node_props = args.node_props
        network_file = args.network_file
        # rand_net = args.rand_net
        propx = args.propx
        propy = args.propy
        top_num = args.top_num
        top_pct_num = args.top_pct_num
        top_pct = args.top_pct
        top_num_per_type = args.top_num_per_type

    outdir = args.outdir
    
    # Correct the path if needed to the output dir
    if not outdir[-1] == "/":
        outdir = outdir + "/"
        
    if not os.path.exists(outdir):
        os.makedirs(outdir)  

    # import the file calc_network_properties creates and create a dictionary to store info into
    with open(node_props, newline='') as txtfile:

        node_props = pd.read_table(txtfile, sep = "\t", index_col=[0])
        node_props = node_props.T
        node_props = node_props[:-1] # Delete last row of df since this is an artifact of adding tabs after every calculation in calc_nw_properties.py
        # print(node_props)

    ###############
    #Plotting
    ###############

    # Unpack the pickle
    p = open(pickle_file, "rb")
    p = pickle.load(p)
    G = p


    #############################################
    # Plot degree distribution
    #############################################
    sns.set_style("white")

    degree_freq = nx.degree_histogram(G)
    degrees = range(len(degree_freq))
    plt.figure(figsize=(12, 8))
    plt.loglog(degrees, degree_freq,'go-', linestyle='None', markersize=10, color='black')
    plt.xlabel('Degree', fontsize=20)
    plt.ylabel('Frequency', fontsize=20)
    plt.savefig(outdir + "degree_distribution_dotplot.png", bbox_inches='tight')
    print("Plot saved:" + outdir + "degree_distribution_dotplot.png")

    plt.clf()


    #####################################################################
    # Plot propx-propy (default degree-BiBC) distribution for all nodes
    #####################################################################

    # Find the data type of each of the nodes in the network
    nw_df = pd.read_csv(network_file)
    nw_w_types = assign_type(nw_df)

    # Merge the data types dict with the node properties df
    node_props = node_props.reset_index()
    node_props['Data_type'] = node_props['index'].map(nw_w_types)
    node_props.dropna(inplace=True)
    unique_data_types = list(node_props.Data_type.unique())


    # list to dump inputs for plot_abundances via pickle
    pick_list = []
    pick_list.append(node_props)

    sns.set(font_scale=1.4, style = 'white')

    # Create plot for all data on one figure
    db_figure = plot_scatter(propx, propy, node_props, node_props.Data_type, None)
    db_figure.figure.savefig(outdir + propx + "_v_" + propy + "_distribution.png", bbox_inches='tight')
    print("Plot saved: " + outdir + propx + "_v_" + propy + "_distribution.png")

    degbibc_fig = db_figure.get_figure()
    degbibc_fig.clf()




    #####################################################################
    # Plot propx-propy (default degree-BiBC) distribution per data type
    #####################################################################

    # counter to hold df number, since the resulting figures need to be colored the same way
    # as the figure with all data types. Just gets used to specify index of sns.color_palette
    color_counter = 0

    # Create plot for each data type
    for i in unique_data_types:

        node_props_one_type = node_props[node_props['Data_type'].isin([i])]

        data_type_figure = plot_scatter(propx, propy, node_props_one_type, None, sns.color_palette()[color_counter])
        data_type_figure.figure.savefig(outdir + propx + "_v_" + propy + "_distribution_" + str(i) + "_nodes_only.png", bbox_inches='tight')
        print("Plot saved: " + outdir + propx + "_v_" + propy + "_distribution_" + str(i) + "_nodes_only.png has been saved")

        dt_fig = data_type_figure.get_figure()
        dt_fig.clf()
        color_counter += 1
        

    #############################################
    # Plot top propx-propy, all nodes
    #############################################

    # Sort values in df, putting the largest prop1 values first, then the largest prop2
    sorted_node_props = node_props.sort_values([propx, propy], ascending=[False, False])
    n_all_rows = sorted_node_props.shape[0]

    # Calculate what the top X% of nodes is if the user wants to plot top nodes based on percent
    if top_pct:
        n_top_rows = math.ceil(n_all_rows * (top_pct_num/100)) # Calculate how many parameters to plot
    else:
        n_top_rows = top_num

    # If there are not at least two top nodes (small network), then then plot the top two nodes
    if n_top_rows > 100:
        print("warning: The number of top nodes is greater than 100")
    elif n_top_rows < 2:
        n_top_rows = 2
    else:
        n_top_rows = n_top_rows


    sorted_node_props_top = sorted_node_props.head(n_top_rows)

    sns.set(font_scale=1.4, style = 'white')

    degbibc = sns.scatterplot(data=sorted_node_props_top, x=propx, y=propy, hue = "Data_type", s=100)

    degbibc.set(xlabel = propx, ylabel = propy)
    degbibc.figure.set_size_inches(8, 6)

    # Add labels to plot
    add_text(sorted_node_props_top, propx, propy, degbibc)

    if top_pct:
        degbibc.figure.savefig(outdir + propx + "_v_" + propy + "_distribution_top_" + str(top_pct_num) + "_percent.png", bbox_inches='tight')
        print("Plot saved: " + outdir + propx + "_v_" + propy + "_distribution_top_" + str(top_pct_num) + "_percent.png")
    else:
        degbibc.figure.savefig(outdir + propx + "_v_" + propy + "_distribution_top_" + str(top_num) + "_nodes.png", bbox_inches='tight')
        print("Plot saved: " + outdir + propx + "_v_" + propy + "_distribution_top_" + str(top_num) + "_nodes.png")

    degbibc_fig_small = degbibc.get_figure()
    degbibc_fig_small.clf()




    ##################################################
    # Plot top propx-propy nodes per data type
    ##################################################
    uniq_top = list(sorted_node_props['Data_type'].unique())

    color_counter = 0

    all_top_nodes = [] # list of all top nodes regardless of type for plotting abundances
    all_top_nodes_dict = {} # will hold key based on type and list of nodes for if they are the top nodes of that type for value

    for i in uniq_top:
        sorted_node_props_sub = sorted_node_props[sorted_node_props['Data_type']  == i]

        top_x_sorted_node_props_sub = sorted_node_props_sub.head(top_num_per_type)

        # for name in top_x_sorted_node_props_sub.iterrows():
        all_top_nodes.append(list(top_x_sorted_node_props_sub['index']))
        all_top_nodes_dict[i] = list(top_x_sorted_node_props_sub['index'])

        nnodes_to_plot = len(top_x_sorted_node_props_sub)

        if nnodes_to_plot < 3:
            print("\nToo few nodes of type " + str(i) + " to plot in " + propx + " v " + propy + " plot. Plot will not be made. Please enter a number 3 or greater.\n")
        else:
            scplot = plot_scatter(propx, propy, top_x_sorted_node_props_sub, None, sns.color_palette()[color_counter])
            scplot = add_text(top_x_sorted_node_props_sub, propx, propy, scplot)

            scplot.figure.savefig(outdir + propx + "_v_" + propy + "_distribution_top_" + str(nnodes_to_plot) + "_nodes_" + str(i) + "_only.png", bbox_inches='tight')


            print("Plot saved: " + outdir + propx + "_v_" + propy + "_distribution_top_" + str(nnodes_to_plot) + "_nodes_" + str(i) + "_only.png")
            dt_fig_top = scplot.get_figure()
            dt_fig_top.clf()

            color_counter += 1

    # Export the top nodes plotted for each node type to be used in plot_abundance.py
    all_top_nodes_flattened = [val for sublist in all_top_nodes for val in sublist]
    pick_list.append(all_top_nodes_flattened) # for plotting abundances
    pick_list.append(all_top_nodes_dict) # for plotting density plots per type
    pick_list.append(outdir)

    pickle.dump(pick_list, open(outdir + "inputs_for_downstream_plots.pickle", "wb")) # dump node_props to pickle for use in abundance code

