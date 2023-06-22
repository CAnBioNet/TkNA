# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 09:59:18 2023

@author: Nolan K Newman <newmanno@oregonstate.edu>
"""

import argparse
import pandas as pd
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import Formatter
import matplotlib as mpl
from functools import partial
import numpy as np
from adjustText import adjust_text
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import sys

class LowHighFormatter(Formatter):
    def format_data(self, value):
        return "Low" if value == 0 else "High"
    def format_ticks(self, values):
        return [self.format_data(value) for value in values]

if __name__ == '__main__':


    def plot_kde(xprop, yprop, col):
        '''
        Function that creates a seaborn kde plot, given two property names and a color palette

        Arguments:
            - xprop: column in data frame of user-defined property to plot on the X-axis
            - yprop: column in data frame of user-defined property to plot on the Y-axis
            - col: name of color palette to use
        '''
        sns.set(font_scale=1.5)
        plt.figure(figsize=(10,8))

        ax = plt.axes()
        #cbaxes = inset_axes(ax, width="100%", height="10%", loc=8)
        cax, kw = mpl.colorbar.make_axes(ax, location="top", orientation="horizontal", pad=0.075)

        # The following works but cbar labels are on top
        fig1 = sns.kdeplot(xprop, yprop, shade=True, ax=ax, cmap=col, cbar=True, cbar_ax=cax, cbar_kws=dict(ticks=[0,9], format=LowHighFormatter(), label='Likelihood of random finding', use_gridspec=False, **kw))

        cax.xaxis.set_ticks_position('bottom')


        #fig1 = sns.kdeplot(xprop, yprop, shade=True, cmap=col)
        #cb = mpl.colorbar.ColorbarBase(fig1, ticks=[0,9], format=LowHighFormatter(), label='Likelihood of random finding', orientation='horizontal')
        #fig1.xaxis.set_ticks_position('top')

        #fig1 = sns.kdeplot(xprop,yprop,shade=True,cmap=col, cbar=True, cbar_kws=dict(ticks = [0,9], label='Likelihood of random finding', use_gridspec = False, location = "top"))
        #fig2 = sns.kdeplot(xprop,yprop,shade=True)


        #cbar = fig1.collections[0].colorbar
        #print(fig1)

        #cbar.set_ticklabels(['Low', 'High'])

        #plt.draw()


        return(fig1)

    def add_text(dic_obs, dic_prob, ax):
        '''
        Function that adds text to kde plot from a dictionary

        Arguments:
            - dic_obs: dictionary keyed on node name that contains a list of the observed Degree and BiBC for each node as values
            - dic_prob: dictionary keyed on node name that contains the likelihood of finding that node randomly as values
        '''
        dots = []
        texts_bibcdeg = []
        bboxes = []

        for k,v in dic_obs.items():
            dot = ax.plot(v[1],v[0], "ko", zorder = 10)[0]
            dots.extend(ax.plot(v[1],v[0], "ko", zorder = 10))

            texts_bibcdeg.append(ax.text(v[1], v[0], k + "\n" + dic_prob[k], ha='center', va='center', fontsize=12))
            bbox = dot.get_window_extent(ax.get_figure().canvas.get_renderer())
            #bb_height_mult = 10000 * bbox.height
            #bb_width_mult = 10000 * bbox.width
            #bbox = mpl.transforms.Bbox([[bbox.xmin-bb_width_mult, bbox.ymin-bb_height_mult], [bbox.xmax+bb_width_mult, bbox.ymax+bb_height_mult]])
            bbox = bbox.expanded(1.1, 1.1)
            bboxes.append(bbox)

        #print(bboxes)
        bkgnd_cols = mpl.cm.get_cmap('Reds')
        bkgnd_col = bkgnd_cols(13)
        print(bkgnd_col)

        ax.set_facecolor(bkgnd_col)
        ax.grid(False)
        adjust_text(texts_bibcdeg, ax=ax, objects = bboxes, time_lim=3)


        return(texts_bibcdeg)

    def check_rand_net_output(df):
        '''
        Function that determines whether a given data frame is suitable for density plotting based on its range of degrees

        There is a limitation to the code, where matplotlib will throw an error if the range in degrees is not at least 2,
        so this function checks if that is the case

        Arguments:
            - df: the random network data frame, pre-filtered to remove any N/As that may have been present
        '''
        max_degree = max(df['Degree'])
        min_degree = min(df['Degree'])
        range_degree = max_degree - min_degree

        if range_degree < 2:
            return(False)
        else:
            return(True)




    testing = False

    if testing:
        rand_net = "random_networks_synthesized_Nolan.csv"
        pick = "inputs_for_downstream_plots_Nolan.pickle"
        # nnodes_density = 4
        # nnodes_density_per_type = 3
        # density_node_types = ['gene']
        bibc_name = "BiBC"
        nodes_to_plot = ['OTU_259506', 'fasting leptin']

    else:
        print("")
        parser = argparse.ArgumentParser(description="Example command: python plot_density.py random_network_condensed.csv inputs_for_downstream_plots.pickle", add_help=False)

        requiredArgGroup  = parser.add_argument_group('Required arguments')
        requiredArgGroup.add_argument("--rand-net", type=str, dest="rand_net", help="Random network file output by synthesize_network_stats.py", required=True)
        requiredArgGroup.add_argument("--pickle", type=str, help="inputs_for_downstream_plots.pickle file output by dot_plots.py", required=True)
        requiredArgGroup.add_argument("--bibc-name", type=str, dest="bibc_name", help="The 'name' of the BiBC calculation performed, from the node_properties.txt file. Example: BiBC_microbe_pheno (if BiBC was calculated between the microbe and pheno groups)", required=True)

        optionalArgGroup  = parser.add_argument_group('Optional arguments')
        optionalArgGroup.add_argument("-h", "--help", action="help", help="Show this help message and exit")
        optionalArgGroup.add_argument("--nodes-to-plot", type=str, nargs='+', dest="nodes_to_plot", help="A list of nodes to plot on the density plot")


        # parser = argparse.ArgumentParser(description="Example command: python plot_density.py random_network_condensed.csv inputs_for_downstream_plots.pickle")
        # parser.add_argument("rand_net", type=str, help="Random network file output by condense_random_networks.py")
        # parser.add_argument("pickle", type=str, help="inputs_for_downstream_plots.pickle file output by dot_plots.py")
        # # parser.add_argument("--nnodes_density", type=str, help="The number of total nodes you want to plot on top of the density plot, regardless of node type")
        # # parser.add_argument("--nnodes_density_per_type", type=str, help="The number of nodes you want to plot on top of the density plot of each node type")
        # # parser.add_argument("--nodes_to_plot", type=list, help="The node types to plot on top of density plot")
        # parser.add_argument("--dens_list_to_plot", type=str, nargs='+', help="A list of nodes to plot on the density plot")

        args = parser.parse_args()

        rand_net = args.rand_net
        pick = args.pickle
        # nnodes_density = args.nnodes_density
        # nnodes_density_per_type = args.nnodes_density_per_type
        # node_props = args.density_node_types
        bibc_name = args.bibc_name
        nodes_to_plot = args.nodes_to_plot


    ########################
    # Import data
    ########################
    rand_nets = pd.read_csv(rand_net, header = 0)
    rand_nets = rand_nets.iloc[:,2:]

    # drop any N/A values from the dataframe, which would be present if its a small
    # and dense network where all nodes of one BiBC group connect to all nodes of the other group
    rand_nets = rand_nets.dropna(axis=0)
    rand_nets_check = check_rand_net_output(rand_nets)

    print(rand_nets)

    # If the random network data frame passes the check:
    if rand_nets_check:

        ##############################################################
        # Plot random network density for top nodes from the dotplots
        ##############################################################
        # Load the node_props pickle file that already contains the data type
        p = open(pick, "rb")
        pload = pickle.load(p)
        node_props = pload[0]
        top_nodes = pload[1]
        top_nodes_per_type = pload[2]
        # plotdir = "./" # location of previously created plots directory

        # !!!!!!!!!!Must uncomment the following before publishing
        plotdir = pload[3] # location of previously created plots directory

        print(node_props)

        # First find how many top nodes have both a degree and BiBC greater than the indicated node
        all_nodes_dict = {}

        BiBC_deg_sorted = node_props.sort_values([bibc_name, 'Node_degrees'], ascending=[False, False])
        sorted_node_props_cut = BiBC_deg_sorted[['index','Node_degrees', bibc_name]]
        sorted_node_props_cut_top = sorted_node_props_cut[sorted_node_props_cut['index'].isin(top_nodes)]

        # populate dictionary so you get node degrees and BiBC as values for each key (node)
        for index, row in sorted_node_props_cut.iterrows():
            all_nodes_dict[row[0]] = [row[1], row[2]]

        prob_dict = {} # Dictionary for plotting probabilities of nodes

        sorted_node_props_top_list = sorted_node_props_cut_top['index'].tolist()

        # Calculate probabilities for all nodes, where chance to find a node randomly is the "square" that can be drawn on the density plot,
        # meaning it counts the number of nodes found randomly in networks that have both a greater degree and BiBC than each calculated node
        with open(plotdir + "probabilities_to_randomly_find_nodes.csv", "w") as file:

            file.write("Node,Observed degree, Observed BiBC, Probability\n")

            for node in BiBC_deg_sorted['index'].tolist():
                greater_than_density_pval_node = rand_nets[(rand_nets['Degree'] > all_nodes_dict[node][0]) & (rand_nets['BiBC'] > all_nodes_dict[node][1])]
                rand_net_pval = round(len(greater_than_density_pval_node)/len(rand_nets) * 100, 2)
                prob_dict[node] = str(rand_net_pval) + "%"
                file.write(node + "," + str(all_nodes_dict[node][0]) + "," + str(all_nodes_dict[node][1]) + "," + str(rand_net_pval) + "%\n")

        top_nodes_dict = {}

        # list_nodes = sorted_node_props_cut_top['index'].tolist()
        sorted_node_props_top_cut = sorted_node_props_cut_top[['index','Node_degrees', bibc_name]]

        # Find top nodes to plot (same nodes as in zoomed in Degree-BiBC plot)
        for index, row in sorted_node_props_top_cut.iterrows():
            top_nodes_dict[row[0]] = [row[1], row[2]]



        kdefig = plot_kde(rand_nets['BiBC'], rand_nets['Degree'], 'Reds')
        texts = add_text(top_nodes_dict, prob_dict, kdefig)


        #kdefig = kdefig.get_figure()
        plt.savefig(plotdir + "density_plot_with_top_nodes_from_dotplots.png", bbox_inches='tight', dpi=600)
        print("Plot saved: " + plotdir + "density_plot_with_top_nodes_from_dotplots.png")
        plt.clf()




        ###############################################################
        # Plot random network density for top nodes per data type
        ###############################################################

        # Only plot the top X nodes, where X is a user-defined number
        for k,v in top_nodes_per_type.items():
            # print(k,v)

            top_dens_nodes_per_type = BiBC_deg_sorted[BiBC_deg_sorted['index'].isin(v)]
            top_dens_nodes_per_type = top_dens_nodes_per_type[['index','Node_degrees', bibc_name]]

            # Convert df to dictionary to plot
            top_dens_nodes_per_type_dict = {}
            for index, row in top_dens_nodes_per_type.iterrows():
                top_dens_nodes_per_type_dict[row[0]] = [row[1], row[2]]

            kdefig = plot_kde(rand_nets['BiBC'], rand_nets['Degree'], 'Reds')
            kdefig = add_text(top_dens_nodes_per_type_dict, prob_dict, kdefig)

            plt.savefig(plotdir + "density_plot_with_top_" + k + "_nodes_only.png", bbox_inches='tight', dpi=600)
            print("Plot saved: " + plotdir + "density_plot_with_top_" + k + "_nodes_only.png")
            plt.clf()




        ###############################################################
        # Plot random network density for user-defined list of nodes
        ###############################################################

        if nodes_to_plot:
            # Check that the nodes the user wants are in the df
            if set(nodes_to_plot).issubset(BiBC_deg_sorted['index'].tolist()) == False:
                raise Exception("ERROR: Please make sure the desired nodes are in the node properties file.")

            user_selected_nodes_dict = {} # dictionary to hold user-specified genes, degree, and bibc

            # Only plot the user-specified nodes from dens_list_to_plot
            for i in nodes_to_plot:
                # print(i)
                selected_node = BiBC_deg_sorted[BiBC_deg_sorted['index']  == i]
                selected_node = selected_node[['index','Node_degrees', bibc_name]]

                for index, row in selected_node.iterrows():
                    user_selected_nodes_dict[row[0]] = [row[1], row[2]]

            kdefig = plot_kde(rand_nets['BiBC'], rand_nets['Degree'], 'Reds')
            kdefig = add_text(user_selected_nodes_dict, prob_dict, kdefig)

            plt.savefig(plotdir + "density_plot_with_user_defined_node_list.png", bbox_inches='tight', dpi=600)
            print("Plot saved: " + plotdir + "density_plot_with_user_defined_node_list.png")
            plt.clf()

    # If network is too small/too dense, don't create figure
    else:
        raise Exception("\nError: random network density cannot be plotted. This is likely because you have too small and too dense of a network. This error will occur if the range of degrees in " + rand_net + " is less than 2, meaning all nodes have too similar of a degree.\n")

