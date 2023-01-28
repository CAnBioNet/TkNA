Reconstruction and analysis of a network
Note that the following must be ran from the network_analysis directory that was created when you cloned the repo to your analysis directory above. 
1. Import the data and metadata for the run
	Usage
		python3 intake_data.py --data-dir <data directory> --out-file <output file>
	
	Example command
		python3 intake_data.py --data-dir ./input/all_genes/ --out-file ./output/all_data_and_metadata.cdf

	Inputs
		--data-dir: Path to the directory containing all experimental file(s), metadata file(s), and config file(s) 
		--out-file: path to file (with .cdf extension) that will be created 
	Outputs
		- A single .cdf file containing most information required for the next step


2. Run the correlations and comparisons, supplying the config file, which tells code which thresholds to apply. For larger datasets, use more cores.

	Usage
		python3 run.py --data-source <file_name> --config-file <config file> --out-file <zip directory>		
	
	Example command
	python3 run.py --data-source ./output/all_genes/all_data_and_metadata.cdf --config-file ./input/all_genes/cervical_cancer_config.json --out-file ./output/all_genes/network_output.zip
	
	Inputs
		--data-source: Path to the .cdf file created using intake_data.py
		--config-file: Path to the config file used for intake_data.py
		--out-file: path to zipped directory that will be created
	Outputs
		- A single zipped directory containing the analysis performed



3. Then convert the output files to csv 
	Usage
		python3 to_csv.py --data-file <zip file> --config-file <config file> --out-dir <output directory>
	
	Example command
		python3 to_csv.py --data-file ./output/all_genes/network_output.zip --config-file ./input/all_genes/cervical_cancer_config.json --out-dir ./output/all_genes/network_output
	
	Inputs
		--data-file: .zip file created with run.py
		--config-file: Path to the config file used for intake_data.py
		--out-dir: Path to the directory to output results to

	Outputs
		- all_comparisons.csv: all comparisons made, no statistical thresholds applied to file
		- correlations_bw_signif_measurables.csv: all correlations between parameters that passed differential change threshold. Correlations in this file are not filtered for statistical or causality criteria, but it contains all p-values, whether each edge is unexpected, and whether each edge makes it into the final network after applying statistical and causality criteria
		- network_output_comp.csv: whole network with nodes/edges under the user-defined statistical thresholds, and edges consistent in direction retained (unless otherwise specified in the configuration file). Unexpected edges are also removed from this file.
		- node_comparisons.csv: comparisons performed and found to be statistically significant, values listed in the name are the values of the thresholds applied in the config file
		- config_values.txt: All the user-specified options for making the network

4. Assess the quality of the reconstructed network
	Usage
		python assess_network.py --file <network file> 
	
	Example command
		python assess_network.py --file ./output/all_genes/network_output/ correlations_bw_signif_measurables.csv 
	
	Inputs
		--file: correlations_bw_signif_measurables.csv file created with to_csv.py
	
	Outputs
		- network_quality_assessment.txt: Contains the calculations (also sent to standard output) on the quality of the reconstructed network. Outputs to the same directory the input file is stored in.
		- network.pickle: A pickled file containing the network. Used as input to future steps. Outputs to the same directory the input file is stored in.

5. Calculate network and node properties of the reconstructed network
	Usage
		python calc_network_properties.py --pickle <file.pickle> --bibc --bibc-groups <choice> --bibc-calc-type <choice> --node-map <file.csv> --node-groups <group 1> <group 2>
	
	Example command
		python calc-network-properties.py --pickle network.pickle --bibc --bibc-groups node_types --bibc-calc-type rbc --node-map map_file.csv --node-groups micro pheno
	
	Inputs and arguments
		--pickle: network.pickle file created with assess_network.py
		--bibc: Flag for whether to compute Bipartite Betweenness Centrality (BiBC). This is highly recommended and also required for future steps
		--bibc-groups: Choice for what to compute BiBC on, either distinct groups (node_types) or on the two most modular regions of the network (found using the Louvain method)
		--bibc-calc-type: Choice for whether to normalize based on amount of nodes in each group (rbc) or not (bibc)?
		--node-map: csv file containing the name of nodes in the first column and the type of the node (gene, phenotype, microbe, etc.) in the second column
		--node-groups: Required if node_types is specified for --bibc-groups. It’s the two groups of nodes to calculate BiBC/RBC on. The types must be present in the --node-map file
	
	Outputs
		- network_properties.txt: Tab-delimited .txt file of calculated network properties
		- subnetwork_properties.txt: Tab-delimited .txt file of calculated subnetwork properties
		- node_properties.txt: Tab-delimited .txt file of calculated node properties

6. Create random networks
	Usage
		python create_random_networks.py --template-network <file.pickle> --networks-file <file.pickle> 
	
	Example command
		python create_random_networks.py --template-network network.pickle --networks-file all_random_nws.pickle 
	
	Inputs
		--template-network: The pickled network file output by assess_network.py
		--networks-file: File to output pickled network list to
	
	Outputs
		- A single pickle file containing all created networks

7. Analyze random networks 
	Usage
		python compute_network_stats.py --networks-file <file.pickle> --bibc-groups <choice> --bibc-calc-type <choice> --stats-file <file.pickle> --node-map <file.csv> --node-groups <group1> <group2>
	
	Example command
		python compute_network_stats.py --networks-file all_random_nws.pickle --bibc-groups node_types --bibc-calc-type rbc --stats-file all_rand_network_results.pickle --node-map map.csv --node-groups gene pheno
	
	Inputs
		--networks-file: pickled file created with create_random_networks.py that contains all random networks previously created
		--bibc-groups: Group nodes for BiBC based on type or modularity
		--bibc-calc-type: Compute raw BiBC or normalize (rbc)
		--stats-file: Pickle file to output the network stats to
		--node-map: CSV file mapping nodes to their types. Required if node_types is specified for --bibc-groups.
		--node-groups: Two types of nodes to use for BiBC grouping. Required if node_types is specified for --bibc-groups.
		
	Outputs
		- A single .pickle file with degree/BiBC results of all random networks

8. Condense random network outputs into one file
	Usage
		python synthesize_network_stats.py --network-stats-file <file.pickle> --synthesized-stats-file <file.csv> 
	
	Example command
		python synthesize_network_stats.py --network-stats-file all_rand_network_results.pickle --synthesized-stats-file top_node_each_random_nw.csv
	
	Inputs
		--network-stats-file: pickled file created with compute_network_stats.py
		--synthesized-stats-file: Name of the .csv file that will be created 
	
	Outputs
		- A single .csv file that contains the top node, sorted first by BiBC and then by Node_degrees (unless otherwise specified with --flip-priority), for each of the random networks

9. Create dot plots for node properties
	Usage
		python dot_plots.py --pickle <file.pickle> --node-props  <file.txt> --network-file <file.csv> --propx BiBC --propy Node_degrees --top-num <integer> --top-num-per-type <integer>
	
	Example command
	python dot_plots.py --pickle network.pickle --node-props node-properties.txt --network-file network_output_comp.csv --propx BiBC --propy Node_degrees --top-num 5 --top-num-per-type 3
	
	Inputs
		--pickle: pickled file created with assess_network.py
		--node-props  : node_properties.txt file created with calc_network_properties.py
		--network-file: network_output_comp.csv file created with to_csv.py 
		--propx: Node property to plot on X-axis
		--propy: Node property to plot on Y-axis.
		--top-num: Number of nodes you want to zoom in to on the property v property plot
		--top-num-per-type: The number of nodes to plot for each data type when zoomed in on the plot

	Default outputs
		- degree_distribution_dotplot.png: Distribution of the number of nodes which each degree in the network
		- <propx>_v_<propy>_distribution.png: A dot plot of user-specified node properties
		- <propx>_v_<propy>_distribution_<node_type>_nodes_only.png: Same as previous plot, but with just the nodes from each data type. There will be one plot produced for each data type
		- <propx>_v_<propy>_distribution_top_<top-num>_nodes.png: Same as the second plot, but zoomed in on the top nodes 
		- <propx>_v_<propy>_distribution_top_<top-num-per-type>_nodes_<data_type>_only.png: same as third plot, but zoomed in on the top nodes per data type. 

10. Plot abundances of nodes
	Usage
		python plot_abundance.py --pickle <file.pickle> --abund-data <list of files> --metadata <list of files> --color-group <choice> --x-axis <choice> 
	
	Example command
		python plot_abundance.py --pickle inputs_for_downstream_plots.pickle --abund-data Expt1.csv Expt2.csv --metadata Expt1_meta.csv Expt2_meta.csv --color-group Treatment --x-axis Experiment 
	
	Inputs
		--pickle: inputs_for_downstream_plots.pickle file output by dot_plots.py
		--abund_data: List of data files containing expressions/abundances
		--metadata: List of metadata files containing Experiment/Treatment columns
		--color-group: Variable to color the plot by
		--x-axis: Variable you wish to group the data by on the x-axis

	Outputs
		- One boxplot for each of the top nodes (found in dot_plots.py) as well as additional plots if specified with the optional argument --nodes-to-plot

11. Plot density/contour plot of distributions of top nodes in a random network and overlay the nodes from the actual reconstructed network
	Usage
		python plot_density.py --rand-net <file.csv> --pickle <file.pickle>
	
	Example command
		python plot_density.py --rand-net synthesized_networks.csv --pickle inputs_for_downstream_plots.pickle
	
	Inputs
		--rand-net: file output by synthesize_network_stats.py
		--pickle: inputs_for_downstream_plots.pickle file output by dot_plots.py

	Default outputs
		- density_plot_with_top_nodes_from_dotplots.png: contour plot with the top nodes (found in dot_plots.py) from the real reconstructed network overlaid on top
		- density_plot_with_top_<data_type>_nodes_only.png: Same as previous, but contains just one data type per output file

Additional (optional) analyses to perform
To identify clusters of nodes in the network, you can use either Infomap or the Louvain method.

Infomap
	Usage
		python infomap_assignment.py --pickle <file.pickle>
		
	Example command
		python infomap_assignment.py --pickle network.pickle
		
	Inputs
		--pickle: network.pickle file output by assess_network.py

	Output
		- network_infomap_partition.csv: .csv file containing the name of the node in column 1 and the subnetwork number it was assigned in column 2. 

Louvain
	Usage
		python infomap_assignment.py --pickle <file.pickle>
	
	Example command
		python Louvain_partition.py --pickle network.pickle
	
	Inputs
		--pickle: network.pickle file output by assess_network.py

	Output
		- network_infomap_partition.csv: .csv file containing the name of the node in column 1 and the subnetwork it was assigned in column 2. 

Find shortest paths between subnetworks
	Usage
		python find_all_shortest_paths_bw_subnets.py --network <file.pickle> --node-map <map.csv> --node-groups <group1> <group2>
	
	Example command
		python find_all_shortest_paths_bw_subnets.py --network <network.pickle> --node-map map.csv --node-groups gene pheno
	
	Inputs
		--network: network.pickle file output by assess_network.py
		--node-map: Mapping file (csv) of nodes and their subnetworks
		--node-groups: The two groups in the mapping file you want to find the shortest paths between

	Output
		- shortest_path_bw_<group1>_and_<group2>_results.csv: .csv file containing the name of each node in each pair in columns 1 and 2, as well as the shortest path length between that pair in column 3 and the number of shortest paths for the pair in column 4

