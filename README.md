# TkNA

## Usage

**NB:** If you have not already set up a conda environment and cloned the repo, please complete the steps in `Creating a conda environment.txt` before running the following commands.

Note that script paths are given relative to the top-level repository directory.

### 1. Normalize data
Data should be normalized prior to running TkNA. Examples of normalization methods can be found in the TkNA manuscript.

### 2. Format data and set statistical thresholds
Data must be formatted in the format specified in the TkNA manuscript.

### 3. Import the data and metadata for the run

#### Usage
```
python ./reconstruction/intake_data.py --data-dir <data directory> --out-file <output file>
```

#### Example command
```
python reconstruction/intake_data.py --data-dir ./project_folder/input/ --out-file ./project_folder/output/all_data_and_metadata.zip
```

#### Inputs
 - `--data-dir`: Path to the directory containing all experimental file(s), metadata file(s), and config file(s)
 - `--out-file`: path to file (with `.zip` extension) that will be created

#### Outputs
A single `.zip` file containing most information required for the next step

### 4. Run the correlations and comparisons, supplying the config file, which tells code which thresholds to apply. For larger datasets, use more cores.

#### Usage
```
python ./reconstruction/run.py --data-source <file_name> --config-file <config file> --out-file <zip file>
```

#### Example command
```
python ./reconstruction/run.py --data-source ./project_folder/output/all_data_and_metadata.zip --config-file ./project_folder/input/config.json --out-file ./project_folder/output/network_output.zip
```

#### Inputs
 - `--data-source`: Path to the `.zip` file created using `intake_data.py`
 - `--config-file`: Path to the config file used for `intake_data.py`
 - `--out-file`: Path to zipped directory that will be created

#### Outputs
 - A single zipped directory containing the analysis performed

### 5. Convert the output files to csv

#### Usage
```
python ./reconstruction/to_csv.py --data-file <zip file> --config-file <config file> --out-dir <output directory>
```

#### Example command
```
python ./reconstruction/to_csv.py --data-file ./project_folder/output/network_output.zip --config-file ./project_folder/input/config.json --out-dir ./project_folder/output/network_output
```

#### Inputs
 - `--data-file`: `.zip` file created with `run.py`
 - `--config-file`: Path to the config file used for `intake_data.py`
 - `--out-dir`: Path to the directory to output results to

#### Outputs
 - `all_comparisons.csv`: all comparisons made, no statistical thresholds applied to file
 - `correlations_bw_signif_measurables.csv`: all correlations between parameters that passed differential change threshold. Correlations in this file are not filtered for statistical or causality criteria, but it contains all p-values, whether each edge is unexpected, and whether each edge makes it into the final network after applying statistical and causality criteria
 - `network_output_comp.csv`: whole network with nodes/edges under the user-defined statistical thresholds, and edges consistent in direction retained (unless otherwise specified in the configuration file). Unexpected edges are also removed from this file.
 - `node_comparisons.csv`: comparisons performed and found to be statistically significant, values listed in the name are the values of the thresholds applied in the config file
 - `config_values.txt`: All the user-specified options for making the network

### 6. Assess the quality of the reconstructed network

#### Usage
```
python ./analysis/assess_network.py --file <network file> --out-dir <directory>
```

#### Example command
```
python ./analysis/assess_network.py --file ./project_folder/output/network_output/correlations_bw_signif_measurables.csv --out-dir ./project_folder/output/network_output/
```

#### Inputs
 - `--file`: `correlations_bw_signif_measurables.csv` file created with `to_csv.py`
 - `--out-dir`: Path to the directory to output results to


#### Outputs
 - `network_quality_assessment.csv`: Contains the network quality statistics of the reconstructed network, calculated per edge type.
 
### 7. Identify clusters of nodes (OPTIONAL)...

### ...using Infomap

#### Usage
```
python ./analysis/infomap_assignment.py --network <file> --network-format <format> --map <file.csv> --out-dir <directory>
```

#### Example command
```
python ./analysis/infomap_assignment.py --network ./project_folder/output/network_output/network_output_comp.csv --network-format csv --map ./project_folder/input/type_map.csv --out-dir ./project_folder/output/network_output/
```

#### Inputs
 - `--network`: `The path to the network file, either in .pickle or .csv format
 - `--network-format`: Format of the network file; Either use 'pickle' or 'csv' (must be in .csv format and have 'partner1' and 'partner2' as the headers for the two node columns, e.g. the network_output_comp.csv from to_csv.py)
 - `--map`: CSV file with the name of the node in the first column and its data type in the second column
 - `--out-dir`: Path to the directory to output results to


#### Output
 - `network_infomap_partition.csv`: CSV file containing the name of the node in column 1 and the subnetwork number it was assigned in column 2. 

### ...using the Louvain method

#### Usage
```
python ./analysis/louvain_partition.py --network <file> --network-format <format> --map <file.csv> --out-dir <directory>	
```

#### Example command
```
python ./analysis/louvain_partition.py --network ./project_folder/output/network_output/network_output_comp.csv --network-format csv --map ./project_folder/input/type_map.csv --out-dir ./project_folder/output/network_output/
```

#### Inputs and arguments
 - `--network`: `The path to the network file, either in .pickle or .csv format
 - `--network-format`: Format of the network file; Either use 'pickle' or 'csv' (must be in .csv format and have 'partner1' and 'partner2' as the headers for the two node columns, e.g. the network_output_comp.csv from to_csv.py)
 - `--map`: CSV file with the name of the node in the first column and its data type in the second column
 - `--out-dir`: Path to the directory to output results to

#### Output
 - `network_louvain_partition.csv`: CSV file containing the name of the node in column 1 and the subnetwork number it was assigned in column 2. 

### 8. Perform functional enrichment analysis for groups of nodes (OPTIONAL)

Functional enrichment of the resulting partitions can be performed using the methods mentioned in the TkNA manuscript.

### 9. Find the distance (shortest path) between two pathways (OPTIONAL)

Pathways closer to one another potentially interact more than those that are further away.

#### Usage
```
python ./analysis/find_all_shortest_paths_bw_subnets.py --network <file.pickle> --network-format <format> --map <map.csv> --node-groups <group1> <group2> --out-dir <directory> 
```

#### Example command
```
python ./analysis/find_all_shortest_paths_bw_subnets.py --network ./project_folder/output/network_output/network.pickle --node-map ./project_folder/input/map_file.csv --node-groups gene pheno --out-dir ./project_folder/output/network_output/
```

#### Inputs and arguments
 - `--network`: The path to the network file, either in .pickle or .csv format
 - `--network-format`: Format of the network file; Either use 'pickle' with the network.pickle file output made by assess_network.py (if network was reconstructed using the TkNA pipeline) or 'csv' if the network was reconstructed using an alternative pipeline (must be in .csv format and have 'partner1' and 'partner2' as the headers for the two node columns
 - `--map`: CSV file with the name of the node in the first column and its data type in the second column
 - `--out-dir`: Path to the directory to output results to

#### Output
 - `shortest_path_bw_<group1>_and_<group2>_results.csv`: CSV file containing the name of each node in each pair in columns 1 and 2, as well as the shortest path length between that pair in column 3 and the number of shortest paths for the pair in column 4.

### 10. Calculate network topology parameters (e.g. degree, BiBC, etc.)

#### Usage
```
python ./analysis/calc_network_properties.py --network <file.csv> --bibc --bibc-groups <choice> --bibc-calc-type <choice> --map <file.csv> --node-groups <group 1> <group 2> --out-dir <directory>
```

#### Example command:
```
python ./analysis/calc_network_properties.py --network ./project_folder/output/network_output/network_output_comp.csv --bibc --bibc-groups node_types --bibc-calc-type rbc --map ./project_folder/input/type_map.csv --node-groups microbe pheno --out-dir ./project_folder/output/network_output/
```

#### Inputs and arguments
 - `--network`: The network file in CSV format containing the reconstructed network. Must have columns called 'partner1' and 'partner2'.
 - `--bibc`: Flag for whether to compute Bipartite Betweenness Centrality (BiBC). This is highly recommended and also required for future steps
 - `--bibc-groups`: Choice for what to compute BiBC on, either distinct groups (node_types) or on the two most modular regions of the network (found using the Louvain method)
 - `--bibc-calc-type`: Choice for whether to normalize based on the number of nodes in each group (rbc) or not (bibc)
 - `--node-map`: CSV file containing the name of nodes in the first column and the type of the node (gene, phenotype, microbe, etc.) in the second column
 - `--node-groups`: Required if node_types is specified for --bibc-groups. It’s the two groups of nodes to calculate BiBC/RBC on. The types must be present in the --node-map file
 - `--out-dir`: Path to the directory to output results to

#### Output
 - `network_properties.txt`: Tab-delimited .txt file of calculated network properties
 - `subnetwork_properties.txt`: Tab-delimited .txt file of calculated subnetwork properties
 - `node_properties.txt`: Tab-delimited .txt file of calculated node properties

### 11. Create random networks

#### Usage
```
python ./random_networks/create_random_networks.py --template-network <file.pickle> --networks-file <file.zip> 	
```

#### Example command
```
python ./random_networks/create_random_networks.py --template-network ./project_folder/output/network_output/network.pickle --networks-file ./project_folder/output/network_output/all_random_nws.zip 
```

#### Inputs
 - `--template-network`: The pickled network file output by `calc_network_properties.py`
 - `--networks-file`: .zip file to output zipped networks to
 - `--num-networks`: optional; default 10,000; number of random networks to create
 
#### Output
 - A single .zip file containing all created networks
 
### 12. Analyze random networks

#### Usage
```
python ./random_networks/compute_network_stats.py --networks-file <file.zip> --bibc-groups <choice> --bibc-calc-type <choice> --stats-file <file.zip> --node-map <file.csv> --node-groups <group1> <group2>
```

#### Example command
```
python ./random_networks/compute_network_stats.py --networks-file ./project_folder/output/network_output/all_random_nws.zip --bibc-groups node_types --bibc-calc-type rbc --stats-file ./project_folder/output/network_output/random_network_analysis.zip --node-map ./project_folder/input/type_map.csv --node-groups microbe pheno
```

#### Inputs
 - `--networks-file`: zip file created with `create_random_networks.py` that contains all random networks previously created
 - `--bibc-groups`: Group nodes for BiBC based on type or modularity
 - `--bibc-calc-type`: Compute raw BiBC or normalize (`rbc`)
 - `--stats-file`: zip file to output the network stats to
 - `--node-map`: CSV file mapping nodes to their types. Required if `node_types` is specified for `--bibc-groups`.
 - `--node-groups`: Two types of nodes to use for BiBC grouping. Required if `node_types` is specified for `--bibc-groups`.

#### Outputs
 - A single zip file with degree/BiBC results of all random networks

### 13. Condense random network outputs into one file

#### Usage
```
python ./random_networks/synthesize_network_stats.py --network-stats-file <file.zip> --synthesized-stats-file <file.csv> 
```

#### Example command
```
python ./random_networks/synthesize_network_stats.py --network-stats-file ./project_folder/output/network_output/random_network_analysis.zip --synthesized-stats-file ./project_folder/output/network_output/random_networks_synthesized.csv
```

#### Inputs
 - `--network-stats-file`: zip file created with `compute_network_stats.py`
 - `--synthesized-stats-file`: Name of the `.csv` file that will be created

#### Outputs
 - A single `.csv` file that contains the top node, sorted first by BiBC and then by Node_degrees (unless otherwise specified with `--flip-priority`), for each of the random networks

### 14. Create dot plots for node properties

#### Usage
```
python ./visualization/dot_plots.py --pickle <file.pickle> --node-props  <file.txt> --network-file <file.csv> --propx BiBC --propy Node_degrees --top-num <integer> --top-num-per-type <inte-ger> --plot-dir <directory> --file-dir <directory>
```

#### Example command
```
python ./visualization/dot_plots.py --pickle ./project_folder/output/network_output/network.pickle --node-props ./project_folder/output/network_output/node_properties.txt --network-file ./project_folder/output/network_output/network_output_comp.csv --propx BiBC_microbe_pheno --propy Node_degrees --top-num 5 --top-num-per-type 3	--plot-dir ./project_folder/output/network_output/plots/ --file-dir ./project_folder/output/network_output/
```

#### Inputs
 - `--pickle`: pickled file created with `assess_network.py`
 - `--node-props`: `node_properties.txt` file created with `calc_network_properties.py`
 - `--network-file`: `network_output_comp.csv` file created with `to_csv.py`
 - `--propx`: Node property to plot on X-axis
 - `--propy`: Node property to plot on Y-axis.
 - `--top-num`: Number of nodes you want to zoom in to on the property v property plot
 - `--top-num-per-type`: The number of nodes to plot for each data type when zoomed in on the plot
 - `--plot-dir`: Path to the directory to output the resulting plots to
 - `--file-dir`: Path to the directory to save the resulting inputs_for_downstream_plots.pickle file to

#### Default outputs
 - `degree_distribution_dotplot.png`: Distribution of the number of nodes which each degree in the network
 - `<propx>_v_<propy>_distribution.png`: A dot plot of user-specified node properties
 - `<propx>_v_<propy>_distribution_<node_type>_nodes_only.png`: Same as previous plot, but with just the nodes from each data type. There will be one plot produced for each data type
 - `<propx>_v_<propy>_distribution_top_<top-num>_nodes.png`: Same as the second plot, but zoomed in on the top nodes
 - `<propx>_v_<propy>_distribution_top_<top-num-per-type>_nodes_<data_type>_only.png`: same as third plot, but zoomed in on the top nodes per data type.
 - `inputs_for_downstream_plots.pickle`: contains information for future commands
 
### 15. Create abundance plots

#### Usage
```
python ./visualization/plot_abundance.py --pickle <file.pickle> --abund-data <list of files> --metadata <list of files> --x-axis <choice> --group-names <list of names> --group-colors <list of color names> 
```

#### Example command
```
python ./visualization/plot_abundance.py --pickle ./project_folder/output/network_output/inputs_for_downstream_plots.pickle --abund-data ./project_folder/input/Experiment1.csv ./project_folder/input/Experiment2.csv --metadata ./project_folder/input/Experiment1_group_map.csv ./project_folder/input/Experiment2_group_map.csv --x-axis Experiment --group-names ileum8wkHFHS ileum8wkNCD stool4wkHFHS stool8wkHFHS --group-colors orange blue red green
```

#### Inputs
 - `--pickle`: `inputs_for_downstream_plots.pickle` file output by `dot_plots.py`
 - `--abund_data`: List of data files containing expressions/abundances
 - `--metadata`: List of metadata files containing Experiment/Treatment columns
 - `--x-axis`: Variable you wish to group the data by on the x-axis
 - `--group-names`: A list of the names of the treatment groups; must match the order of names in --group-colors
 - `--group-colors`: A list of the names of the colors to use for each specified group; must match the order of colors in --group-names. Accepted colors can be found at https://matplotlib.org/stable/gallery/color/named_colors.html

#### Outputs
 - One boxplot for each of the top nodes (found in `dot_plots.py`) as well as additional plots if specified with the optional argument `--nodes-to-plot`

### 16. Plot density/contour plot of distributions of top nodes in a random network and overlay the nodes from the actual reconstructed network

#### Usage
```
python ./visualization/plot_density.py --rand-net <file.csv> --pickle <file.pickle> --bibc-name <name>
```

#### Example command
```
python ./visualization/plot_density.py --rand-net ./project_folder/output/network_output/random_networks_synthesized.csv --pickle ./project_folder/output/network_output/inputs_for_downstream_plots.pickle --bibc-name BiBC_microbe_pheno
```

#### Inputs
 - `--rand-net`: file output by `synthesize_network_stats.py`
 - `--pickle`: `inputs_for_downstream_plots.pickle` file output by `dot_plots.py`
 - `--bibc-name`: The 'name' of the BiBC calculation performed, from the node_properties.txt file. Example: BiBC_microbe_pheno (if BiBC was calculated between the microbe and pheno groups)
 
#### Default outputs
 - `density_plot_with_top_nodes_from_dotplots.png`: contour plot with the top nodes (found in `dot_plots.py`) from the real reconstructed network overlaid on top
 - `density_plot_with_top_<data_type>_nodes_only.png`: Same as previous, but contains just one data type per output file
