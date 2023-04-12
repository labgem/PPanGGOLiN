#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import time
import os
from itertools import combinations, product
from collections.abc import Callable
from collections import defaultdict 
from multiprocessing.pool import Pool
from itertools import islice
import time
import sys
import numpy as np

# installed libraries
from tqdm import tqdm
import networkx as nx

import plotly.express as px
from scipy.stats import gaussian_kde
import pandas as pd
from sklearn.metrics.cluster import adjusted_mutual_info_score, mutual_info_score, adjusted_rand_score, rand_score, normalized_mutual_info_score



# local libraries
from ppanggolin.genome import Organism, Contig
from ppanggolin.pangenome import Pangenome
from ppanggolin.region import Region
from ppanggolin.formats import check_pangenome_info, write_pangenome, erase_pangenome
from ppanggolin.utils import restricted_float, mk_outdir


from itertools import chain



def compute_weigthed_metrics(G, partitions, metric, node_to_identical_rgp_count=None):
    
    if not node_to_identical_rgp_count:
        # count identical rgps is not taken into account
        node_to_identical_rgp_count = defaultdict(lambda:1)
    # else: each identical rgp is taken into account
    
    pairs_in_same_cluster = chain(*(combinations(cluster, 2) for cluster in partitions))

    # pairs_same_clstr_same_spot ==> true positive
    a = 0
    a_metric = 0

    # pairs_same_clstr_diff_spot ==> False Positive
    b = 0   
    b_metric = 0

    no_edge =0

    for rgp1, rgp2 in pairs_in_same_cluster:
        identical_edges = node_to_identical_rgp_count[rgp1] * node_to_identical_rgp_count[rgp2]
        assert identical_edges >= 1
        try:
            edge_data = G[rgp1][rgp2]
        except KeyError:
            no_edge+=1
            edge_data = {metric:0}
        if G.nodes[rgp1]['spot_id'] == G.nodes[rgp2]['spot_id']:
            a += 1 * identical_edges
            a_metric += edge_data[metric]* identical_edges
        else:
            # same cluster diff spot
            b += 1* identical_edges
            b_metric += (1 - edge_data[metric]) * identical_edges

            
    # if we count identical rgps then we should add pair of identical rgp in the True Positive metric (a)
    if node_to_identical_rgp_count:
        for rgp in G.nodes:
            identical_rgp_count = node_to_identical_rgp_count[rgp]

            a += int((identical_rgp_count ** 2 -  identical_rgp_count)/2)
            a_metric += int((identical_rgp_count ** 2 -  identical_rgp_count)/2)
    
    
    #logging.info(f'{no_edge=}')
    #logging.info(f'pairs_same_clstr_same_spot {a=} {a_metric=}')

    #logging.info(f'pairs_same_clstr_diff_spot {b=} {b_metric=}')

    pairs_in_diff_cluster = chain( * (product(cluster1, cluster2)  for cluster1, cluster2 in combinations(partitions, 2)))

    # pairs_diff_clstr_same_spot ==> False Negative
    c = 0
    c_metric = 0
    # pairs_diff_clstr_diff_spot ==> True Negative
    d = 0
    d_metric = 0

    no_edge =0

    for rgp1, rgp2 in pairs_in_diff_cluster:
        identical_edges = node_to_identical_rgp_count[rgp1] * node_to_identical_rgp_count[rgp2]
        
        try:
            edge_data = G[rgp1][rgp2]
        except KeyError:
            no_edge+=1
            edge_data = {metric:0}

        if G.nodes[rgp1]['spot_id'] == G.nodes[rgp2]['spot_id']:
            # same spot diff cluster
            c += 1 * identical_edges
            c_metric += edge_data[metric] * identical_edges
        else:
            # diff cluster and diff spot
            d += 1 * identical_edges
            d_metric += 1 - edge_data[metric] * identical_edges

    #logging.info(f'pairs_diff_clstr_same_spot {c=} {c_metric=}')

    #logging.info(f'pairs_diff_clstr_diff_spot {d=} {d_metric=}')
    ri = (a + d)/(a + b + c + d)

    weigthed_ri = (a_metric + d_metric)/(a_metric + d_metric + c_metric + b_metric )


    # sensitivity or recall =  True Positive / (True Positive + False Negative)
    sensitivity = a / ( a + c)  
    
    weighted_sensitivity = a_metric / (a_metric + c_metric)

    

    #  specificity = True Negative / (True negative + False Positive)  
    specificity = d / ( d + b)  
    
    weighted_specificity = d_metric / (d_metric + b_metric)

    
    #  precision = True positive / (True Positive + False Positive)

    try:
        precision = a/(a+b)
        weighted_precision = a_metric/(a_metric+b_metric)
    except ZeroDivisionError:

        precision = 0
        weighted_precision = 0 
    
    f1_score = (2*a) / (2*a + b + c)

    #logging.info(f"{ri=}")
    #logging.info(f"{weigthed_ri=}")

    metrics = {"sensitivity":sensitivity, 
               "weighted_sensitivity":weighted_sensitivity,
               "specificity":specificity,
               "weighted_specificity":weighted_specificity,
               "precision":precision,
               "weighted_precision":weighted_precision,
               "F1_score":f1_score,
                "weigthed_rand_index":weigthed_ri,
                "homemade_rand_index":ri,
               "True Negative": d, 
               "weighted True Negative":d_metric, 
               "False Positive":b,
               "weighted False Positive":b_metric,
               "True Positive":a,
               "weighted True Positive":a_metric,
               "False Negative":c,
               "weighted False Negative":c_metric
               }


    return metrics


def add_edges_to_identical_rgp_with_diff_spot(G):
    """
    """
    new_edges_count = 0
    for u, v, d in G.edges(data=True):
        # deal with rgp that have the flag different_spot
        if not 'different_spot' in d:
            continue
            
        identical_rgp, main_rgp = (u, v) if G.nodes[u]['identical_rgp'] else (v, u)

        # replicate all edges that connect main rgp to all identical rgps
        edges_to_add = []
        for i, connected_rgp in enumerate(G.neighbors(main_rgp)):
            edge_data = G[main_rgp][connected_rgp]
            edge_data['added_edge'] = True
            if identical_rgp != connected_rgp:
                edges_to_add.append((identical_rgp, connected_rgp, edge_data))

        G.add_edges_from(edges_to_add)
        new_edges_count += i+1

    logging.info(f"{new_edges_count} new edges, {G}")




def init_logging(verbose, debug):
    """Initialise logging."""
    if debug:
        level = logging.DEBUG
    elif verbose:
        level = logging.INFO
    else:
        level = logging.WARNING

    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s - %(message)s",
        datefmt="[%Y-%m-%d %H:%M:%S]",
    )

    logging.info("Program started")
    logging.info(
        f'command line: {" ".join(sys.argv)}',
    )


def get_node_to_cluster(partitions):
    node2cluster = {}
    for i, cluster in enumerate(partitions):
        node2cluster.update({node: i for node in cluster})

    return node2cluster       


def compute_grr_density(G, x_start=0, smoothing_parameter=0.1,  node_to_identical_rgp_count=None):
    
    if not node_to_identical_rgp_count:
        # count identical rgps is not taken into account
        node_to_identical_rgp_count = defaultdict(lambda:1)
        

    metric_counter = {"max_grr":defaultdict(int),
                    "min_grr":defaultdict(int),
                    "jaccard_index":defaultdict(int)}
    
    multiplier = 1000
    for u, v, d in G.edges(data=True):
        identical_edges = node_to_identical_rgp_count[u] * node_to_identical_rgp_count[v]
        for metric, grr_value2count  in metric_counter.items():
            grr = int(d[metric]  * multiplier)
            grr_value2count[grr] += 1 * identical_edges
            
            if node_to_identical_rgp_count:
                # add edges between identical rgp 
                # identical rgp has a grr of 1 so here the value of the multiplier
                grr_value2count[multiplier] += int((node_to_identical_rgp_count[u]**2 - node_to_identical_rgp_count[u])/2)
                grr_value2count[multiplier] += int((node_to_identical_rgp_count[v]**2 - node_to_identical_rgp_count[v])/2)
             

    df_count_list = []
    for metric, grr_value2count  in metric_counter.items():
        df_count_metric = pd.DataFrame([{"metric":metric, 'value':value, 'count':count} for value, count in grr_value2count.items()])
        df_count_list.append(df_count_metric)
    df_count = pd.concat(df_count_list)
    df_count["value"] = df_count["value"]/multiplier
    
    x_vals = np.linspace(int(x_start*multiplier),multiplier,multiplier) # Specifying the limits of our data
    density_df_samples = []
    #logging.info(f'Computing density for each metric')

    for metric, grr_value2count in metric_counter.items():
        
        counter = [grr_value2count[i+1] for i in range(int(x_start*multiplier), multiplier)]
        
        
        density = gaussian_kde(range(int(x_start*multiplier)+
                                     1,multiplier +1), weights = counter)
        density.covariance_factor = lambda : smoothing_parameter #Smoothing parameter
        density._compute_covariance()

        density_info = {"metric_value":x_vals, "Density":density(x_vals), 'metric':metric}


        df_density_metric = pd.DataFrame(data=density_info)
        df_density_metric["metric_value"] = df_density_metric["metric_value"]/multiplier
        density_df_samples.append(df_density_metric)

    df_density = pd.concat(density_df_samples)
    return df_density, df_count




def compute_cluster_metrics(G, clustering_weight, threshold, node_spot_vector, node_spot_vector_unmerged, node_to_identical_rgp_count, graph_nodes):
    """
    """
    unmerge_rgp = False

    edges_to_rm = [(u,v) for u,v,e in G.edges(data=True) if e[clustering_weight] < threshold]

    G_filt = nx.restricted_view(G, nodes=[], edges=edges_to_rm )

    #logging.info(f"Filtering graph edges with {edge_metric}<{threshold}: {grr_graph_filtered}")

    #print([d for u,v,d in G.edges(data=True)])
    if len(G_filt.edges) == 0:
        # if no egde then all node are in a separate cluster
        partitions = [{n} for n in G]
    else:
        partitions = nx.algorithms.community.louvain_communities(G_filt, weight=clustering_weight)

    node_to_clusterid = get_node_to_cluster(partitions)
    
    node_cluster_vector = []
    node_cluster_vector_unmerged = []
    for n in graph_nodes:
        node_cluster_vector.append(node_to_clusterid[n])
        node_cluster_vector_unmerged += [node_to_clusterid[n]] * node_to_identical_rgp_count[n]
        
    #print("metric computation", flush=True)
    ami = adjusted_mutual_info_score(node_spot_vector, node_cluster_vector)
    mi = normalized_mutual_info_score(node_spot_vector, node_cluster_vector)
    arand = adjusted_rand_score( node_spot_vector, node_cluster_vector)
    rand = rand_score(node_cluster_vector, node_spot_vector)
    
    all_rgp_ami = adjusted_mutual_info_score(node_spot_vector_unmerged, node_cluster_vector_unmerged)
    all_rgp_mi = normalized_mutual_info_score(node_spot_vector_unmerged, node_cluster_vector_unmerged)
    all_rgp_arand = adjusted_rand_score( node_spot_vector_unmerged, node_cluster_vector_unmerged)
    all_rgp_rand = rand_score(node_cluster_vector_unmerged, node_spot_vector_unmerged)

    w_metrics= compute_weigthed_metrics(G, partitions, clustering_weight)
    if unmerge_rgp:
        all_rgp_w_metrics= compute_weigthed_metrics(G, partitions, clustering_weight, node_to_identical_rgp_count)
        
        all_rgp_w_metrics = {f'all_rgp_{m}':v  for m, v in all_rgp_w_metrics.items() }
    

    #modularity_info = {}
    #for resolution in [0.5, 1, 2, 5]:
    #    modularity_info[f"modularity_resolution={resolution}"] = nx.algorithms.community.modularity(G, partitions, resolution=resolution, weight=None)
    #    modularity_info[f"weigthed_modularity_resolution={resolution}"] = nx.algorithms.community.modularity(G, partitions, weight=clustering_weight, resolution=resolution)

    #print("Done computing metrics", flush=True)
    metric = {"adjusted_mutual_info_score":ami, 
               "normalized_mutual_info_score":mi,
               "adjusted_rand_score":arand,
               "rand_score":rand,
              "all_rgp_adjusted_mutual_info_score":all_rgp_ami, 
               "all_rgp_normalized_mutual_info_score":all_rgp_mi,
               "all_rgp_adjusted_rand_score":all_rgp_arand,
               "all_rgp_rand_score":all_rgp_rand,
               "threshold":threshold,
               "metric":clustering_weight}
    metric.update(w_metrics)
    if unmerge_rgp:
        metric.update(all_rgp_w_metrics)
    #metric.update(modularity_info)
    
    return metric
    


def add_max_val(fig, df, y):
    df_max_list = []
    for metric in ['min_grr', 'max_grr', 'jaccard_index']:
        max_y = df.loc[(df['metric'] == metric), y].max()

        df_max = df.loc[(df['metric'] == metric) & (df[y] == max_y)]
        df_max_list.append(df_max)
    
    df_max = pd.concat(df_max_list)
    fig_max = px.scatter(df_max, x="threshold", y=y, color="metric",
                 )
    fig.add_traces(list(fig_max.select_traces()))  
      


def parse_arguments():
    """Parse script arguments."""
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter,)

    parser.add_argument("--graph", help="graphml file", required=True)

    parser.add_argument("-o", "--outdir", help="output dir", default='.')
    
    parser.add_argument("-c", "--cpu", help="number of cpus", type=int, default=1)

    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")

    parser.add_argument("--debug", help="active debug mode", action="store_true")

    parser.add_argument("--resume", help="active resume mode", action="store_true")

    
    args = parser.parse_args()
    return args


def main():
    "Orchestrate the execution of the program"

    args = parse_arguments()

    init_logging(args.verbose, args.debug)

    color_discrete_map={metric:color for metric, color in zip(["max_grr", 'min_grr', 'jaccard_index'], px.colors.qualitative.Vivid)}
    color_discrete_map

    graphml_file = args.graph
    outfile = os.path.join(args.outdir, "plot_metric.html")
    
    logging.info('Readind graph...')
    G = nx.read_graphml(graphml_file, node_type=int)
    logging.info(f'Graph: {G}')
    
    logging.info('Adding edges to identical rgp with different rgp...')
    add_edges_to_identical_rgp_with_diff_spot(G)
    
    logging.info(f"Graph before filtering {G}")
   
    invalid_nodes = [n for n, d in G.nodes(data=True) if d['spot_id'] == "No spot"]
    logging.info(f"{len(invalid_nodes)} nodes have no spot assciated")
    
    G.remove_nodes_from(invalid_nodes)

    logging.info(f"Graph after filtering {G}")

    node_to_identical_rgp_count = {n:d.get('identical_rgp_fam_and_spot', 1) for n,d in G.nodes(data=True)}
    df_density_derep, df_count_derep = compute_grr_density(G, 0, 0.09)
    df_density_derep['RGP dereplication'] = True
    df_count_derep['RGP dereplication'] = True

    df_density, df_count = compute_grr_density(G, 0, 0.09, node_to_identical_rgp_count)
    df_density['RGP dereplication'] = False
    df_count['RGP dereplication'] = False

    title="GRR distribution"

    df_density = pd.concat([df_density_derep, df_density])
    df_count = pd.concat([df_count_derep, df_count])
    #logging.info(f'Plot creation')
    fig = px.line(df_density, x='metric_value', y='Density', color='metric',  facet_row='RGP dereplication',
                title=title, color_discrete_map=color_discrete_map)
    fig.update_traces(line=dict(width=2), opacity=.7)

    #logging.info(f'Writing output html in {output}')


    fig_histo = px.histogram(df_count, x="value", color="metric", y="count", marginal=None, nbins=50, facet_row='RGP dereplication',
                    title=title,  facet_col="metric", color_discrete_map=color_discrete_map)

    #logging.info(f'Writing output html in {output}')
    #fig.write_html(output)
    with open(outfile, 'w') as f:
        f.write(fig.to_html(full_html=False, include_plotlyjs=True))
        f.write(fig_histo.to_html(full_html=False, include_plotlyjs=False))


    graph_nodes = list(G.nodes) 
    node_to_spotid = {n:d['spot_id'] for n, d in G.nodes(data=True)}
    node_spot_vector = [node_to_spotid[n] for n in graph_nodes]
    node_spot_vector_unmerged = []

    for n in graph_nodes:
        node_spot_vector_unmerged += [node_to_spotid[n]]*node_to_identical_rgp_count[n]


    thresholds = [i/100 for i in range(0,101,5)]
    clustering_weights = ["min_grr", "max_grr", 'jaccard_index']

    fct_args = ((G, w, t, node_spot_vector, node_spot_vector_unmerged, node_to_identical_rgp_count, graph_nodes) for t, w in product(thresholds, clustering_weights))


    cpu = 8
    metric_values = []
    with Pool(processes=cpu) as p:
        
        for result in tqdm(p.starmap(compute_cluster_metrics, fct_args), total=len(thresholds)*len(clustering_weights)):
        
            metric_values.append(result)

        
        
    figs= []
    figs_all_rgps = []
    df = pd.DataFrame(metric_values)
    for y in df.columns:
        if y in ['threshold', 'metric']:
            continue
        
        fig = px.line(df, x="threshold", y=y, color="metric", title=f"{y}", color_discrete_map=color_discrete_map)
        if df[y].max() <= 1:
            fig.update_layout(yaxis_range=[0,1])
        add_max_val(fig, df, y)
        if "all_rgp" in y:
            figs_all_rgps.append(fig)
        else:
            figs.append(fig)
    
    df['False Positive rate'] = 1 - df['specificity']
    df['True Positive rate'] = df['sensitivity']
    fig_roc = px.line(df, y='True Positive rate', x='False Positive rate', color="metric", title=f"ROC curve", text="threshold", color_discrete_map=color_discrete_map)
    fig_roc.update_layout(yaxis_range=[0,1], xaxis_range=[0,1])
    fig_roc.update_traces(textposition="bottom right")

    with open(outfile, 'a') as f:
        for i, fig in enumerate(figs):
            f.write(fig.to_html(full_html=False, include_plotlyjs=False))
        f.write(fig_roc.to_html(full_html=False, include_plotlyjs=False)) 

        
    # df['all_rgp_False Positive rate'] = 1 - df['all_rgp_specificity']
    # df['all_rgp_True Positive rate'] = df['all_rgp_sensitivity']
    # fig_roc = px.line(df, y='all_rgp_True Positive rate', x='all_rgp_False Positive rate', color="metric", title=f"ROC curve", text="threshold", color_discrete_map=color_discrete_map)
    # fig_roc.update_layout(yaxis_range=[0,1], xaxis_range=[0,1])
    # fig_roc.update_traces(textposition="bottom right")
        

    logging.info(f'Plots are written in {outfile}')


    outfile_allrgps =  os.path.join(args.outdir, "all_rgp_plot_metric.html")
    with open(outfile_allrgps, 'w') as f:
        for i, fig in enumerate(figs_all_rgps):
            include = True if i == 0 else False

            
            f.write(fig.to_html(full_html=False, include_plotlyjs=include))

        # f.write(fig_roc.to_html(full_html=False, include_plotlyjs=False)) 
    logging.info(f'Plots are written in {outfile_allrgps}')

# If this script is run from the command line then call the main function.
if __name__ == "__main__":
    main()
