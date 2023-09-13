
To cluster RGPs (Regions of Genome Plasticity) based on their gene families, you can use the command `panggolin rgp_cluster`.
The panggolin rgp_cluster command performs the following steps to cluster RGPs (Regions of Genome Plasticity) based on their gene families:

1. Calculation of GRR (Gene Repertoire Relatedness): The command calculates the GRR values for all pairs of RGPs. The GRR metric evaluates the similarity between two RGPs by assessing their shared gene families.
2. Graph Construction: The command constructs a graph representation of the RGPs, where each RGP is represented as a node in the graph. The edges between the nodes are weighted using the GRR values, indicating the strength of the relationship between the RGPs.
3. Filtering GRR Values: GRR values below the `--grr_cutoff` threshold (default 0.8) are filtered out to remove noise from the analysis.
4. Louvain Communities Clustering: The Louvain communities clustering algorithm is then applied to the graph. This algorithm identifies clusters of RGPs with similar gene family relationships.

There are three modes available for calculating the GRR value: `min_grr`, `max_grr`, or `incomplete_aware_grr`.
- `min_grr` mode: This mode computes the number of gene families shared between two RGPs and divides it by the smaller number of gene families among the two RGPs.
- `max_grr` mode: In this mode, the number of gene families shared between two RGPs is calculated and divided by the larger number of gene families among the two RGPs.
- `incomplete_aware_grr` (default) mode: If at least one RGP is considered incomplete, which typically happens when it is located at the border of a contig, the `min_grr` mode is used. Otherwise, the `max_grr` mode is applied. This mode is useful to correctly cluster incomplete RGPs.


The resulting RGP clusters are stored in a tsv file with the folowing columns:

| column  | description                  |
|---------|------------------------------|
| RGP     | The unique region identifier |
| cluster | The cluster id of the RGP    |
| spot_id    | the spot ID of the RGP       |


The command also generates an RGP graph in the gexf format, which can be utilized to explore the RGP clusters along with their spots of insertion. In this graph identical RGPs with the same family content and with the same spot are merged into a single node to simplify the graph representation. This feature can be disable with the parameter `--no_identical_rgp_merging`.
