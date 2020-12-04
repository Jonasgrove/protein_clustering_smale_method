# protein_clustering_smale_method
Implementation of Stephen Smale's method for calculating protein distances from sequences. "Towards a Mathematical Foundation of Immunology and Amino Acid Chains"

This program uses the distances calculated using Smale's method to construct a minimum spanning tree using Kruskal's algorithm, in which nodes are proteins and edges are protein distances. K clusters are fored by iteratively removing the most expensive (greatest distance) K-1 edges. This clusters the chosen proteins based on their sequence similarity, as definied by Smale's method.
