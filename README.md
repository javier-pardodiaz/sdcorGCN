# sdcorGCN
Robust gene coexpression networks using signed distance correlation

Preprint:
Supplementary Infomration:

Use the generate_sdcor_GCN.R script to construct GCN using signed distance correlation (network NS(dS)). This file also contains the code to generate GCN using Pearson correlations (network NP(dP)). The optimal edge densities for both correlations is estimated using COGENT (https://github.com/lbozhilova/COGENT). It also saves the similarity and scores for each tested threshold.

Use the generate_matching_CGN.R script to construct GCN with a different edge density (NS(dP) and NP(dS). This file also contains the code to analyse the GCN: number of edges, size of largest connected component, global clustering coefficient.

The script spearman_GCN.R constructs the Spearman GCNs NR(dS) and NR(dP).

To perform the biological evaluation of the networks using STRING, use STRING_evaluation.R. The required format for the input matrices with the confidence scores is described here. To obtain the matrices, use combine_subscores.py and follow the instructions in the Supplementary Information (Section 6).

The data folder contains the expression matrices ready to be used and the matrices with the confidence scores form STRING for each organism.

The results can be found in 
