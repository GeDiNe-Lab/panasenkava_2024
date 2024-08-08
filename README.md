# General notes
"Condition1 VS Condition2" in DE means "Expression in **Condition1** compare to **Condition2**"
A gene with a FoldChange > 0 is **MORE** expressed in **Condition1** compare to **Condition2**
A gene with a FoldChange < 0 is **LESS** expressed in **Condition1** compare to **Condition2**

# Figure 1

## Known literature marker
### F1_1_inter_marker_HM
Heatmap of the known intermediate marker genes.

### F1_1_marker_HM
Heatmap of the known genes grouped by Pluripotency markers, Forebrain markers, Ventral markers, Dorsal markers, & Posterior markers

## PCA
### F1_2_PCA
PCA on all genes and PCs variation percentage (F1_2a_percentVar.png)

## Dorso/ventral differentially expressed genes
### F1_3_DE_HM
Heatmap of the differentially expressed genes (absolute log2FC >= 1, pvalue < 0.01) with 3 levels of clustering

Associated tables : dorsal_VS_ventral_DEGs.csv

### F1_DE_GO
(First level of clustering : clusters 1 & 2) top 15 GO biological process term for each DE genes clusters

Associated tables : GO_enrichment_cluster_n.csv

# Figure 2

## PCA
### F2A_1_PCA
PCA colored by dorso/ventral patterning (type) or lineage (line), with all or the top 500 variable genes, for PC1/2 or PC2/3

## Kinetic Heatmap
### F2A_DE_HM
Heatmap of the top 1000 genes correlated with PC1 and PC2 (absolute Pearson correlation) that or not in the top 1000 genes correlated with PC3.
PC1 is associated with kinetic, PC2 is associated with dorso/ventral patterning, PC3 is associated with lineage difference

There is 2 level of clustering subclustering has been made using gene from each cluster of the first level respectively

Associated tables : genes_cluster.csv

### F2A_DE_GO
(First level of clustering : clusters 1 to 5) top 15 GO biological process term for each DE genes clusters

Associated tables : GO_enrichment_cluster_n.csv

## Other Differential expression tables
Highly differentilly expressed genes : absolute log2FC >= 2 and pvalue < 0.01

### DE between day
tables : DEGs_dayX_Y_dorsal.csv and DEGs_dayX_Y_ventral.csv
images: in the volcano_plots folder

### DE between type
tables : DEGs_DV_dayX.csv
iùages: in the volcano_plots folder

# Figure 3
## Genes correlated with cyclopamin dosage
Genes correlated with cyclopamin dosage value in uL (absolute Pearson correlation >= 0.5
### F3_cyclo_genes_HM
Heatmap of the genes correlated with cyclopamin dosage.
There is 2 level of clustering : 

Associated tables : DEG_cyclo_vAN_VS_high.csv

### GO_enrichment_cluster_n
(First level of clustering : clusters 1 & 2) top 15 GO biological process term for each DE genes clusters

Associated tables : GO_enrichment_cluster_n.csv

### Volcano plot
Volcano plots for highly DE genes between the different dose of cyclopamin (abs(log2FC) >= 2, padj < 0.01)

## Coexpressed genes
### Cytoscape_SHH
genes retained in the coexpression analysis with either 0.8 or 0.85 absolute pearson correlation with one the genes used (SHH, NKX2.1, PAX6)

# Figure 4
## PCA_type
PCA on the WTC CRISPR lineage for vAN,dAN, vAN+dAN

## volcano plots
Volcano plots for the DEGs with 3 group of samples : vAN, dAN and vAN+dAN
contrast are :
 * control_vs_homo
 * control_vs_het
 * het_vs_homo

## Markers barplots
Barplots of the scaled normalized expression for 4 marker genes for vAN and dAN samples
(values are the normalized readcounts rescaled so the minimum normalized value (accounting for 0 read) is equal to 0)

## Coexpressed genes barplots
Barplots of the scaled normalized expression for genes found through co-expression analysis with both dAN and vAN samples
