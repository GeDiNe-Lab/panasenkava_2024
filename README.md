# General notes
"Condition1 VS Condition2" in DE means "Expression in **Condition1** compare to **Condition2**"
A gene with a FoldChange > 0 is **MORE** expressed in **Condition1** compare to **Condition2**
A gene with a FoldChange < 0 is **LESS** expressed in **Condition1** compare to **Condition2**

Differentially expressed genes used for Heatmap are not filtered out if they have no Symbol associated to the Ensemble ID.

Log2FoldChange threshold is 2 for the Volcano plots in order to vizualize genes that are "highly" differentially expressed, otherwise it is 1.

# Figure 1

## Known literature marker
### F1_1_marker_HM
Heatmap of the known marker genes.

### F1_1_marker_HM
Heatmap of the known genes grouped by Pluripotency markers, Forebrain markers, Ventral markers, Dorsal markers, & Posterior markers

## PCA
### F1_2_PCA
PCA on top 3000 most variable genes and PCs variation percentage (F1_2a_percentVar.png)
PC_covariate_correlation indicates the correlation between PC and covariates, the pvalue of the ANOVA associated are in tables/Figure_1/F1_PC_covariate_ANOVA.csv

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
PCA  of the top 3000 most variable genes colored by lineage (line), symbols shows the differentiation days
PC_covariate_correlation indicates the correlation between PC and covariates, the pvalue of the ANOVA associated are in tables/Figure_2A/F2_PC_covariate_ANOVA_3000genes.csv

## Kinetic Heatmap
### F2A_DE_HM
Heatmap of the top 1000 genes correlated with PC1 and PC2 (absolute Pearson correlation) that or not in the top 1000 genes correlated with PC3.
PC1 is associated with kinetic, PC2 is associated with dorso/ventral patterning, PC3 is associated with lineage difference

There is 2 level of clustering subclustering has been made using gene from each cluster of the first level respectively

Associated tables : genes_cluster.csv

### F2A heatmap functional enrichment
(First level of clustering : clusters 1,2,3, and 4) top 15 GO biological process term for each DE genes clusters

Associated tables : GO_enrichment_cluster_n.csv

Also did the same for LON only and WTC only clustering, images adn tables have _LON or _WTC at the end accordingly

## Other Differential expression tables
Highly differentilly expressed genes : absolute log2FC >= 2 and pvalue < 0.01

### DE between day
tables : DEGs_dayX_Y_dorsal.csv and DEGs_dayX_Y_ventral.csv
images: in the volcano_plots folder

### DE between type
tables : DEGs_DV_dayX.csv
iùages: in the volcano_plots folder

# Figure 3
## PCA
### F3_PCA
PCA of the top 3000 most variable genes colored by ventral/cyclopamine type, shape relative to cyclopamine dosage
PC_covariate_correlation indicates the correlation between PC and covariates, the pvalue of the ANOVA associated are in tables/Figure_3/F3_PC_covariate_ANOVA.csv 

## Table recapitulating DE, WGCNA and heatmap clustering
the file cyclo_genes_df.csv recapitulates DE results for the 3 following contrast :
HvsN : High Cylopamin dosage (1/0.5) VS No cyclopamin (ventral samples)
HvsL : High Cylopamin dosage (1/0.5) VS Low Cyclopamin dosage (0.25/0.125)
LvsN : Low Cylopamin dosage (0.25/0.125) VS No cyclopamin (ventral samples)

Column H,I, and J indicates if the genes pass the DE thresholds (log2FC>=1 & padj<0.01) for the 3 contrasts

The WGCNA column indicates if the genes is in the SHH WGCNA cluster.
The cluster and the sub_cluster column are related to the genes cluster on the heatmap.

### F3_cyclo_genes_HM
Heatmap of the genes within the same WGCNA cluster than SHH.
There is 2 level of clustering

### GO_enrichment_cluster_n
(First level of clustering : clusters 1 & 2) top 15 GO biological process term for each DE genes clusters

Associated tables : GO_enrichment_cluster_n.csv

### Volcano plot
DE genes information can be found in tables/Figure_3/cyclo_genes_df.csv
Volcano plots for highly DE genes between the different dose of cyclopamin (abs(log2FC) >= 2, padj < 0.01)

## STRING network
Networks built by passing the top 500 most correlated genes with SHH that are also in the same WGCNA cluster to STRING, as STRING give information about protein relationship it is not taking in account non-coding RNA

# Figure 4
## PCA_type
PCA of the top 3000 most variable genes on the WTC CRISPR lineage for vAN only
PC_covariate_correlation indicates the correlation between PC and covariates, the pvalue of the ANOVA associated are in tables/Figure_4/
F4_PC_covariate_vAN_ANOVA.csv


## volcano plots
Volcano plots for the DEGs with 3 group of samples : vAN
contrast are :
 * control_vs_homo
 * control_vs_het
 * het_vs_homo

## Markers barplots
Barplots of the scaled normalized expression for 4 marker genes for vAN samples
(values are the normalized readcounts rescaled so the minimum normalized value (accounting for 0 read) is equal to 0)

## Coexpressed genes barplots
Barplots of the scaled normalized expression for genes within the SHH WGCNA cluster for the vAN samples
