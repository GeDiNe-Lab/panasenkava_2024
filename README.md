# Figure 1

## Known literature marker
### F1_1_inter_marker_HM
Heatmap of the known intermediate marker genes.

### F1_1_marker_HM
Heatmap of the known genes grouped by Pluripotency markers, Forebrain markers, Ventral markers, Dorsal markers, & Posterior markers

## PCA
### F1_2_PCA
PCA on all genes

## Dorso/ventral differentially expressed genes
### F1_3_DE_HM
Heatmap of the differentially expressed genes (absolute log2FC >= 1, pvalue < 0.01) with 3 levels of clustering

Associated tables : dorsal_VS_ventral_DEGs.csv

### F1_DE_GO
(First level of clustering : clusters 1 & 2) top 15 GO biological process term for each DE genes clusters

Associated tables : GO_enrichment_cluster_n.csv

### F1_4_LON_WTC_DE
(log2FC : 0.5, 1, 1.5, 2 and pvalue < 0.01) Venn diagramm of the overlap of dorso/ventral diffenrentially expressed genes between LON and WTC lineage

Associated tables : dorsal_VS_ventral_DEGs_LON.csv & dorsal_VS_ventral_DEGs_WTC.csv


# Figure 2

## PCA
### F2A_1_PCA
PCA colored by dorso/ventral patterning (type) or lineage (line)

## Kinetic Heatmap
### F2A_DE_HM
Heatmap of the differentially expressed genes between dorsal and ventral, (absolute log2FC >= 1, pvalue < 0.01). There is 2 level of clustering subclustering has been made using gene from each cluster of the first level respectively

Associated tables : DEG_kinetic_dorsal_VS_ventral.csv

### F1_DE_GO
(First level of clustering : clusters 1 to 5) top 15 GO biological process term for each DE genes clusters

Associated tables : GO_enrichment_cluster_n.csv

## Other Differential expression tables
(absolute log2FC >= 1 and pvalue < 0.01
### DE between day
tables : DEGs_dayX_Y_dorsal.csv and DEGs_dayX_Y_ventral.csv
images: in the volcano_plots folder

### DE between type
tables : DEGs_DV_dayX.csv
iùages: in the volcano_plots folder

# Figure 3
## Differentially expressed genes between no cylo and high clyclo dosage
### F3_cyclo_DE_HM
Heatmap of the differentially expressed genes between no cyclopamin and high dosage of cyclopamin (absolute log2FC >= 1, pvalue < 0.01).
There is 3 level of clustering : 

Associated tables : DEG_cyclo_vAN_VS_high.csv

### GO_enrichment_cluster_n
(First level of clustering : clusters 1 & 2) top 15 GO biological process term for each DE genes clusters

Associated tables : GO_enrichment_cluster_n.csv

## Coexpressed genes
### Cytoscape_SHH
genes retained in the coexpression analysis with either 0.8 or 0.85 absolute pearson correlation with one the genes used (SHH, NKX2.1, PAX6), tehre is an image for the genes positively co-expressed with SHH and one for the genes negatively co-expressed with SHH, color are : gold -> already known genes, red -> gene selected through literature curation and potentially tested in mice.




