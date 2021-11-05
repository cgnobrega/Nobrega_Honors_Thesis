# Nobrega Honors Thesis
Repository for the code &amp; data involved in the Honors Thesis project of Celeste Nobrega

## Methods
This study will be working towards identifying specific cell type responses in cortisol carrier adults bulk RNA seq data by using the gene lists for specific cell types from Hou et. al. scRNA seq data, and pulling those same genes out of the Hartig et. al. bulk RNA seq data, and seeing if those genes are differentially expressed between cortisol and control, suggesting that specific cell types are differentially impacted by cortisol treatment.

The first goal of the project will be to reproduce the results from the Hou et. al. (2020) single cell analysis to fully understand any assumptions made during that study that may be imperative to the interpretation of our results. This will be done by using the data available on the NCBI GEO database and the analysis tools available on MDI Biological Laboratory’s remote cluster.  Then, the next step will be to use RNA-Sieve or another tool for integrating scRNA seq with bulk RNA seq to combine the data between Hou et. al. (2020) and Hartig et. al (2016) and identify candidate cell types that are most affected in the adults that were grown in cortisol based on gene expression data. 

As a proof of concept, the first step will be to do a simple comparison of the existing gene lists available from the Coffman lab and the supplementary tables in Hou et. al. (2020) to find shared genes between the two studies and further look at their expression values to make any preliminary observations and direct the next steps in the project. 

Hartig, E. I., Zhu, S., King, B. L., & Coffman, J. A. (2016). Cortisol-treated zebrafish embryos develop into pro-inflammatory adults with aberrant immune gene regulation. *Biology open*, 5(8), 1134-1141.

Hou, Y., Lee, H. J., Chen, Y., Ge, J., Osman, F. O. I., McAdow, A. R., ... & Wang, T. (2020). Cellular diversity of the regenerating caudal fin. *Science advances*, 6(33), eaba2084.
