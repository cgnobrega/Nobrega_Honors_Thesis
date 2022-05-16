# Data

```
data
│   README.md
│   Danio_rerio.GRCz10.91.gene.counts.txt
│   Hartig_0dpa_ctrl_vs_treated.csv
│   Hartig_0dpa_vs_2dpa_ctrl.csv
│   Hartig_0dpa_vs_2dpa_interaction.csv
│   Hartig_0dpa_vs_2dpa_treated.csv
│   Hartig_2dpa_ctrl_vs_treated.csv
│   Hartig_4dpa_ctrl_vs_treated.csv
│   Hou_gene_list_full.csv
│   Hou_gene_list_top50.csv
└───Hartig_Fin_Data_raw/
│   │   Fin_0dpa_CTRL_vs_0dpa_Treated.xlsx
│   │   Fin_0dpa_CTRL_vs_2dpa_CTRL.xlsx
│   │   Fin_0dpa_TREATED_vs_2dpa_TREATED.xlsx
│   │   Fin_2dpa_CTRL_vs_2dpa_Treated.xlsx
│   │   Fin_4dpa_CTRL_vs_4dpa_Treated.xlsx
│   │   Fin_interaction_0dpa_vs_2dpa.xlsx
└───Hou_Data_raw
│   │   GSE137971_cells.csv
│   │   SupplementaryTableS2.xlsx
│   │   aba2084_tables_s1_to_s7.xlsx
└───RNASieve_input_2dpa
│   │   bulk_count_2dpa.csv.gz
│   │   bulk_design_2dpa.csv.gz
│   │   single_count_2dpa.csv.gz
│   │   single_design_2dpa.csv.gz
└───RNASieve_output
│   │   rnasieve_jcoffman007.csv
│   │   rnasieve_jcoffman007_0dpa.csv
│   │   rnasieve_jcoffman007_2dpa.csv
│   │   rnasieve_jcoffman007_4dpa.csv
```

Preliminary data storage

### Hartig_Fine_Data_raw

Excel spreadsheets of the tailfin RNA-seq analysis that identified differentially expressed genes  of the Hartig 2016 paper. The differentially expressed genes are highlighted, green for genes downregulated in cort-treated fish, red for upregulated genes.

### Hou_Data_raw

`aba2084_tables_s1_to_s7.xlsx` : Excel spreadsheet of supplemental tables S1 through S7 from Hou 2020. 

`GSE137971_cells.csv` : assignment of cell identifiers to clusters

`GSM4095[393-400]_[samp1,samp2,1dpa,2dpa,4dpa]_DGEmatrix.csv.gz` : expression matrices from Hou data, pulled from NCBI GEO database supplementary file. 

`SupplementaryTableS2.xlsx` : Excel spreadsheet of full supplementary table 2.

### Hartig_gene_list.csv

`Hartig_gene_list.csv` : Initial gene list of differentially expressed genes. Add factor column titled `DiffExp` with values [`upregulated`, `downregulated`, `none`] to reflect the highlighting in `Fin_interaction_0dpa_vs_2dpa.xlsx`.

`Hartig_[0,2,4]dpa_ctrl_vs_treated` : gene lists of differentially expressed genes from Hartig 2016 paper for each dpa. 

`Hartig_0dpa_vs_2dpa_[ctrl,treated].csv` : gene list of differentially expressed genes for the interaction between 0 and 2 dpa. 

### Hou_gene_list.csv

`Hou_gene_list.csv` : Supplementary Table S2 Genes differentially enriched in major cell types **(top 50)**. 

`Hou_gene_list_full.csv` : **Full** Supplementary Table S2
