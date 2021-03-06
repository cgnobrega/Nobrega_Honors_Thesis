# Code

```
code
│   README.md
│   Running RNA-seq pipeline using singularity image.md    
│   gene_lists_parse.ipynb
│   matrix_cleaning.ipynb
│   set_matching.ipynb
│   subsetting_matrices.ipynb
│   all_box_plots.png
│   all_std_mad_plots.png
│   all_violin_plots.png
└───RNASieve_scripts
│   │   rnasieve_jcoffman007_0dpa.ipynb
│   │   rnasieve_jcoffman007_2dpa.ipynb
│   │   rnasieve_jcoffman007_4dpa.ipynb
│   │   rnasieve_jcoffman007.ipynb
└───singularity_param_files
│   │   kallisto_index_parameter_template.txt
│   │   kallisto_tximport_to_genes.R
│   │   v_Danio_rerio.GRCz10.91.tx2gene.txt
│   └───fastqc_param
│   │   │  # param files for each bulk sample
│   └───trimgalore_param
│   │   │  # param files for each bulk sample
│   └───kallisto_quant_param
│   │   │  # param files for each bulk sample
└───boxplot
│   │   # 6 boxplot .png files
└───violin_plot
│   │   # 6 violin plot .png files
└───std_mad_barplot
│   │   # 6 bar plot .png files
└───summary_tables
│   │   # 6 data table .png files

```

### Boxplot

Box plots for each of the samples grouped by cell type 
cluster. Values plotted are the log fold-change according
to the Hartig data. Data points are also shown, marked 
with differing colors to show different expression 
levels. Red denotes upregulated genes, blue represents 
downregulated genes, and grey denotes genes that are not
differentially expressed. 

### std_mad_barplot

Standard deviation (std) and mean absolute deviation (MAD) plotted as 
bar plots for each sample, grouped by cell type cluster. 
This plot is to show two different measures of deviation. 

### summary_tables

Tables for each sample that contain summary statistics about
the log fold change value (from the Hartig paper). The 
information is grouped into columns based on cell type
cluster. The summary statistics included are: count (number
of genes in this cluster), mean (avg logFC), std (standard
deviation), min, max, quartiles, MAD (mean absolute deviation). 

### violin_plot

Violin plots for each of the samples grouped by cell type 
cluster. Values plotted are the log fold-change according
to the Hartig data. Data points are also shown, marked 
with differing colors to show different expression 
levels. Red denotes upregulated genes, blue represents 
downregulated genes, and grey denotes genes that are not
differentially expressed. 

### all_[box,std_mad,violin]_plots.png

Each file shows a combined image of all the plots for the 
individual samples. This is to compare plots across samples. 
We can also see a pattern between 0, 2, and 4 dpa and the
deviation. Unclear if this pattern has significance or meaning. 

### gene_lists_parse.ipynb

Jupyter Notebook used to generate these plots, and to perform
all preliminary data comparisons and analyses. Mainly used for
set matching to create a combined list of genes that were found in 
both the Hartig and Hou studies. 

### Running RNA-seq pipeline using singularity image.md

A file that describes the process of how RNA-seq analysis was
run. Using singularity images, on an Amazon AWS machine, and bash
scripting to automate the process.