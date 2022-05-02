#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
#!pip install anndata --target D:/usr/local/lib/python3.8/dist-packages
import anndata as ad
from collections import defaultdict
from rnasieve.preprocessing import model_from_raw_counts
print(ad.__version__)


# This file uses guidance from https://github.com/theislab/anndata-tutorials/blob/master/getting-started.ipynb
#  to get the jcoffman008 data as an AnnData structure so that it can be used for RNAsieve as in https://github.com/songlab-cal/rna-sieve
# 

# In[5]:


##-----loading in jcoffman007 expression and design data -----
#load in jcoffman007 count data (raw counts)
widejcoffman007_count = pd.read_csv("all_bulk_count_file.csv.gz", compression="gzip")
#changing GeneID_genename format to just GeneID
# widejcoffman007_count["gene_id"]=widejcoffman007_count["gene_id"].str[:18]
#load in design file. note that each row is a sample that matches up with the "data" dataframe  
jcoffman007_design= pd.read_csv("all_bulk_design_file.csv.gz", compression="gzip",index_col=0)


##-----loading in atlas expression and design data----- 
#uploading the expression data for the reference (atlas)
wideatlas_count = pd.read_csv("fill_sorted_matched_Hou_shared_gene_counts.csv.gz", compression="gzip")
#changing GeneID|genename format to just GeneID
# wideatlas_count["gene"]=wideatlas_count["gene"].str[:18]
#uploading the meta/"obs" for the reference (atlas)
''' 
atlas_meta_alldpf = pd.read_csv("meta.tsv", sep="\t",index_col=0)
#restricting meta to only 5dpf 
sample_names_5dpf= ['5a_5dpf', '5b_5dpf']
atlas_meta= atlas_meta_alldpf[atlas_meta_alldpf['sample_name'].isin(sample_names_5dpf)]
'''
atlas_design = pd.read_csv("sorted_labeled_GSE137971_cells.csv.gz", compression="gzip",index_col=0)


# In[6]:


#altering both expression datasets to contain identical sets of genes 
jcoffman007_genes= widejcoffman007_count["genes"].tolist()
atlas_genes= wideatlas_count["genes"].tolist()
genes_in_both= list(set(jcoffman007_genes) & set(atlas_genes))
print("there are",len(genes_in_both),"genes that are in both the atlas and jcoffman007 count data") 
widejcoffman007_count=widejcoffman007_count[widejcoffman007_count["genes"].isin(genes_in_both)]
wideatlas_count=wideatlas_count[wideatlas_count["genes"].isin(genes_in_both)]
#sorting to be in same order
widejcoffman007_count=widejcoffman007_count.sort_values(by=['genes'])
wideatlas_count=wideatlas_count.sort_values(by=['genes'])


# In[7]:


#------making annData structure that holds the counts and meta for jcoffman007-----
#convert to having genes as columns, samples as rows (first have to make the gene ID the row label)
widejcoffman007_count= widejcoffman007_count.set_index("genes")
jcoffman007_count= widejcoffman007_count.transpose()
jcoffman007 = ad.AnnData(jcoffman007_count, obs=jcoffman007_design)

#-----making annData structure that holds the counts and meta for atlas----- 
wideatlas_count=wideatlas_count.set_index("genes")
atlas_count= wideatlas_count.transpose()
atlas = ad.AnnData(atlas_count, obs=atlas_design)


# In[8]:


widejcoffman007_count


# In[9]:


#modeling off of the raw counts prep section from example.ipynb from song lab github

# grouping the cells in the reference data by cluster 
print('Aggregating by cluster...')
counts_by_cluster = {}
for i in range(len(atlas)):
    sc = atlas[i]
    if len(sc.obs['major.cl']) == 0:
        continue
    cell_cluster = sc.obs['major.cl'][0]
    if cell_cluster not in counts_by_cluster:
        counts_by_cluster[cell_cluster] = np.empty((sc.X.shape[1], 0), dtype=np.float32)
    counts_by_cluster[cell_cluster] = np.hstack(
        (counts_by_cluster[cell_cluster], sc.X.toarray().reshape(-1, 1)))
    

# Bulk prep
print('Aggregating bulks by name...')
G = jcoffman007.n_vars
bulk_by_time = defaultdict(list)
for i in range(len(jcoffman007)):
    bulk = jcoffman007[i]
    if len(bulk.obs['name']) == 0:
        continue
    time = bulk.obs['name'][0]
    bulk_by_time[time].append(bulk.X.toarray().reshape(-1, 1))

bulk_labels = []
psis = np.empty((G, 0), dtype=np.float32)
for name in sorted(bulk_by_time.keys()):
    bulks = bulk_by_time[name]
    for i in range(len(bulks)):
        bulk_labels.append("{} name, subject {}".format(name, i))
        psis = np.hstack((psis, bulks[i]))
        
print('Done!')


# In[10]:


counts_by_cluster['Intermediate Epithelial'].shape


# In[11]:


model, cleaned_psis = model_from_raw_counts(counts_by_cluster, psis[:, :18])


# In[12]:


cleaned_psis


# In[13]:


output_table= model.predict(cleaned_psis) #this takes a minute ... 
output_table.to_csv('resorted_rnasieve_jcoffman007.csv')


# In[14]:


# In this example, the intervals at a significance level of 0.05 indicate the estimate is poor.
# We set sig=0.9999 for the sake of visualization.
model.compute_marginal_confidence_intervals(sig=0.2)


# In[15]:


model.plot_proportions('bar').properties(title="Bar visualization ")


# In[16]:


model.plot_proportions('stacked').properties(title="Stacked visualization")


# In[ ]:




