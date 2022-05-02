#!/usr/bin/env python
# coding: utf-8

# # Re-doing set matching & plot generation
# This is all code from last semester, just being re-documented in a cleaner & more efficient way

# In[1]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# In[2]:


# create a data frame for the Hou gene list

# read in full marker file as Hou
hou = pd.read_csv("Hou_gene_list_full.csv")
hou = hou.loc[:,['cluster','gene', 'avg_logFC']]


# In[3]:


#create a data frame for the Hartig gene list 0dpa ctrl vs treated
hartig_0dpa = pd.read_csv("Hartig_0dpa_ctrl_vs_treated.csv")


# In[4]:


#create a data frame for the Hartig gene list 2dpa ctrl vs treated
hartig_2dpa = pd.read_csv("Hartig_2dpa_ctrl_vs_treated.csv")


# In[5]:


#create a data frame for the Hartig gene list 4dpa ctrl vs treated
hartig_4dpa = pd.read_csv("Hartig_4dpa_ctrl_vs_treated.csv")


# In[6]:


#create a data frame for the Hartig gene list 0dpa vs 2dpa ctrl
hartig_0v2_ctrl = pd.read_csv("Hartig_0dpa_vs_2dpa_ctrl.csv")


# In[7]:


#create a data frame for the Hartig gene list 0dpa vs 2dpa treat
hartig_0v2_treat = pd.read_csv("Hartig_0dpa_vs_2dpa_treated.csv")


# In[8]:


#create a data frame for the Hartig gene list 0dpa vs 2dpa interaction
hartig_0v2_inter = pd.read_csv("Hartig_gene_list.csv")


# In[9]:


# create a new data frame that contains only those rows that have matching values in both of the original data frames
# based on the 'gene' column in Hou and 'symbol' column in Hartig
merged_gene_list_0dpa = pd.merge(left=hou, right=hartig_0dpa, left_on='gene', right_on='Symbol')

# re-order & re-name columns as necessary
merged_gene_list_0dpa = merged_gene_list_0dpa[['gene', 'Symbol', 'GeneID', 'cluster', 'DiffExp', 'avg_logFC', 'logFC']]
merged_gene_list_0dpa.columns = ["gene_hou", "gene_hartig", "GeneID", "cluster", "DiffExp", "logFC_hou", "logFC_hartig"]


# In[10]:


# create a new data frame that contains only those rows that have matching values in both of the original data frames
# based on the 'gene' column in Hou and 'symbol' column in Hartig
merged_gene_list_2dpa = pd.merge(left=hou, right=hartig_2dpa, left_on='gene', right_on='Symbol')

# re-order & re-name columns as necessary
merged_gene_list_2dpa = merged_gene_list_2dpa[['gene', 'Symbol', 'GeneID', 'cluster', 'DiffExp', 'avg_logFC', 'logFC']]
merged_gene_list_2dpa.columns = ["gene_hou", "gene_hartig", "GeneID", "cluster", "DiffExp", "logFC_hou", "logFC_hartig"]


# In[11]:


# create a new data frame that contains only those rows that have matching values in both of the original data frames
# based on the 'gene' column in Hou and 'symbol' column in Hartig
merged_gene_list_4dpa = pd.merge(left=hou, right=hartig_4dpa, left_on='gene', right_on='Symbol')

# re-order & re-name columns as necessary
merged_gene_list_4dpa = merged_gene_list_4dpa[['gene', 'Symbol', 'GeneID', 'cluster', 'DiffExp', 'avg_logFC', 'logFC']]
merged_gene_list_4dpa.columns = ["gene_hou", "gene_hartig", "GeneID", "cluster", "DiffExp", "logFC_hou", "logFC_hartig"]


# In[12]:


# create a new data frame that contains only those rows that have matching values in both of the original data frames
# based on the 'gene' column in Hou and 'symbol' column in Hartig
merged_gene_list_0v2_ctrl = pd.merge(left=hou, right=hartig_0v2_ctrl, left_on='gene', right_on='Symbol')

# re-order & re-name columns as necessary
merged_gene_list_0v2_ctrl = merged_gene_list_0v2_ctrl[['gene', 'Symbol', 'GeneID', 'cluster', 'DiffExp', 'avg_logFC', 'logFC']]
merged_gene_list_0v2_ctrl.columns = ["gene_hou", "gene_hartig", "GeneID", "cluster", "DiffExp", "logFC_hou", "logFC_hartig"]


# In[13]:


# create a new data frame that contains only those rows that have matching values in both of the original data frames
# based on the 'gene' column in Hou and 'symbol' column in Hartig
merged_gene_list_0v2_treat = pd.merge(left=hou, right=hartig_0v2_treat, left_on='gene', right_on='Symbol')

# re-order & re-name columns as necessary
merged_gene_list_0v2_treat = merged_gene_list_0v2_treat[['gene', 'Symbol', 'GeneID', 'cluster', 'DiffExp', 'avg_logFC', 'logFC']]
merged_gene_list_0v2_treat.columns = ["gene_hou", "gene_hartig", "GeneID", "cluster", "DiffExp", "logFC_hou", "logFC_hartig"]


# In[14]:


# create a new data frame that contains only those rows that have matching values in both of the original data frames
# based on the 'gene' column in Hou and 'symbol' column in Hartig
merged_gene_list_0v2_inter = pd.merge(left=hou, right=hartig_0v2_inter, left_on='gene', right_on='Symbol')

# re-order & re-name columns as necessary
merged_gene_list_0v2_inter = merged_gene_list_0v2_inter[['gene', 'Symbol', 'GeneID', 'cluster', 'DiffExp', 'avg_logFC', 'logFC']]
merged_gene_list_0v2_inter.columns = ["gene_hou", "gene_hartig", "GeneID", "cluster", "DiffExp", "logFC_hou", "logFC_hartig"]


# In[52]:


sns.set(style="whitegrid")
sns.set(font_scale = 1.5)
#plt.figure(figsize=(8,4))

ax = sns.boxplot(x="cluster", y="logFC_hartig", data=merged_gene_list_0dpa, showfliers = False, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic']).set_title("0 dpa")

ax = sns.stripplot(x="cluster", y="logFC_hartig", hue="DiffExp", palette=["red","grey","blue"], data=merged_gene_list_0dpa, color="0.5", alpha=0.25, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'])

ax.set_xticks(range(6))
ax.set_xticklabels(['Basal \nEpithelial','Intermediate \nEpithelial','Superficial \nEpithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], rotation=80)
ax.set(xlabel=None)
ax.axhline(y=0, linewidth=4, color='black')
ax.legend([],[], frameon=False)
#ax.legend(bbox_to_anchor=(1.02, 0.55), loc='upper left', borderaxespad=0)

ax.set(ylim=(-4, 4))

plt.savefig("boxplot_0dpa", bbox_inches='tight', transparent=True)
plt.show()


# In[53]:


#plt.figure(figsize=(8,4))

ax = sns.boxplot(x="cluster", y="logFC_hartig", data=merged_gene_list_2dpa, showfliers = False, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic']).set_title("2 dpa")
ax = sns.stripplot(x="cluster", y="logFC_hartig", hue="DiffExp", palette=["grey","blue"], data=merged_gene_list_2dpa, color="0.5", alpha=0.25, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'])

ax.set_xticks(range(6))
ax.set_xticklabels(['Basal \nEpithelial','Intermediate \nEpithelial','Superficial \nEpithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], rotation=80)
ax.set(xlabel=None)
ax.axhline(y=0, linewidth=4, color='black')
#ax.legend(bbox_to_anchor=(1.02, 0.55), loc='upper left', borderaxespad=0)
ax.legend([],[], frameon=False)

ax.set(ylim=(-4, 4))

plt.savefig("boxplot_2dpa", bbox_inches='tight', transparent=True)
plt.show()


# In[54]:


#plt.figure(figsize=(8,4))

ax = sns.boxplot(x="cluster", y="logFC_hartig", data=merged_gene_list_4dpa, showfliers = False, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic']).set_title("4 dpa")
ax = sns.stripplot(x="cluster", y="logFC_hartig", hue="DiffExp", palette=["grey","red"], data=merged_gene_list_4dpa, color="0.5", alpha=0.25, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'])

ax.set_xticks(range(6))
ax.set_xticklabels(['Basal \nEpithelial','Intermediate \nEpithelial','Superficial \nEpithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], rotation=80)
ax.set(xlabel=None)
ax.axhline(y=0, linewidth=4, color='black')
#ax.legend(bbox_to_anchor=(1.02, 0.55), loc='upper left', borderaxespad=0)
ax.legend([],[], frameon=False)

ax.set(ylim=(-4, 4))

plt.savefig("boxplot_4dpa", bbox_inches='tight', transparent=True)
plt.show()


# In[55]:


#plt.figure(figsize=(8,4))

ax = sns.boxplot(x="cluster", y="logFC_hartig", data=merged_gene_list_0v2_ctrl, showfliers = False, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic']).set_title("Control \n(vehicle)")
ax = sns.stripplot(x="cluster", y="logFC_hartig", hue="DiffExp", palette=["grey","blue","red"], data=merged_gene_list_0v2_ctrl, color="0.5", alpha=0.25, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'])

ax.set_xticks(range(6))
ax.set_xticklabels(['Basal \nEpithelial','Intermediate \nEpithelial','Superficial \nEpithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], rotation=80)
ax.set(xlabel=None)
ax.axhline(y=0, linewidth=4, color='black')
#ax.legend(bbox_to_anchor=(1.02, 0.55), loc='upper left', borderaxespad=0)
ax.legend([],[], frameon=False)

ax.set(ylim=(-5, 6))

plt.savefig("boxplot_0v2_ctrl", bbox_inches='tight', transparent=True)
plt.show()


# In[56]:


#plt.figure(figsize=(8,4))

treat = sns.boxplot(x="cluster", y="logFC_hartig", data=merged_gene_list_0v2_treat, showfliers = False, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic']).set_title("Treated \n(cortisol)")
treat = sns.stripplot(x="cluster", y="logFC_hartig", hue="DiffExp", palette=["grey","blue","red"], data=merged_gene_list_0v2_treat, color="0.5", alpha=0.25, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'])

treat.set_xticks(range(6))
treat.set_xticklabels(['Basal \nEpithelial','Intermediate \nEpithelial','Superficial \nEpithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], rotation=80)
treat.set(xlabel=None)
treat.axhline(y=0, linewidth=4, color='black')
#treat.legend(bbox_to_anchor=(1.02, 0.55), loc='upper left', borderaxespad=0)
treat.legend([],[], frameon=False)

treat.set(ylim=(-5, 6))

plt.savefig("boxplot_0v2_treat", bbox_inches='tight', transparent=True)
plt.show()


# In[57]:


#plt.figure(figsize=(8,4))

inter= sns.boxplot(x="cluster", y="logFC_hartig", data=merged_gene_list_0v2_inter, showfliers = False, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic']).set_title("Interaction \n(difference)")
inter= sns.stripplot(x="cluster", y="logFC_hartig", hue="DiffExp", palette=["grey","blue","red"], data=merged_gene_list_0v2_inter, color="0.5", alpha=0.25, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'])

inter.set_xticks(range(6))
inter.set_xticklabels(['Basal \nEpithelial','Intermediate \nEpithelial','Superficial \nEpithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], rotation=80)
inter.set(xlabel=None)
inter.axhline(y=0, linewidth=4, color='black')
#inter.legend(bbox_to_anchor=(1.02, 0.55), loc='upper left', borderaxespad=0)
inter.legend([],[], frameon=False)

inter.set(ylim=(-5, 6))

plt.savefig("boxplot_0v2_inter", bbox_inches='tight', transparent=True)
plt.show()


# In[ ]:


fig, axes = plt.subplots(nrows=3, ncols=1, dpi = 100, constrained_layout=True, figsize=(8,8))
sns.set(style="whitegrid")

sns.boxplot(x="cluster", y="logFC_hartig", data=merged_gene_list_0v2_inter, showfliers = False, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], ax=axes[0]).set_title("0dpa v. 2dpa \ninteraction")
sns.stripplot(x="cluster", y="logFC_hartig", hue="DiffExp", palette=["silver","blue","red"], data=merged_gene_list_0v2_inter, color="0.5", alpha=0.25, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], ax=axes[0])
axes[0].legend(bbox_to_anchor=(1.02, 0.55), loc='upper left', borderaxespad=0)
axes[0].set(xticklabels=[])
axes[0].set(xlabel=None)

sns.boxplot(x="cluster", y="logFC_hartig", data=merged_gene_list_0v2_ctrl, showfliers = False, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], ax=axes[1]).set_title("Veh")
sns.stripplot(x="cluster", y="logFC_hartig", hue="DiffExp", palette=["silver","blue","red"], data=merged_gene_list_0v2_ctrl, color="0.5", alpha=0.25, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], ax=axes[1])
axes[1].legend(bbox_to_anchor=(1.02, 0.55), loc='upper left', borderaxespad=0)
axes[1].set(xticklabels=[])
axes[1].set(xlabel=None)

sns.boxplot(x="cluster", y="logFC_hartig", data=merged_gene_list_0v2_treat, showfliers = False, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], ax=axes[2]).set_title("Cort")
sns.stripplot(x="cluster", y="logFC_hartig", hue="DiffExp", palette=["silver","blue","red"], data=merged_gene_list_0v2_treat, color="0.5", alpha=0.25, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], ax=axes[2])
axes[2].legend(bbox_to_anchor=(1.02, 0.55), loc='upper left', borderaxespad=0)
axes[2].set_xticklabels(axes[2].get_xticklabels(),rotation = 90)

plt.savefig("0v2_boxplots", bbox_inches='tight')


# In[ ]:


fig, axes = plt.subplots(nrows=3, ncols=1, dpi = 100, constrained_layout=True, figsize=(8,8))
sns.set(style="whitegrid")

sns.boxplot(x="cluster", y="logFC_hartig", data=merged_gene_list_0dpa, showfliers = False, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], ax=axes[0]).set_title("Veh vs. Cort \n0dpa")
sns.stripplot(x="cluster", y="logFC_hartig", hue="DiffExp", palette=["red","silver","blue"], data=merged_gene_list_0dpa, color="0.5", alpha=0.25, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], ax=axes[0])
axes[0].legend(bbox_to_anchor=(1.02, 0.55), loc='upper left', borderaxespad=0)
axes[0].set(xticklabels=[])
axes[0].set(xlabel=None)

sns.boxplot(x="cluster", y="logFC_hartig", data=merged_gene_list_2dpa, showfliers = False, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], ax=axes[1]).set_title("2dpa")
sns.stripplot(x="cluster", y="logFC_hartig", hue="DiffExp", palette=["silver","blue"], data=merged_gene_list_2dpa, color="0.5", alpha=0.25, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], ax=axes[1])
axes[1].legend(bbox_to_anchor=(1.02, 0.55), loc='upper left', borderaxespad=0)
axes[1].set(xticklabels=[])
axes[1].set(xlabel=None)

sns.boxplot(x="cluster", y="logFC_hartig", data=merged_gene_list_4dpa, showfliers = False, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], ax=axes[2]).set_title("4dpa")
sns.stripplot(x="cluster", y="logFC_hartig", hue="DiffExp", palette=["silver","red"], data=merged_gene_list_4dpa, color="0.5", alpha=0.25, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], ax=axes[2])
axes[2].legend(bbox_to_anchor=(1.02, 0.55), loc='upper left', borderaxespad=0)
axes[2].set_xticklabels(axes[2].get_xticklabels(),rotation = 90)

plt.savefig("0_2_4_boxplots", bbox_inches='tight')


# In[ ]:


fig, axes = plt.subplots(nrows=3, ncols=1, dpi = 100, constrained_layout=True, figsize=(8,8))
sns.set(style="whitegrid")

sns.violinplot(x="cluster", y="logFC_hartig", inner=None, color=".8", data=merged_gene_list_0dpa, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], ax=axes[0]).set_title("Veh vs. Cort \n0dpa")
sns.stripplot(x="cluster", y="logFC_hartig", hue="DiffExp", palette=["red","grey","blue"], data=merged_gene_list_0dpa, color="0.5", alpha=0.25, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], ax=axes[0])
axes[0].legend(bbox_to_anchor=(1.02, 0.55), loc='upper left', borderaxespad=0)
axes[0].set(xticklabels=[])
axes[0].set(xlabel=None)

sns.violinplot(x="cluster", y="logFC_hartig", inner=None, color=".8", data=merged_gene_list_2dpa, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], ax=axes[1]).set_title("2dpa")
sns.stripplot(x="cluster", y="logFC_hartig", hue="DiffExp", palette=["grey","blue"], data=merged_gene_list_2dpa, color="0.5", alpha=0.25, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], ax=axes[1])
axes[1].legend(bbox_to_anchor=(1.02, 0.55), loc='upper left', borderaxespad=0)
axes[1].set(xticklabels=[])
axes[1].set(xlabel=None)

sns.violinplot(x="cluster", y="logFC_hartig", inner=None, color=".8", data=merged_gene_list_4dpa, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], ax=axes[2]).set_title("4dpa")
sns.stripplot(x="cluster", y="logFC_hartig", hue="DiffExp", palette=["grey","red"], data=merged_gene_list_4dpa, color="0.5", alpha=0.25, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], ax=axes[2])
axes[2].legend(bbox_to_anchor=(1.02, 0.55), loc='upper left', borderaxespad=0)
axes[2].set_xticklabels(axes[2].get_xticklabels(),rotation = 90)

plt.savefig("0_2_4_violinplots", bbox_inches='tight')


# In[ ]:


fig, axes = plt.subplots(nrows=3, ncols=1, dpi = 100, constrained_layout=True, figsize=(8,8))
sns.set(style="whitegrid")

sns.violinplot(x="cluster", y="logFC_hartig", inner=None, color=".8", data=merged_gene_list_0v2_inter, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], ax=axes[0]).set_title("0dpa v. 2dpa \ninteraction")
sns.stripplot(x="cluster", y="logFC_hartig", hue="DiffExp", palette=["grey","blue"], data=merged_gene_list_0v2_inter, color="0.5", alpha=0.25, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], ax=axes[0])
axes[0].legend(bbox_to_anchor=(1.02, 0.55), loc='upper left', borderaxespad=0)
axes[0].set(xticklabels=[])
axes[0].set(xlabel=None)

sns.violinplot(x="cluster", y="logFC_hartig", inner=None, color=".8", data=merged_gene_list_0v2_ctrl, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], ax=axes[1]).set_title("Veh")
sns.stripplot(x="cluster", y="logFC_hartig", hue="DiffExp", palette=["grey","blue","red"], data=merged_gene_list_0v2_ctrl, color="0.5", alpha=0.25, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], ax=axes[1])
axes[1].legend(bbox_to_anchor=(1.02, 0.55), loc='upper left', borderaxespad=0)
axes[1].set(xticklabels=[])
axes[1].set(xlabel=None)

sns.violinplot(x="cluster", y="logFC_hartig", inner=None, color=".8", data=merged_gene_list_0v2_treat, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], ax=axes[2]).set_title("Cort")
sns.stripplot(x="cluster", y="logFC_hartig", hue="DiffExp", palette=["grey","blue","red"], data=merged_gene_list_0v2_treat, color="0.5", alpha=0.25, order=['Basal Epithelial','Intermediate Epithelial','Superficial Epithelial', 'Mucosal-like', 'Mesenchymal', 'Hematopoietic'], ax=axes[2])
axes[2].legend(bbox_to_anchor=(1.02, 0.55), loc='upper left', borderaxespad=0)
axes[2].set_xticklabels(axes[2].get_xticklabels(),rotation = 90)

plt.savefig("0v2_violinplots", bbox_inches='tight')


# In[ ]:




