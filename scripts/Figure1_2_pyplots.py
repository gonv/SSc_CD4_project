#!/usr/bin/env python
# coding: utf-8

# # Case study: CD4+ cells from systemic sclerosis patients and healthy controls
# 
# The greatest and noblest pleasure which men can have in this world is to discover new truths; and the next is to shake off old prejudices.
# 
#                            Frederick II of Prussia
# 
# At the heart of science is an essential balance between two seemingly contradictory attitudes—an openness to new ideas, no matter how bizarre or counterintuitive, and the most ruthlessly skeptical scrutiny of all ideas, old and new. This is how deep truths are winnowed from deep nonsense.
# 
#                                   Carl Sagan

# ## 1. Loading deck

# In[2]:


import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
from matplotlib import rc_context
from matplotlib import pyplot as plt, cm as mpl_cm
from cycler import cycler
import matplotlib_inline.backend_inline

import seaborn as sb
#from gprofiler import gprofiler
#import rpy2.rinterface_lib.callbacks
#import logging
#import os

import time
from datetime import timedelta
t1 = time.time() 

import scanpy.external as sce
#import scvi
from sklearn.metrics import silhouette_score, silhouette_samples
#import celltypist
#from celltypist import models
import os #to set wd
#import scirpy as ir


# In[3]:


# Ignore R warning messages
#Note: this can be commented out to get more verbose R output

plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=75, dpi_save=300, facecolor='white')
sc.logging.print_header()


# In[4]:


print("Current working directory: {0}".format(os.getcwd()))
os.chdir('/lab-share/IM-Gutierrez-e2/Public/gvillanueva/00_CD4/sc_ob')
print("Current working directory: {0}".format(os.getcwd()))
#


# In[5]:


CD4_pallete = {
    'Tn': '#a9dfbb',
    'Th1_GZMK': '#9B59B6',
    'Th1_GZMB': '#2ECC71',
    'Th2': '#3faccd',
    'Th17':'#FF5733',
    'Treg': '#fc6daf',
    'Thf':'#F1C40F',
         'Tn': '#a9dfbb',
         'Tcm': '#6cd0c4',
         'Treg': '#fc6daf',
         'Tem': '#3faccd',
         'Temra': '#ff0000',
    '0':"#CDB4DB",'1': "#F78181","2": "#82CFFD",
    "6": "#FFB347","4":"#AEB6BF","5" :"#F7ACB4",
    "7": "#7FFFD4","10":"#FFA07A","8":"#F4D03F","9":'#F8D5E1',
    "3":"#8FCFD1","11":"#2ECC71"
    
}


# In[6]:


adata3 = sc.read('scRNAseq_SSc_CD4_allMerge_20PCs_TCR_it92_2.h5ad')
#


# In[6]:


adata3


# In[7]:


adata3.uns["neighbors"]["params"]


# ## 2. Density plot

# In[8]:


sc.tl.embedding_density(adata3, basis='umap', key_added="condition_dens",groupby='condition')
#


# In[9]:


sc.pl.embedding_density(adata3, 
                        key= "condition_dens", 
                        wspace=0.3,
                        save='scRNAseq_CD4_SSc_UMAP_DensityMap_byCON_it92.png')
#


# ## 3. DEA

# ### 3.1 Cluster markers

# In[10]:


sc.tl.rank_genes_groups(adata3, 'R09', method='wilcoxon', key_added='ClusterMarkers')


# In[11]:


resultCM = adata3.uns['ClusterMarkers']
groupsCM = resultCM['names'].dtype.names
ClusterMarkers = pd.DataFrame({group + '_' + key[:15]: resultCM[key][group]
                        for group in groupsCM for key in ['names','scores','pvals_adj','logfoldchanges']})
ClusterMarkers.head(10)


# ### 3.2 Other markers

# In[10]:


MarkersCluster = ["FOS", "PRDX2", "PNISR", "UCP2", "TPT1", "IL32", "TIGIT", "LEF1", "LRPAP1", "HLA-DRB5", "FCRL3", "GZMH",
         "IFITM1", "PNRC1", "SRP14", "EVL", "SEPTIN1", "LSM7", "LIMS1", "RIPOR2", "DCTN2", "HLA-DRA", "TMED2", "PLEK",
         "JUNB", "TRIM22", "ARHGDIB", "SSBP1", "SSR2", "UQCC2", "CMTM7", "C12orf57", "SYTL1", "H3F3A", "ARHGDIA", "GNLY",
         "LTB", "TSHZ2", "ADTRP", "NDUFB9", "CDK5RAP3", "DARS", "LINC00892", "DGKA", "XRCC6", "GZMA", "ARL6IP1", "PEPD",
         "DUSP1", "CYCS", "LSP1", "TBC1D10C", "KCNAB2", "UBE2V1", "CD2", "NSA2", "GAS5", "HMGB1", "CX3CR1"]


# In[11]:


# Perform the dotplot with the modified order of cell types and custom color scale
dp = sc.pl.dotplot(adata3, MarkersCluster, 'R09',standard_scale='var' ,return_fig=True, dendrogram=False)
dp.style(cmap= "Blues" ,dot_edge_color='white', dot_edge_lw=0.5).show()


# In[65]:


MarkersCluster2 = pd.Series(["CCR7", "EEF1B2", "TUBA1A", "CD7", "NACA", "MALAT1", "LEF1",
                             "SELL", "FOSB", "FOS", "IER2","JUN","KLRB1", "LTB", "AQP3",  
                             "SARAF", "ZNF331", "EIF1", "TSC22D3", "VIM","ANXA1", "ITGB1", 
                             "FTH1", "JUNB", "SRGN", "NR4A2", "GZMK", "LYAR", "TPT1", 
                             "HIST1H1E", "HOOK2", "DDX5", "EEF2","LMNA", "CRIP1", "ZFP36", 
                             "S100A10",  "FOXP3", "RTKN2", "TIGIT", "RGS1", "IL10RA", "S100A4",
                             "SH3BGRL3", "NKG7", "CCL5","GZMA", "GZMH", "CST7", "GNLY"])
MarkersCluster3 = MarkersCluster2.unique()


# In[66]:


MarkersCluster3


# In[67]:


MarkersCluster3 is MarkersCluster2


# In[53]:


# Define the custom order of R09 categories
custom_order = [0, 1, 2, 4, 8, 3, 5, 6, 9, 7, 11, 10]


# In[47]:


# Perform the dotplot with the modified order of cell types and custom color scale
dp = sc.pl.dotplot(adata3, MarkersCluster3, 'R09',standard_scale='var' ,return_fig=True, dendrogram=False)
dp.style(cmap= "Blues" ,dot_edge_color='white', dot_edge_lw=0.5).show()


# In[68]:


# Perform the dotplot with the modified order of cell types and custom color scale
dp = sc.pl.dotplot(adata3, MarkersCluster3, 'R09',swap_axes=True,standard_scale='var',
                   return_fig=True, dendrogram=True, save='fiugres/CD4_clusterMarkers_v2.pdf')
dp.style(cmap= "viridis" ,dot_edge_color='white', dot_edge_lw=0.5).show()


# In[69]:


# Generar el dotplot
dp = sc.pl.dotplot(
    adata3, 
    MarkersCluster2, 
    groupby='R09', 
    swap_axes=True, 
    standard_scale='var',
    return_fig=True, 
    dendrogram=False
)

# Aplicar estilo y mostrar el plot
dp.style(cmap="viridis", dot_edge_color='white', dot_edge_lw=0.5)

# Guardar la figura
dp.savefig('figures/CD4_clusterMarkers_v2.pdf')

# Mostrar la figura
dp.show()


# In[80]:


# Generar el dotplot
dp = sc.pl.dotplot(
    adata3, 
    MarkersCluster2, 
    groupby='R09', 
    swap_axes=True, 
    standard_scale='var',
    return_fig=True, 
    dendrogram=False
)

# Aplicar estilo y mostrar el plot
dp.style(cmap="viridis", dot_edge_color='white', dot_edge_lw=0.5)

# Guardar la figura
dp.savefig('figures/CD4_clusterMarkers_v2.pdf')

# Mostrar la figura
dp.show()


# In[84]:


# Generar el dotplot
dp = sc.pl.dotplot(
    adata3, 
    MarkersCluster2, 
    groupby='R09', 
    swap_axes=True, 
    standard_scale='var',
    return_fig=True, 
    dendrogram=False
)

# Aplicar estilo y mostrar el plot
dp.style(cmap="viridis", dot_edge_color='white', dot_edge_lw=0.5)

axes_dict = dp.get_axes()


# ## 4. Ploting UMAPs

# In[20]:


cluster2annotation4 = {
    '0': '0_CCR7_EEF1B2_TUBA1A',
    '1': '1_MALAT1_EEF1B2_CCR7',
    '2': '2_FOSB_FOS_IER2',
    '3': '3_S100A4_KLRB1_LTB',
    '4': '4_YPEL5_EIF1_SARAF',
    '5': '5_VIM_ANXA1_S100A4',
    '6': '6_FTH1_JUNB_SRGN',
    '7': '7_SRGN_GZMK_ZFP36',
    '8': '8_TPT1_DDX5_EEF2',
    '9':'9_S100A4_LMNA_S100A10',
    '10':'10_FOXP3_RTKN2_TIGIT',
    '11':'11_NKG7_CCL5_GZMH',
}
adata3.obs['R09_2'] = adata3.obs['R09'].map(cluster2annotation4).astype('category')


# In[21]:


CD4_pallete = {
    'Tn': '#a9dfbb',
    'Th1_GZMK': '#9B59B6',
    'Th1_GZMB': '#2ECC71',
    'Th2': '#3faccd',
    'Th17':'#FF5733',
    'Treg': '#fc6daf',
    'Thf':'#F1C40F',
         'Tn': '#a9dfbb',
         'Tcm': '#6cd0c4',
         'Treg': '#fc6daf',
         'Tem': '#3faccd',
         'Temra': '#ff0000',
    '0':"#CDB4DB",'1': "#F78181","2": "#82CFFD",
    "6": "#FFB347","4":"#AEB6BF","5" :"#F7ACB4",
    "7": "#7FFFD4","10":"#FFA07A","8":"#F4D03F","9":'#F8D5E1',
    "3":"#8FCFD1","11":"#2ECC71",
    
        '0_CCR7_EEF1B2_TUBA1A':"#CDB4DB",'1_MALAT1_EEF1B2_CCR7': "#F78181","2_FOSB_FOS_IER2": "#82CFFD",
    "6_FTH1_JUNB_SRGN": "#FFB347","4_YPEL5_EIF1_SARAF":"#AEB6BF","5_VIM_ANXA1_S100A4" :"#F7ACB4",
    "7_SRGN_GZMK_ZFP36": "#7FFFD4","10_FOXP3_RTKN2_TIGIT":"#FFA07A","8_TPT1_DDX5_EEF2":"#F4D03F","9_S100A4_LMNA_S100A10":'#F8D5E1',
    "3_S100A4_KLRB1_LTB":"#8FCFD1","11_NKG7_CCL5_GZMH":"#2ECC71"
}


# In[22]:


with rc_context({'figure.figsize': (5, 6)}):
    sc.pl.umap(adata3, color = ['R09_2','R09'],ncols=2 ,add_outline=True,
               size=10,legend_fontsize=12, legend_fontoutline=1,frameon=True, wspace=0.55,palette= CD4_pallete)


# In[23]:


with rc_context({'figure.figsize': (5, 6)}):
    sc.pl.umap(adata3, color = ['R09_2','R09'],ncols=2 ,add_outline=True,
               size=10,legend_fontsize=12, legend_fontoutline=1,frameon=True, wspace=0.55,palette= CD4_pallete)


# In[ ]:


sc.pl.embedding_density(adata3, key= "condition_dens", wspace=0.3,save='scRNAseq_CD4_SSc_UMAP_DensityMap_byCON_it92.pdf')
#


# In[ ]:


sc.pl.embedding_density(adatanoCTRL, key= "clinical_dens", wspace=0.3,save='scRNAseq_CD4_SSc_UMAP_DensityMap_bySub_it92.pdf')
#


# In[24]:


# Define matplotlib Axes
# Number of Axes & plot size
ncols = 2
nrows = 1
figsize = 4
wspace = 0.5
fig, axs = plt.subplots(
    nrows=nrows,
    ncols=ncols,
    figsize=(ncols * figsize + figsize * wspace * (ncols - 1), nrows * figsize),
)
plt.subplots_adjust(wspace=wspace)
# This produces two Axes objects in a single Figure
print("axes:", axs)

# We can use these Axes objects individually to plot on them
# We need to set show=False so that the Figure is not displayed before we
# finished plotting on all Axes and making all plot adjustments
sc.pl.umap(adata3, color = ['R09'],ncols=2 ,add_outline=True,
               size=10,legend_fontsize=12, legend_fontoutline=1,frameon=True, wspace=0.55,palette= CD4_pallete, ax=axs[0], show=False)
# Example zoom-in into a subset of louvain clusters
sc.pl.umap(adata3, color = ['cell_type'],ncols=2 ,add_outline=True,
               size=10,legend_fontsize=12, legend_fontoutline=1,frameon=True, wspace=0.55,palette= CD4_pallete, ax=axs[1])


# In[ ]:


import math

# Get the angles from 0 to 2 pie (360 degree) in narray object
X = np.arange(0, math.pi*2, 0.05)

# Using built-in trigonometric function we can directly plot
# the given cosine wave for the given angles
Y1 = np.sin(X)
Y2 = np.cos(X)
Y3 = np.tan(X)
Y4 = np.tanh(X)

# Initialise the subplot function using number of rows and columns
figure, axis = plt.subplots(2, 2)

# For Sine Function
axis[0, 1].plot(X, Y1)
axis[0, 1].set_title("Sine Function")

# For Cosine Function
axis[0, 0].plot(X, Y2)
axis[0, 0].set_title("Cosine Function")

# For Tangent Function
axis[1, 0].plot(X, Y3)
axis[1, 0].set_title("Tangent Function")

# For Tanh Function
axis[1, 1].plot(X, Y4)
axis[1, 1].set_title("Tanh Function")

# Combine all the operations and display
plt.show()


# ## 5. Subclustering

# In[7]:


adataTh = sc.read('../sc_ob/scRNAseq_SSc_CD4_ThSubclustering.h5ad')
#


# In[8]:


ThMarkers =["NCAM1","IL22","NOTCH2","RORA","RORC", "TBX21", "GATA3", "STAT4", "STAT3", "STAT1","STAT5A","STAT5B","STAT6", "IRF4", "HLA-DRA", "CCR1","CCR2","CCR3","CCR4","CCR5","CCR6","CCR8", "CCR10","CXCR5","ICOS","CXCR3","FOXO4","AHR","FOXP3","FCRL3"]
#


# In[9]:


ThMarkers =["NCAM1","IL22","NOTCH2","RORA","RORC", "TBX21", "GATA3", "STAT4", "STAT3", "STAT1","STAT5A","STAT5B","STAT6",
            "IRF4", "HLA-DRA", "CCR1","CCR2","CCR3","CCR4","CCR5","CCR6","CCR8", "CCR10","CXCR5","ICOS","CXCR3","FOXO4",
            "AHR","FOXP3","FCRL3"]
#


# In[10]:


ThMarkers =["GATA3","STAT4","STAT3","STAT1","STAT5A","STAT5B","STAT6","RORC","TBX21","CCR4","CCR6","CXCR5","ICOS","CXCR3","FOXO4","AHR"]
#


# In[11]:


#Visualization (UMAP) && Clustering {Leiden&Louvain}
with rc_context({'figure.figsize': (4, 5)}):
    sc.pl.umap(adataTh, color = ['R09','Th'], add_outline=True, legend_loc='on data',
               legend_fontsize=12, ncols=1,legend_fontoutline=1,frameon=True, wspace=0.3,size=10,title=['Subclustering name by cluster','Subclustering name by annotation'],
               save='scRNAseq_CD4_SSc_UMAP_Subclustering_it92.pdf')


# In[12]:


dp = sc.pl.dotplot(adataTh, ThMarkers, 'Th',standard_scale='var',return_fig=True, dendrogram=False,
                   save='scRNAseq_CD4_SSc_Dotplot_Subclustering_it92.pdf')
dp.style(cmap= "Blues" ,dot_edge_color='white', dot_edge_lw=0.5).show()


# In[14]:


# Genera el dotplot con los argumentos deseados
dotplot = sc.pl.dotplot(adataTh, ThMarkers, groupby='Th', standard_scale='var', return_fig=True, dendrogram=False, swap_axes=True)

# Modifica el estilo del dotplot
dotplot.style(cmap='Blues', dot_edge_color='white', dot_edge_lw=0.5)

# Guarda la figura si es necesario
#dotplot.savefig('scRNAseq_CD4_SSc_Dotplot_Subclustering_v3_it92.png')
dotplot.savefig('../analysis/scRNAseq_CD4_SSc_Dotplot_Subclustering_v3_it92.png')

# Muestra la figura
dotplot.show()


# ## 6. Annotation

# In[31]:


Soskic_makers_list = ['CCR7','FHIT','CREM','SELL','KLRB1','ITGB1','S100A4','S100A6','LYAR','GZMA','CCL5','NKG7','GNLY','GZMH','GZMB','TIGIT','FOXP3','IKZF2','IL2RA']


# In[32]:


dp = sc.pl.dotplot(adata3, Soskic_makers_list, 'R09',standard_scale='var',return_fig=True, dendrogram=False,swap_axes=True)
dp.style(cmap= "Blues" ,dot_edge_color='white', dot_edge_lw=0.5).show()
#dotplot.savefig('scRNAseq_CD4_SSc_Dotplot_Annotation_it92.png')


# In[35]:


sc.settings.set_figure_params(figsize=(4, 6),dpi=100)

dp = sc.pl.dotplot(
    adata3,
    Soskic_makers_list,
    groupby="R09",
    standard_scale="var",
    return_fig=True,
    dendrogram=False,
    swap_axes=True
)

dp.style(
    cmap="Blues",
    dot_edge_color="white",
    dot_edge_lw=0.5,
    dot_max=0.6,
    dot_min=0.1
)

dp.legend(show=True,)

dp.show()

dp.savefig('../sc_ob/figures/scRNAseq_CD4_SSc_Dotplot_Annotation_v2_it92.png')


# In[36]:


dp.show()


# In[52]:


fig, ax = plt.subplots(figsize=(3.25,7))
dp = sc.pl.dotplot(adata3, 
                   Soskic_makers_list, 
                   'R09',
                   standard_scale='var',
                   return_fig=True, 
                   ax=ax, 
                   dendrogram=False,
                   swap_axes=True, 
                   show=False)
dp.style(cmap= "Blues" ,dot_edge_color='white', dot_edge_lw=0.5)
dp.legend(show=False,)
ax.tick_params(axis="x", labelrotation=90)
dp.show()


# In[ ]:




