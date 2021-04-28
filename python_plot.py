#!/usr/bin/env python
#
# python_plot.py
# author: H. Kim
# date created: 2017, Aug.
# date last modified: 2017, Aug.
#
# input:
#
#   type_plot:
#   fname_in:
#   fname_out:
#
# options:
#
#   --x: x column name in df
#   --y: y column name in df
#   --h: hue column name in df
#   --cmap: color map
#     Possible values are: Accent, Accent_r, Blues, Blues_r, BrBG, BrBG_r, BuGn, BuGn_r, BuPu, BuPu_r, CMRmap, CMRmap_r, Dark2, Dark2_r, GnBu, GnBu_r, Greens, Greens_r, Greys, Greys_r, OrRd, OrRd_r, Oranges, Oranges_r, PRGn, PRGn_r, Paired, Paired_r, Pastel1, Pastel1_r, Pastel2, Pastel2_r, PiYG, PiYG_r, PuBu, PuBuGn, PuBuGn_r, PuBu_r, PuOr, PuOr_r, PuRd, PuRd_r, Purples, Purples_r, RdBu, RdBu_r, RdGy, RdGy_r, RdPu, RdPu_r, RdYlBu, RdYlBu_r, RdYlGn, RdYlGn_r, Reds, Reds_r, Set1, Set1_r, Set2, Set2_r, Set3, Set3_r, Spectral, Spectral_r, Vega10, Vega10_r, Vega20, Vega20_r, Vega20b, Vega20b_r, Vega20c, Vega20c_r, Wistia, Wistia_r, YlGn, YlGnBu, YlGnBu_r, YlGn_r, YlOrBr, YlOrBr_r, YlOrRd, YlOrRd_r, afmhot, afmhot_r, autumn, autumn_r, binary, binary_r, bone, bone_r, brg, brg_r, bwr, bwr_r, cool, cool_r, coolwarm, coolwarm_r, copper, copper_r, cubehelix, cubehelix_r, flag, flag_r, gist_earth, gist_earth_r, gist_gray, gist_gray_r, gist_heat, gist_heat_r, gist_ncar, gist_ncar_r, gist_rainbow, gist_rainbow_r, gist_stern, gist_stern_r, gist_yarg, gist_yarg_r, gnuplot, gnuplot2, gnuplot2_r, gnuplot_r, gray, gray_r, hot, hot_r, hsv, hsv_r, inferno, inferno_r, jet, jet_r, magma, magma_r, nipy_spectral, nipy_spectral_r, ocean, ocean_r, pink, pink_r, plasma, plasma_r, prism, prism_r, rainbow, rainbow_r, seismic, seismic_r, spectral, spectral_r, spring, spring_r, summer, summer_r, terrain, terrain_r, viridis, viridis_r, winter, winter_r
#   --xlabel: xlabel
#   --ylabel: ylabel
#   --title: title
#   --fontsize_xlabel: (default=18)
#   --fontsize_ylabel: (default=18)
#   --fontsize_title: (default=20)
#   --figsize_x: figsize x (default=5)
#   --figsize_y: figsize y (default=5)
#   --fname_map_fname2barcode: mapping file name for maping vcf file fnames to barcodes
#   --type_output: {'all','short'}
#     
# output:
#
#   a pdf file
#
# test:
#
# git clone https://github.com/hyunsoo77/demo.git
# cd demo
# python_plot.py --xlabel xx --ylabel yy --title title --figsize_y 5 heatmap input/table_mtx.txt output/table_mtx_heatmap.pdf
# acroread table_mtx_clustermap.pdf
#
# python_plot.py --xlabel xx --ylabel yy --title title --fname_row_colors input/table_row_colors.txt --fname_col_colors input/table_col_colors.txt clustermap input/table_mtx.txt output/table_mtx_clustermap.pdf
# acroread table_mtx_clustermap.pdf
#
# python_plot.py --xlabel xx --ylabel yy --title title --fname_row_colors input/table_row2_colors.txt --fname_col_colors input/table_col2_colors.txt clustermap input/table_mtx.txt output/table_mtx_clustermap.pdf
# acroread table_mtx_clustermap.pdf
#
#
import os
import sys
import psutil
import re
import ssl
import math
import time
import datetime
import glob
#import pickle
import ujson
import shelve
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from natsort import natsorted, ns
from optparse import OptionParser
from collections import namedtuple
from collections import defaultdict
from itertools import compress
# comment out the following line
#from print_log import print_log

# main routine
usage = "usage: %prog [options] type_plot fname_in fname_out"
parser = OptionParser(usage)
parser.add_option("-l", "--n_log_level", type="int", default=0, help='n_log_level')
parser.add_option("--n_log", type="int", default=1, help='n_log')
parser.add_option("--n_debug", type="int", default=0, help='n_debug')
parser.add_option("--x", default="", help="x column name in df")
parser.add_option("--y", default="", help="y column name in df")
parser.add_option("--hue", default="", help="hue column name in df")
parser.add_option("--colseq", default="", help="colsep column name in df")
parser.add_option("--cmap", default="BuGn", help="color map")
parser.add_option("--xlabel", default="", help="xlabel")
parser.add_option("--ylabel", default="", help="ylabel")
parser.add_option("--ylabel_cbar", default="", help="ylabel in colorbar")
parser.add_option("--title", default="", help="title")
parser.add_option("--fontsize_xlabel", type="int", default=18, help="xlabel")
parser.add_option("--fontsize_ylabel", type="int", default=18, help="ylabel")
parser.add_option("--fontsize_title", type="int", default=20, help="title")
parser.add_option("--figsize_x", type="float", default=4, help='figsize x')
parser.add_option("--figsize_y", type="float", default=4, help='figsize y')
parser.add_option("--xlim", default="", help="xlim")
parser.add_option("--ylim", default="", help="ylim")
parser.add_option("--row_cluster", type="int", default=1, help='row cluster')
parser.add_option("--col_cluster", type="int", default=1, help='col cluster')
parser.add_option("--fname_row_colors", default="", help="file name of row colors")
parser.add_option("--fname_col_colors", default="", help="file name of col colors")
parser.add_option("--prefix_table", default="", help="prefix of table for deconstructsigs")
(options, args) = parser.parse_args()
if len(args) != 3:
  parser.error("incorrect number of arguments")

type_plot=args[0]
fname_in=args[1]
fname_out=args[2]

pfilename=os.path.basename(__file__)
options.n_log_level=options.n_log_level+1
args_log=vars(options); args_log['n_log_level']=options.n_log_level+1;

if (options.n_log > 0):
  #print_log(pfilename,"%s %s %s" % (type_plot,fname_in,fname_out),args_log)
  print("%s %s %s" % (type_plot, fname_in, fname_out))

df=pd.read_table(fname_in,index_col=0)

# df.index for rows
# df.columns
# remove df.index.name
df.index.name=''
# remove df.column.name
#df.columns.name=''

# convert object to float
df=df[df.columns].astype(float) 

sns.set(font="Helvetica")
sns.set(font_scale=1.0)
sns.set_style("ticks",{"xtick.direction": u'out',"ytick.direction":u'out',"xtick.major.size": 1.0,"ytick.major.size": 1.0, "xtick.minor.size": -1.0, "ytick.minor.size": -1.0})
#sns.set_style("whitegrid")
#sns.axes_style("darkgrid"):

### color palette
# palette=pkmn_type_colors
pkmn_type_colors = ['#78C850',  # Grass
                    '#F08030',  # Fire
                    '#6890F0',  # Water
                    '#A8B820',  # Bug
                    '#A8A878',  # Normal
                    '#A040A0',  # Poison
                    '#F8D030',  # Electric
                    '#E0C068',  # Ground
                    '#EE99AC',  # Fairy
                    '#C03028',  # Fighting
                    '#F85888',  # Psychic
                    '#B8A038',  # Rock
                    '#705898',  # Ghost
                    '#98D8D8',  # Ice
                    '#7038F8',  # Dragon
                   ]


### plot
ax=plt.axes()
if type_plot=='lmplot':
  sns.lmplot(x=options.x, y=options.y, data=df, fit_reg=False, hue=options.hue)
  # hue='col3': display with different colors when df.col3 has values in [1,2,3].
elif type_plot=='distplot':
  sns.distplot(x=options.x, data=df)

elif type_plot=='countplot':
  sns.countplot(x=options.x, data=df)
  plt.xticks(rotation=-45)

elif type_plot=='barplot':
  sns.barplot(x=options.x, y=options.y, data=df, hue=options.hue)

elif type_plot=='hist':
  plt.hist(x=options.x, data=df, alpha=.3)
  sns.rugplot(x=options.x, data=df);

elif type_plot=='factorplot':
  cg = sns.factorplot(x=options.x, y=options.y, data=df, hue=options.hue, col=options.colsep, kind='swarm')
  cg.set_xticklabels(rotation=-45)

elif type_plot=='kdeplot':
  sns.kdeplot(df.y)
  #sns.kdeplot(df.col1, df.col2)

elif type_plot=='jointplot':
  # joint distribution plos combine information from scatter plots and histograms to give you detailed information for bi-variate distributions.
  # it displays pearson's corr and pvalue
  sns.jointplot(x=options.x, y=options.y, data=df)

elif type_plot=='boxplot':
  #sns.boxplot([df.y, df.x])
  sns.boxplot(data=df, palette="deep")

elif type_plot=='violinplot':
  sns.violinplot(x=options.x, y=options.y, data=df)

elif type_plot=='swarmplot':
  sns.swarmplot(x=options.x, y=options.y, data=df)
  # split=True: separte points by hue

elif type_plot=='swarm_violin':
  sns.violinplot(x=options.x, y=options.y, data=df, inner=None)
  # inner=None: remove the bars inside the violins
  sns.swarmplot(x=options.x, y=options.y, data=df, color='k', alpha=0.7) 

elif type_plot=='heatmap':
  ax = sns.heatmap(df,cmap=options.cmap,ax=ax)
  #ax=sns.heatmap(df, vmin=df.values.min(), vmax=1, square=True, cmap="YlGnBu", linewidths=0.1, annot=True, annot_kws={"size":8}) 
  #ax=sns.heatmap(correlation, vmax=1, square=True,annot=True,cmap='cubehelix',xticklabels=True,yticklabels=True)
  # turn the axis label
  #for item in ax.get_yticklabels():
  #  item.set_rotation(0)
  #for item in ax.get_xticklabels():
  #  item.set_rotation(90)
  ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
  ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

elif type_plot=='clustermap':
  # seaborn.clustermap(data, pivot_kws=None, method='average', metric='euclidean', z_score=None, standard_scale=None, figsize=None, cbar_kws=None, row_cluster=True, col_cluster=True, row_linkage=None, col_linkage=None, row_colors=None, col_colors=None, mask=None, **kwargs)

  f_row_cluster=True; f_col_cluster=True;
  if (options.row_cluster==0): f_row_cluster=False
  if (options.col_cluster==0): f_col_cluster=False
  
  # http://cloford.com/resources/colours/500col.htm
  df_row_colors=None; df_col_colors=None;
  if (len(options.fname_row_colors) > 0):
    df_row_colors=pd.read_table(options.fname_row_colors,index_col=0)
    df_row_colors.index.name=''
  if (len(options.fname_col_colors) > 0):
    df_col_colors=pd.read_table(options.fname_col_colors,index_col=0)
    df_col_colors.index.name=''

  (m1,n1)=df_row_colors.shape
  (m2,n2)=df_col_colors.shape
  if (n1 > 1) or (n2 > 1):
    sns.set(font_scale=0.6)
  cg = sns.clustermap(df.fillna(0), mask=df.isnull(), cmap=options.cmap, cbar_kws={"label": options.ylabel_cbar}, row_cluster=f_row_cluster, col_cluster=f_col_cluster, row_colors=df_row_colors, col_colors=df_col_colors)

  #cg.dendrogram_col.linkage # linkage matrix for columns
  #cg.dendrogram_row.linkage # linkage matrix for rows

elif type_plot=='subplot_sinplot':
  with sns.axes_style("darkgrid"):
    plt.subplot(211)
    sinplot()
  plt.subplot(212)
  sinplot(-1)

#sns.despine()
#sns.despine(offset=10, trim=True);

fig=plt.gcf();

#plot.axis('off')

if (len(options.xlim) > 0):
  v=options.xlim.split(',')
  try:
    xlim1=float(v[0])
  except ExceptionToCheck as e:
    xlim1=None
   
  try:
    xlim2=float(v[1])
  except ExceptionToCheck as e:
    xlim2=None
  plt.xlim(xlim1,xlim2)
  #fig.set_xlim([1,10])

if (len(options.ylim) > 0):
  v=options.ylim.split(',')
  try:
    ylim1=float(v[0])
  except ExceptionToCheck as e:
    ylim1=None

  try:
    ylim2=float(v[1])
  except ExceptionToCheck as e:
    ylim2=None
  plt.ylim(ylim1,ylim2)
  #fig.set_ylim([1,12])

# place legend to the right
#plt.legend(bbox_to_anchor=(1, 1), loc=2)

fig=plt.gcf();

if type_plot=='clustermap':

  ## autoadjust a clustermap figure size
  # matrix
  (m,n)=df.shape; size_x=0.15*n; size_y=0.15*m; 
  # ticklabels
  max_xlabel=0;
  for ticklabel in cg.ax_heatmap.xaxis.get_majorticklabels():
     len1=len(ticklabel.get_text())
     max_xlabel=max(max_xlabel,len1)
  max_ylabel=0;
  for ticklabel in cg.ax_heatmap.yaxis.get_majorticklabels():
     len1=len(ticklabel.get_text())
     max_ylabel=max(max_ylabel,len1)
  size_x += max_ylabel*0.15; size_y += max_xlabel*0.15
  # row_colors,col_colors
  (m1,n1)=df_row_colors.shape; size_x += 0.15*n1;
  (m2,n2)=df_col_colors.shape; size_y += 0.15*n2;
  # dendrogram
  size_x += 0.1; size_y += 0.1;
  # title
  if (len(options.title) > 0): size_y += 0.15;
  # expand if a figure was too small
  v=max(size_x,size_y)
  if (v < 5):
    ratio = 5/v; size_x *= ratio; size_y *= ratio;
  print('size_x={0} size_y={1}'.format(size_x,size_y))
  fig.set_size_inches(size_x, size_y)

  s_title=options.fontsize_title
  s_xlabel=options.fontsize_xlabel
  s_ylabel=options.fontsize_ylabel

  #print dir(cg)
  #print dir(cg.dim_ratios)
  #print dir(cg.ax_col_colors)
  #ax = cg.ax_heatmap 
  #ax.plot([5, 30], [5, 5], 'k-', lw = 10)

  cg.ax_col_dendrogram.set_title(options.title,fontsize=s_title)
  plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=90, fontsize=s_xlabel)
  plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, fontsize=s_ylabel )

  #fig.tight_layout() does not work for clustermap
  # The returned object has a savefig method that should be used if you want to save the figure object without clipping the dendrograms. :from the sns.clustermap documentation.
  # perfrom similar actions to tight_layout 
  cg.savefig(fname_out, dpi=100)

else:
  fig.set_size_inches(options.figsize_x, options.figsize_y)
  if (len(options.title) > 0):
    plt.title(options.title,fontsize=options.fontsize_title)
    #fig.subplots_adjust(top=-0.2)
  if (len(options.xlabel) > 0):
    plt.xlabel(options.xlabel,fontsize=options.fontsize_xlabel)
    fig.subplots_adjust(bottom=0.2)
  if (len(options.ylabel) > 0):
    plt.ylabel(options.ylabel,fontsize=options.fontsize_ylabel)
    fig.subplots_adjust(bottom=0.2)
  # this should be called after all exes have been added
  fig.tight_layout()
  # savefig
  fig.savefig(fname_out, dpi=100)


