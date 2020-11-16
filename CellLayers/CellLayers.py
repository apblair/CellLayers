import collections
import numpy as np
import pandas as pd
import scipy.io as sci
from itertools import islice

import seaborn as sns

import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import plotly
import chart_studio.plotly as py
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def CellLayers(exp_df, meta_df,
               modularity=None,
               silhouette=None,
               genes=None, exp_color=None, 
               coexpressed_genes=None, coexp_color=None, 
               tri_coexpressed_genes=None,
               edge_cutoff=None):
    while genes == None:
        print('Please add gene(s) for Cell Layers.')
        break
    else:
        return CellLayersCompute(exp_df, meta_df,
                                 modularity,
                                 silhouette,
                                 genes, exp_color,
                                 coexpressed_genes, coexp_color,
                                 tri_coexpressed_genes,
                                 edge_cutoff).compute()