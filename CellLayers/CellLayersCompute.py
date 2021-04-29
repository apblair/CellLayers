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

class CellLayersCompute:
    def __init__(self, 
                 exp_df, meta_df, 
                 modularity,
                 silhouette,
                 genes, exp_color,
                 coexpressed_genes, coexp_color,
                 tri_coexpressed_genes,
                 edge_cutoff):
        
        self.exp_df = exp_df
        self.meta_df = meta_df.astype(str)
        self.sankey_dict = {
            'data':None,
            'sankey_flow_count': None, 
            
            'genes' : genes,
            'exp_dict' : {gene:[] for gene in genes},
            'exp_color' : exp_color,
            'exp_colorbar':{},

            'edge_cutoff':edge_cutoff,
            'resolutions':[data for data in list(self.meta_df) if data.split('.')[0] == 'res']}
        
        if silhouette is not None:
            self.sankey_dict['silhouette'] = {data[0]:data[1] for data in silhouette.values}

        if modularity is not None:
            self.sankey_dict['modularity'] = {data[0]:data[1] for data in modularity.values}
            
        if coexpressed_genes is not None:
            self.sankey_dict['coexp_genes'] = coexpressed_genes
            self.sankey_dict['coexp_dict'] = {tuple(coexp_set): [] for coexp_set in coexpressed_genes}
            self.sankey_dict['coexp_color'] = coexp_color
            self.sankey_dict['coexp_colorbar'] = {}
        
            
        if tri_coexpressed_genes is not None:
            self.sankey_dict['tri_coexp_genes'] = tri_coexpressed_genes
            self.sankey_dict['tri_coexp_dict'] = {tuple(tri_coexp_set): [] for tri_coexp_set in tri_coexpressed_genes}
            self.sankey_dict['tri_coexp_color'] = {tuple(tri_coexp_set): [] for tri_coexp_set in tri_coexpressed_genes}
        
        
            
        self.sankey_dict['node_labels'] = [res for sublist in 
                                           [list(map(( lambda x: res + '_' + str(x)), list(set(self.meta_df[res]))))
                                            for res in self.sankey_dict['resolutions']] for res in sublist]
        
        self._count_flow_by_flow()
        
        self.sankey_dict['data'] = pd.DataFrame([[res[0], res[1], 
                                                  self.sankey_dict['node_labels'].index(res[0]), 
                                                  self.sankey_dict['node_labels'].index(res[1]), cell_count] 
                                                 for res, cell_count in self.sankey_dict['sankey_flow_count'].items()],
                                                
                                                columns=['source_label', 'target_label', 'source', 'target', 'value'])
        
        self.sankey_dict['data']['source_res'] = np.array([x.split('_') for x in self.sankey_dict['data']['source_label']])[:,0].tolist()
        self.sankey_dict['data']['source_cluster'] = np.array([x.split('_') for x in self.sankey_dict['data']['source_label']])[:,1].tolist()
        
        self.sankey_dict['data']['target_res'] = np.array([x.split('_') for x in self.sankey_dict['data']['target_label']])[:,0].tolist()
        self.sankey_dict['data']['target_cluster'] = np.array([x.split('_') for x in self.sankey_dict['data']['target_label']])[:,1].tolist()
        
    def _count_flow_by_flow(self):
        """"Count the number of cells flowing from resolution (n) community (i) to resolution (m) community (j)"""
        self.sankey_dict['sankey_flow_count'] = dict(collections.Counter(tuple(['_'.join([self.sankey_dict['resolutions'][res], str(dag[0])]), 
                                                                                '_'.join([self.sankey_dict['resolutions'][res+1], str(dag[1])])])
              for res in range(0, len(self.sankey_dict['resolutions'])-1)
              for dag in self.meta_df[[self.sankey_dict['resolutions'][res], self.sankey_dict['resolutions'][res+1]]].values))
        
    def _follow_the_flow(self, seq, n=2):
        """
        Yield a sliding window of width n over an iterable 
        
        Parameters
        ----------
        seq : list
            List of strings to be sliced
        n : int, default=2
            Window size slice
        
        Yields
        -------
        result : list
            List of sliced strings of length n
        """
        it = iter(seq)
        result = tuple(islice(it, n))
        if len(result) == n:
            yield result
        for elem in it:
            result = result[1:] + (elem,)
            yield result 
            
    def _count_flow_by_flow(self):
        """"Count the number of cells flowing from resolution (n) community (i) to resolution (m) community (j)"""
        self.sankey_dict['sankey_flow_count'] = dict(collections.Counter(tuple(['_'.join([self.sankey_dict['resolutions'][res], str(dag[0])]), 
                                                                                '_'.join([self.sankey_dict['resolutions'][res+1], str(dag[1])])])
              for res in range(0, len(self.sankey_dict['resolutions'])-1)
              for dag in self.meta_df[[self.sankey_dict['resolutions'][res], self.sankey_dict['resolutions'][res+1]]].values))
    
    def _avg_expression(self, cell_ids):
        for gene in self.sankey_dict['genes']:
            self.sankey_dict['exp_dict'][gene].append(self.exp_df.loc[cell_ids][gene].mean())
    
    def _create_expression_hex_color(self):
        for k,v in self.sankey_dict['exp_dict'].items():
            norm = matplotlib.colors.Normalize(vmin=min(v), vmax=max(v), clip=True)
            mapper = cm.ScalarMappable(norm=norm, cmap=cm.Purples)
            rgba_list = [mapper.to_rgba(exp) for exp in v]
            hex_list = [mcolors.to_hex(color) for color in rgba_list]
            self.sankey_dict['data'][k+'_hex'] = hex_list
            self._create_expression_colorbar(k,v, hex_list)
    
    def _create_expression_colorbar(self, k,v, hex_list):
            color_bar_df = pd.DataFrame(self.sankey_dict['exp_dict'][k])
            color_bar_df['hex'] = hex_list
            color_bar_df = color_bar_df.sort_values(0)
            self.sankey_dict['exp_colorbar'][k] = dict(colorscale=color_bar_df['hex'].tolist(),showscale=True, cmin=min(v),cmax=max(v))
            
    def _tri_coexpression(self, cell_ids):
        for genes in self.sankey_dict['tri_coexp_genes']:
            total_sum = sum(self.exp_df.loc[cell_ids][genes].sum().tolist())
            if total_sum > 0:
                dec_percentile = [x/total_sum for x in self.exp_df.loc[cell_ids][genes].sum().tolist()]
                self.sankey_dict['tri_coexp_dict'][tuple(genes)].append(dec_percentile)
                self.sankey_dict['tri_coexp_color'][tuple(genes)].append(matplotlib.colors.to_hex(dec_percentile))
                
            else:
                pass
    
    def _compute_genomic_data(self):
        """Compute the average expression of the user defined genes for each flow transition"""
        for cluster_data in self.sankey_dict['data'][['source_res', 'source_cluster', 'target_res', 'target_cluster']].values:
            cell_ids = self.meta_df[(self.meta_df[cluster_data[0]]==cluster_data[1]) & (self.meta_df[cluster_data[2]]==cluster_data[3])].index.tolist()
            self._avg_expression(cell_ids)
#             self._tri_coexpression(cell_ids)
        self._create_expression_hex_color()
    
    def _create_node_data(self):
        node_df = pd.DataFrame(self.sankey_dict['node_labels'], columns=['node_labels'])
        node_df['res'] = [x.split('_')[0] for x in node_df['node_labels'].tolist()]
        node_df['modularity'] = [round(self.sankey_dict['modularity'][x],2) for x in node_df['res']]
        node_df['silhoutte_score'] = [round(self.sankey_dict['silhouette'][x],2) for x in node_df['node_labels']]
        node_df['new_label'] = ['Modularity '+str(items[-2:][0]) + '<br />' + 'Silhouette ' + str(items[-2:][1]) for items in node_df.values]
        self.sankey_dict['node_data']=node_df
        self._normalize_silhouette()
        self._color_silhouette()
        
    def _normalize_silhouette(self):
        node_list = []
        for items in self.sankey_dict['resolutions']:
            node_temp = self.sankey_dict['node_data'][self.sankey_dict['node_data']['res'] == items]
            node_temp['silhoutte_norm_by_res'] = [abs(float(i)/max(node_temp['silhoutte_score'].tolist())) 
                                                  for i in node_temp['silhoutte_score'].tolist()]
            node_list.append(node_temp)
        norm_sil_df = pd.concat(node_list)
        self.sankey_dict['node_data']['silhoutte_norm_by_res'] = norm_sil_df['silhoutte_norm_by_res']
    
    def _color_silhouette(self):
        sil_list = [x/10.0 for x in list(range(-10,11,1))]
        norm = matplotlib.colors.Normalize(vmin=min(sil_list), vmax=max(sil_list), clip=True) # renornmalize per module
        mapper = cm.ScalarMappable(norm=norm, cmap=cm.RdYlBu_r)
        rgba_list = [mapper.to_rgba(exp) for exp in self.sankey_dict['node_data']['silhoutte_norm_by_res'].tolist()]
        self.sankey_dict['node_data']['silhoutte_hex'] = [mcolors.to_hex(color) for color in rgba_list] # update to modularity
        self.sankey_dict['silhouette_list'] = sil_list
        self.sankey_dict['silhouette_mapper'] = mapper
    
    def compute(self):
        self._compute_genomic_data()
        self._create_node_data()
        return self.sankey_dict