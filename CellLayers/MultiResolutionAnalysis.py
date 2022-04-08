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

class MultiResolutionAnalysis:
    """
    Class for multi-resolution single cell analysis

    Parameters
    ----------
        exp_df: DataFrame (Pandas)
            Cell barcode x gene dataframe
        meta_df : DataFrame (Pandas)
            Cell barcode x metadata attributes
        modularity: DataFrame (Pandas)
            Cluster resolution value by modularity score
        silhouette: DataFrame (Pandas)
            Cluster resolution and community assignment by silhouette score
        genes: List[str]
            A list of strings that are gene names. The gene names must be present in exp_df.
        exp_color: str

        coexpressed_genes: List[List[str]]
            Nested lists of strings of le
    """
    def __init__(self, 
                 exp_df, 
                 meta_df, 
                 modularity,
                 silhouette,
                 genes, 
                 exp_color,
                 coexpressed_genes):
        
        self._exp_df = exp_df
        self._meta_df = meta_df.astype(str)
        self._sankey_dict = {'data':None,
                            'genes' : genes,
                            'exp_dict' : {gene:[] for gene in genes},
                            'exp_color' : exp_color,
                            'exp_colorbar':{},
                            'resolutions':[data for data in list(self._meta_df) if data.split('.')[0] == 'res']}
        self._check_args(silhouette, modularity, coexpressed_genes)
        self._create_node_labels()
        self._count_flow_by_flow()
        self._wrangle_source_to_target()
        
    def _check_args(self, silhouette, modularity, coexpressed_genes):
        """Check additional user arguments and add data to sankey_dict"""

        if silhouette is not None:
            self._sankey_dict['silhouette'] = {data[0]:data[1] for data in silhouette.values}
            
        if modularity is not None:
            self._sankey_dict['modularity'] = {data[0]:data[1] for data in modularity.values}
            
        if coexpressed_genes is not None:
            self._sankey_dict['coexp_genes'] = coexpressed_genes
            self._sankey_dict['coexp_dict'] = {tuple(coexp_set): [] for coexp_set in coexpressed_genes}
            self._sankey_dict['coexp_color'] = {tuple(coexp_set): [] for coexp_set in coexpressed_genes}

    def _create_node_labels(self):
        """Create a list of resolution values and their communities"""
        self._sankey_dict['node_labels'] = [res for sublist in 
                                           [list(map(( lambda x: res + '_' + str(x)), list(set(self._meta_df[res]))))
                                            for res in self._sankey_dict['resolutions']] for res in sublist]

    def _count_flow_by_flow(self):
        """"Count the number of cells flowing from resolution (n) community (i) to resolution (m) community (j)"""
        self._sankey_dict['sankey_flow_count'] = dict(collections.Counter(tuple(['_'.join([self._sankey_dict['resolutions'][res], str(dag[0])]), '_'.join([self._sankey_dict['resolutions'][res+1], str(dag[1])])])
              for res in range(0, len(self._sankey_dict['resolutions'])-1)
              for dag in self._meta_df[[self._sankey_dict['resolutions'][res], self._sankey_dict['resolutions'][res+1]]].values))
    
    def _wrangle_source_to_target(self):
        """Create a source to target mapping for each resolution (n) community (i) to resolution (m) community (j)"""
        self._sankey_dict['data'] = pd.DataFrame([[res[0], res[1], 
                                                self._sankey_dict['node_labels'].index(res[0]), 
                                                self._sankey_dict['node_labels'].index(res[1]), cell_count] 
                                                for res, cell_count in self._sankey_dict['sankey_flow_count'].items()],
                                                columns=['source_label', 'target_label', 'source', 'target', 'value'])
        
        self._sankey_dict['data']['source_res'] = np.array([x.split('_') for x in self._sankey_dict['data']['source_label']])[:,0].tolist()
        self._sankey_dict['data']['source_cluster'] = np.array([x.split('_') for x in self._sankey_dict['data']['source_label']])[:,1].tolist()
        
        self._sankey_dict['data']['target_res'] = np.array([x.split('_') for x in self._sankey_dict['data']['target_label']])[:,0].tolist()
        self._sankey_dict['data']['target_cluster'] = np.array([x.split('_') for x in self._sankey_dict['data']['target_label']])[:,1].tolist()

        
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
    
    def _avg_expression(self, cell_ids):
        """
        Compute the average expression of a gene for a subset of cell ids
        
        Parameters
        ----------
        cell_ids : list
            List of strings to be spliced from expression dataframe
        """
        for gene in self._sankey_dict['genes']:
            self._sankey_dict['exp_dict'][gene].append(self._exp_df.loc[cell_ids][gene].mean())
    
    def _create_expression_hex_color(self):
        """
        Create a hex code color range for a gene's expression in the sankey network streams
        """
        for gene, expression in self._sankey_dict['exp_dict'].items():
            norm = matplotlib.colors.Normalize(vmin=min(expression), vmax=max(expression), clip=True)
            mapper = cm.ScalarMappable(norm=norm, cmap=cm.Purples)
            rgba_list = [mapper.to_rgba(exp) for exp in expression]
            hex_list = [mcolors.to_hex(color) for color in rgba_list]
            self._sankey_dict['data'][gene+'_hex'] = hex_list
            self._create_expression_colorbar(gene, expression, hex_list)
    
    def _create_expression_colorbar(self, gene, expression, hex_list):
        """
        Create a gene expression bar

        Parameters
        ----------
        gene : str
            A string that is a gene name
        expression : list
            A list of floats for each cell's gene expression value in the sankey stream
        hex_list : list
            A list of strings that are hex codes
        """
        gene_color_bar_df = pd.DataFrame(self._sankey_dict['exp_dict'][gene])
        gene_color_bar_df['hex'] = hex_list
        gene_color_bar_df.sort_values(0, inplace=True)
        self._sankey_dict['exp_colorbar'][gene] = dict(colorscale=gene_color_bar_df['hex'].tolist(), showscale=True, cmin=min(expression), cmax=max(expression))
            
    def _coexpression(self, cell_ids):
        """
        Parameters
        ----------
        cell_ids : list
            A list of strings that are cell barcodes
        """
        for genes in self._sankey_dict['coexp_genes']:
            total_sum = sum(self._exp_df.loc[cell_ids][genes].sum().tolist())
            if total_sum > 0:
                dec_percentile = [x/total_sum for x in self._exp_df.loc[cell_ids][genes].sum().tolist()]
                self._sankey_dict['coexp_dict'][tuple(genes)].append(dec_percentile)
                self._sankey_dict['coexp_color'][tuple(genes)].append(matplotlib.colors.to_hex(dec_percentile))
            else:
                pass
    
    def _compute_genomic_data(self):
        """Compute the average expression of the user defined genes for each flow transition"""
        for cluster_data in self._sankey_dict['data'][['source_res', 'source_cluster', 'target_res', 'target_cluster']].values:
            cell_ids = self._meta_df[(self._meta_df[cluster_data[0]]==cluster_data[1]) & (self._meta_df[cluster_data[2]]==cluster_data[3])].index.tolist()
            self._avg_expression(cell_ids)
            # TODO: Add coexpression check
            self._coexpression(cell_ids)
        self._create_expression_hex_color()
    
    def _create_node_data(self):
        """ """
        node_df = pd.DataFrame(self._sankey_dict['node_labels'], columns=['node_labels'])
        node_df['res'] = [x.split('_')[0] for x in node_df['node_labels'].tolist()]
        node_df['modularity'] = [round(self._sankey_dict['modularity'][x],2) for x in node_df['res']]
        node_df['silhoutte_score'] = [round(self._sankey_dict['silhouette'][x],2) for x in node_df['node_labels']]
        node_df['new_label'] = ['Modularity '+str(items[-2:][0]) + '<br />' + 'Silhouette ' + str(items[-2:][1]) for items in node_df.values]
        self._sankey_dict['node_data']=node_df
        self._normalize_silhouette()
        self._color_silhouette()
        
    def _normalize_silhouette(self):
        """
        """
        node_list = []
        for items in self._sankey_dict['resolutions']:
            node_temp = self._sankey_dict['node_data'][self._sankey_dict['node_data']['res'] == items].copy()
            node_temp['silhoutte_norm_by_res'] = [abs(float(i)/max(node_temp['silhoutte_score'].tolist())) 
                                                  for i in node_temp['silhoutte_score'].tolist()]
            node_list.append(node_temp)
        norm_sil_df = pd.concat(node_list)
        self._sankey_dict['node_data']['silhoutte_norm_by_res'] = norm_sil_df['silhoutte_norm_by_res']
    
    def _color_silhouette(self):
        """
        """
        silhouette_list = [x/10.0 for x in list(range(-10,11,1))]
        norm = matplotlib.colors.Normalize(vmin=min(silhouette_list), vmax=max(silhouette_list), clip=True) # renornmalize per module
        mapper = cm.ScalarMappable(norm=norm, cmap=cm.RdYlBu_r)
        rgba_list = [mapper.to_rgba(exp) for exp in self._sankey_dict['node_data']['silhoutte_norm_by_res'].tolist()]
        self._sankey_dict['node_data']['silhoutte_hex'] = [mcolors.to_hex(color) for color in rgba_list] # update to modularity
        self._sankey_dict['silhouette_list'] = silhouette_list
        self._sankey_dict['silhouette_mapper'] = mapper
    
    def compute(self):
        """
        
        """
        self._compute_genomic_data()
        self._create_node_data()
        return self._sankey_dict