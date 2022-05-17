import collections
import numpy as np
import pandas as pd

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
    meta_df: DataFrame (Pandas)
        Cell barcode x metadata attributes
    genes: List[str]
        List of strings that are gene names. The gene names must be present in exp_df.
    coexpressed_genes: List[List[str]] (optional)
        Nested lists of strings of length three that are gene names. The gene names must be present in exp_df.
    exp_color: str
        String that is a matplotlib colormap continuous color.
    modularity: DataFrame (Pandas) (optional)
        Cluster resolution value by modularity score
    silhouette: DataFrame (Pandas) (optional)
        Cluster resolution and community assignment by silhouette score

    Attributes
    ----------
    sankey_dict: dict
        Dictionary containing the multi-resolution cluster analysis for building a Sankey Network
    """
    #TODO: Add node hovertemplate text
    #TODO: Cleanup normalize silhouette
    #TODO: Add class parameters for silhouette and expression coloring

    def __init__(self, 
                 exp_df, 
                 meta_df, 
                 genes,
                 coexpressed_genes,
                 exp_color,
                 modularity,
                 silhouette): #TODO: add additional args for user metadata
        
        if genes == None:
            print('Exit Error: Please add gene(s).')

        self._exp_df = exp_df
        self._meta_df = meta_df.astype(str)
        self._build_sankey_dict(genes, exp_color, silhouette, modularity, coexpressed_genes)
      
    def _build_sankey_dict(self, genes, exp_color, silhouette, modularity, coexpressed_genes):
        """
        Build a sankey dictionary for the multi-resolution cluster analysis

        Parameters
        ----------
        genes: List[str]
            List of strings that are gene names. The gene names must be present in exp_df.
        exp_color: str
            String that is a matplotlib colormap continuous color.
        silhouette: DataFrame (Pandas) (optional)
            Cluster resolution and community assignment by silhouette score
        modularity: DataFrame (Pandas) (optional)
            Cluster resolution value by modularity score
        coexpressed_genes: List[List[str]] (optional)
            Nested lists of strings of length three that are gene names. The gene names must be present in exp_df.
        """
        
        # Check user defined gene list
        check_gene_input = list(set(genes) - set(list(self._exp_df)))
        if len(check_gene_input) > 0:
            print('Warning the following gene(s) are missing: ', *check_gene_input, sep = ", ")

        self.sankey_dict  = {'data' : None,
                    'genes' : genes,
                    'exp_dict' : {gene:[] for gene in genes},
                    'exp_color' : exp_color,
                    'exp_colorbar' : {},
                    'coexp_genes' : None,
                    'resolutions' : [data for data in list(self._meta_df) if data.split('.')[0] == 'res']}

        self._wrangle_node_data()
        self._check_args(coexpressed_genes, silhouette, modularity)
        self._count_flow_by_flow()
        self._wrangle_source_to_target()


    def _wrangle_node_data(self):
        """Create a list of resolution values and their communities"""
        self.sankey_dict['node_labels'] = [res for sublist in 
                                           [list(map(( lambda x: res + '_' + str(x)), list(set(self._meta_df[res]))))
                                            for res in self.sankey_dict['resolutions']] for res in sublist]
        self.sankey_dict['node_data'] = pd.DataFrame(self.sankey_dict['node_labels'], columns=['node_labels'])
        self.sankey_dict['node_data']['res'] = [x.split('_')[0] for x in self.sankey_dict['node_data']['node_labels'].tolist()]

    def _check_args(self, coexpressed_genes, silhouette, modularity):
        """
        Check additional user arguments and add data to sankey_dict

        Parameters
        ----------
        coexpressed_genes: List[List[str]] (optional)
            Nested lists of strings of length three that are gene names. The gene names must be present in exp_df.
        silhouette: DataFrame (Pandas) (optional)
            Cluster resolution and community assignment by silhouette score
        modularity: DataFrame (Pandas) (optional)
            Cluster resolution value by modularity score
        """

        if silhouette is not None:
            self.sankey_dict['silhouette'] = {data[0]:data[1] for data in silhouette.values}
            self.sankey_dict['node_data']['silhoutte_score'] = [round(self.sankey_dict['silhouette'][x],2) for x in self.sankey_dict['node_data']['node_labels']]
            self._create_silhouette_colorbar()
            
        if modularity is not None:
            self.sankey_dict['modularity'] = {data[0]:data[1] for data in modularity.values}
            self.sankey_dict['node_data']['modularity'] = [round(self.sankey_dict['modularity'][x],2) for x in self.sankey_dict['node_data']['res']]

        if coexpressed_genes is not None:

            check_coexp_gene_input = list(set([genes for gene_sublist in coexpressed_genes for genes in gene_sublist]) - set(list(self._exp_df)))
            if len(check_coexp_gene_input) > 0:
                print('Exiting Error: The following gene(s) are missing: ', *check_coexp_gene_input, sep = ", ")
                exit()

            self.sankey_dict['coexp_genes'] = coexpressed_genes
            self.sankey_dict['coexp_dict'] = {tuple(coexp_set): [] for coexp_set in coexpressed_genes}
            self.sankey_dict['coexp_color'] = {tuple(coexp_set): [] for coexp_set in coexpressed_genes}
    
    def _color_mapper(self, vmin, vmax, values_to_map, cmap):
        """
        Map a list of values to a hex color codes

        Parameters
        ----------
        vmin: float
            minimum value for normalization
        vmax: float
            max value for normalization
        values_to_map: list
            list of floats for normalization

        Returns
        -------
        hex_list: list[str]
            List of strings that are color hex codes
        mapper: class object
            Class object from matplotlib.cm ScalarMappable
        """
        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
        mapper = cm.ScalarMappable(norm=norm, cmap=cmap)
        rgba_list = [mapper.to_rgba(val) for val in values_to_map]
        hex_list = [mcolors.to_hex(color) for color in rgba_list]
        return hex_list, mapper

    def _normalize_silhouette(self):
        """Normalize silhouette score within each cluster resolution"""
        node_list = []
        for items in self.sankey_dict['resolutions']:
            node_temp = self.sankey_dict['node_data'][self.sankey_dict['node_data']['res'] == items].copy()
            node_temp['silhoutte_norm_by_res'] = [abs(float(i)/max(node_temp['silhoutte_score'].tolist())) 
                                                  for i in node_temp['silhoutte_score'].tolist()]
            node_list.append(node_temp)
        norm_sil_df = pd.concat(node_list)
        self.sankey_dict['node_data']['silhoutte_norm_by_res'] = norm_sil_df['silhoutte_norm_by_res']
    
    def _create_silhouette_colorbar(self):
        """Create Silhouette color bar"""
        self._normalize_silhouette()
        silhouette_list = [x/10.0 for x in list(range(-10,11,1))] # define silhouette score range
        silhouette_hex_list, mapper = self._color_mapper(min(silhouette_list), max(silhouette_list), self.sankey_dict['node_data']['silhoutte_norm_by_res'], 'RdYlBu_r')
        self.sankey_dict['node_data']['silhoutte_hex'] = silhouette_hex_list # update to modularity
        self.sankey_dict['silhouette_list'] = silhouette_list
        self.sankey_dict['silhouette_mapper'] = mapper
        self.sankey_dict['silhouette_colorbar'] = [mcolors.to_hex(self.sankey_dict['silhouette_mapper'].to_rgba(x)) 
                                         for x in self.sankey_dict['silhouette_list']]

    def _count_flow_by_flow(self):
        """"Count the number of cells flowing from resolution (n) community (i) to resolution (m) community (j)"""
        self.sankey_dict['sankey_flow_count'] = dict(collections.Counter(tuple(['_'.join([self.sankey_dict['resolutions'][res], str(dag[0])]), '_'.join([self.sankey_dict['resolutions'][res+1], str(dag[1])])])
              for res in range(0, len(self.sankey_dict['resolutions'])-1)
              for dag in self._meta_df[[self.sankey_dict['resolutions'][res], self.sankey_dict['resolutions'][res+1]]].values))
    
    def _wrangle_source_to_target(self):
        """Create a source to target mapping for each resolution (n) community (i) to resolution (m) community (j)"""
        self.sankey_dict['data'] = pd.DataFrame([[res[0], res[1], 
                                                self.sankey_dict['node_labels'].index(res[0]), 
                                                self.sankey_dict['node_labels'].index(res[1]), cell_count] 
                                                for res, cell_count in self.sankey_dict['sankey_flow_count'].items()],
                                                columns=['source_label', 'target_label', 'source', 'target', 'value'])
        
        self.sankey_dict['data']['source_res'] = np.array([x.split('_') for x in self.sankey_dict['data']['source_label']])[:,0].tolist()
        self.sankey_dict['data']['source_cluster'] = np.array([x.split('_') for x in self.sankey_dict['data']['source_label']])[:,1].tolist()
        
        self.sankey_dict['data']['target_res'] = np.array([x.split('_') for x in self.sankey_dict['data']['target_label']])[:,0].tolist()
        self.sankey_dict['data']['target_cluster'] = np.array([x.split('_') for x in self.sankey_dict['data']['target_label']])[:,1].tolist()
    
    def _create_expression_colorbar(self):
        """Create a gene expression bar"""
        # Create a hex code color range for a gene's expression in the sankey network streams
        for gene, expression in self.sankey_dict['exp_dict'].items():
            exp_hex_list, mapper = self._color_mapper(min(expression), max(expression), expression, 'Purples')
            self.sankey_dict['data'][gene+'_hex'] = exp_hex_list
            self._sort_hex_colorbar(gene, expression, exp_hex_list)
    
    def _sort_hex_colorbar(self, gene, expression, hex_list):
        """
        Sort the hex color code by gene expression in ascending order

        Parameters
        ----------
        gene: str
            String that is a gene name
        expression: list
            List of floats for each cell's gene expression value in the sankey stream
        hex_list: list
            List of strings that are hex codes
        """
        gene_color_bar_df = pd.DataFrame(self.sankey_dict['exp_dict'][gene])
        gene_color_bar_df['hex'] = hex_list
        gene_color_bar_df.sort_values(0, inplace=True)
        self.sankey_dict['exp_colorbar'][gene] = dict(colorscale=gene_color_bar_df['hex'].tolist(), showscale=True, cmin=min(expression), cmax=max(expression))
            
    def _coexpression(self, cell_ids):
        """
        Compute the co-expression of user defined genes for each batch of cell ids

        Parameters
        ----------
        cell_ids : list
            List of strings that are cell barcodes
        """
        #TODO: Check code
        for genes in self.sankey_dict['coexp_genes']:
            # gene_sums = self._exp_df.loc[cell_ids, genes].sum(axis = 1).tolist()
            gene_sum = sum(self._exp_df.loc[cell_ids][genes].sum().tolist())
            # if sum(gene_sums) > 0:
            if gene_sum > 0:
                dec_percentile = [x/gene_sum for x in self._exp_df.loc[cell_ids][genes].sum().tolist()]
                # dec_percentile = [x/sum(gene_sums) for x in gene_sums]
                self.sankey_dict['coexp_dict'][tuple(genes)].append(dec_percentile)
                self.sankey_dict['coexp_color'][tuple(genes)].append(matplotlib.colors.to_hex(dec_percentile))
            else:
                pass

    def _avg_expression(self, cell_ids):
        """
        Compute the average expression of a gene for a subset of cell ids
        
        Parameters
        ----------
        cell_ids: list
            List of strings to be spliced from expression dataframe
        """
        for gene in self.sankey_dict['genes']:
            self.sankey_dict['exp_dict'][gene].append(self._exp_df.loc[cell_ids][gene].mean())

    def compute(self):
        """Compute the average expression of the user defined genes for each flow transition"""
        for cluster_data in self.sankey_dict['data'][['source_res', 'source_cluster', 'target_res', 'target_cluster']].values:
            cell_ids = self._meta_df[(self._meta_df[cluster_data[0]]==cluster_data[1]) & (self._meta_df[cluster_data[2]]==cluster_data[3])].index.tolist()
            self._avg_expression(cell_ids)
            if self.sankey_dict['coexp_genes'] is not None:
                self._coexpression(cell_ids)
        self._create_expression_colorbar()