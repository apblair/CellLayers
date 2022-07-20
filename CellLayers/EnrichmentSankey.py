from .EnrichmentAnalysis import *
class EnrichmentSankey:
    """
    Class for building a Sankey network for multi-resolution single cell enrichment analysis  

    Parameters
    ----------
    sankey_dict: dict
        Dictionary containing the multi-resolution cluster analysis for building an enrichment Sankey network
    geneset_oi: list
        List of strings that are the genesets of interest
    enrichment_df: DataFrame (Pandas)
        DataFrame containing the multi-resolution enrichment analysis

    Attributes
    ----------
    sankey_dict: dict
        Dictionary containing the multi-resolution cluster analysis for building a Sankey network
    """
    def __init__(self, 
                 sankey_dict,
                 geneset_oi,
                 genes,
                 enrichment_df
                ):
        self.sankey_dict = sankey_dict
        self.starter_gene = genes[0]
        self.geneset_oi = geneset_oi
        self.starter_geneset = geneset_oi[-1]
        self.enrichment_df = enrichment_df
    
    def _create_node_labels(self):
        """Create node hovertemplate text"""
        label_list = []
        # temporary fix
        go_term = self.starter_geneset.split(' ')[-1].replace('(', '').replace(')','')
        go_title = self.starter_geneset.split(' ')
        go_title.pop(-1)
        go_title = ' '.join(go_title)
        for data in self.sankey_dict['node_data'][['modularity', 'silhoutte_norm_by_res', self.starter_geneset + '_combined.score', 'top_genes']].values:
            label = 'Modularity Score: ' + str(data[0]) + \
                '<br />' +'Normalized Silhouette Score: ' + str(round(data[1],2)) + \
                '<br />' + go_term + ': ' + go_title + \
                '<br />' + 'Gene Set Score: ' + str(round(data[2],2)) + \
                '<br />' + 'Top Genes: ' + str(data[-1])
            label_list.append(label)
        self.sankey_dict['node_data']['label'] = label_list
    
    def _create_sankey(self):
        """Create a Sankey Plotly graph object"""
        return go.Figure(data=[dict(type='sankey', orientation='h', 
                           node = dict(pad = 10, thickness=10,
                                       label=self.sankey_dict['node_data']['node_labels'],
                                       customdata = self.sankey_dict['node_data']['label'],
                                       hovertemplate= '%{customdata}',
                                       color = self.sankey_dict['node_data'][self.starter_geneset + '_hex']),
                                    
                           link = dict(source = self.sankey_dict['data']['source'],
                                       target = self.sankey_dict['data']['target'],
                                       color = self.sankey_dict['data'][self.starter_gene+'_hex'],
                                       value = self.sankey_dict['data']['value'].tolist()),
                                       textfont=dict(color='black',size=15))])
    
    def _create_exp_colorbar(self, fig):
        """
        Create an expression color bar for a gene

        Parameters
        ----------
        fig: plotly.graph_objs._figure.Figure
            Plotly Figure subclass from the graph_objects class
        """
        fig.add_trace(go.Scatter(x=[None],
                         y=[None],
                         mode='markers',
                         visible=True,
                         marker=self.sankey_dict['exp_colorbar'][self.starter_gene]))
        
    def _create_enrichment_colorbar(self, fig):
        """
        Create an enrichment colorbar for a geneset

        Parameters
        ----------
        fig: plotly.graph_objs._figure.Figure
            Plotly Figure subclass from the graph_objects class
        """
        fig.add_trace(go.Scatter(x=[None],
                         y=[None],
                         mode='markers',
                         visible=True,
                         marker={'colorscale':self.sankey_dict['node_data'][[self.starter_geneset + '_combined.score', self.starter_geneset+'_hex']].sort_values(self.starter_geneset + '_combined.score')[self.starter_geneset+'_hex'].tolist(), 
                                 'showscale':True, 
                                 'cmin':min(self.enrichment_df[self.enrichment_df['gene.set'].isin([self.starter_geneset])]['combined.score'].tolist()),
                                 'cmax':max(self.enrichment_df[self.enrichment_df['gene.set'].isin([self.starter_geneset])]['combined.score'].tolist()), 
                                 'colorbar': {'x':1.1}}))
    
        
    def _create_dropdown_menus(self, fig):
        """
        Create dropdown menus and update figure layout

        Parameters
        ----------
        fig: plotly.graph_objs._figure.Figure
            Plotly Figure subclass from the graph_objects class
        """
        fig.update_layout({'plot_bgcolor': 'rgba(0, 0, 0, 0)',
                           'paper_bgcolor': 'rgba(0, 0, 0, 0)'})
        fig.update_xaxes(showticklabels=False) # hide all the xticks
        fig.update_yaxes(showticklabels=False) # hide all the xticks
        fig['layout']['showlegend'] = False
        fig['layout']['xaxis']['showgrid'] = False
        fig['layout']['yaxis']['showgrid'] = False
        fig.update_layout(xaxis_zeroline=False, yaxis_zeroline=False)
        fig.update_layout(
            updatemenus=[
                
                dict(y=0.9, buttons=[dict(label=self.starter_gene, 
                                           method='update', args=[{"visible":True}])],font=dict(size=15)), # temporary fix
                
                dict(y=0.5, buttons=[dict(label='GO:0002480', 
                                           method='update', args=[{"visible":True}])],font=dict(size=15)), # temporary fix

                dict(x=0, y=1.30,
                     buttons=[dict(label='Snap',method='restyle', args=['arrangement', 'snap']),
                              dict(label='Perpendicular', method='restyle',args=['arrangement', 'perpendicular']),
                              dict(label='Freeform', method='restyle',args=['arrangement', 'freeform']),
                              dict(label='Fixed', method='restyle',args=['arrangement', 'fixed'])],font=dict(size=15)),
                
                dict(x=0.2, y=1.30,
                      buttons=[dict(label='Light', method='relayout', args=['paper_bgcolor', 'white']),
                               dict(label='Dark', method='relayout', args=['paper_bgcolor', 'black'])],font=dict(size=15)),
                               
                dict(x=0.4, y=1.30,
                buttons=[dict(label='Thin', method='restyle',args=['node.thickness', 8]),
                        dict(label='Thick',method='restyle',args=['node.thickness', 15])],font=dict(size=15)),

                dict(x=0.6, y=1.30,
                      buttons=[dict(label='Small gap',method='restyle',args=['node.pad', 15]),
                               dict(label='Large gap',method='restyle',args=['node.pad', 20])],font=dict(size=15)),

                 dict(x=0.8, y=1.30,
                      buttons=[dict(label='Horizontal', method='restyle', args=['orientation', 'h']),
                               dict(label='Vertical',method='restyle',args=['orientation', 'v'])],font=dict(size=15))
                        ])
        
    def build(self):
        """Build the enrichment Sankey"""
        self._create_node_labels()
        fig = self._create_sankey()
        self._create_exp_colorbar(fig)
        self._create_enrichment_colorbar(fig)
        self._create_dropdown_menus(fig)
        return fig, self.sankey_dict