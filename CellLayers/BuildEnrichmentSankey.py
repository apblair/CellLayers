class BuildEnrichmentSankey:
    def __init__(self, 
                 sankey_dict,
                 starter_gene,
                 geneset_oi,
                 enrichment_df
                ):
                """
        Keyword arguments:
        - sankey_dict
        - starter_gene
        - enrichment_df
        - enrichment_df
        """
        self.sankey_dict = sankey_dict
        self.starter_gene = starter_gene
        self.geneset_oi = geneset_oi
        self.starter_geneset = geneset_oi[-1]
        self.enrichment_df = enrichment_df
    
    def _create_node_label(self):
        """
        """
        label_list = []
        for data in self.sankey_dict['node_data'][['modularity', 'silhoutte_score', self.starter_geneset + '_combined.score', 'top_genes']].values:
            label = '<br />' + self.starter_geneset + '<br />' + 'Gene Set Score: ' + str(data[2]) + '<br />' + 'Top Genes: ' + str(data[-1])
            label_list.append(label)
        self.sankey_dict['node_data']['label'] = label_list
    
    def _create_sankey_obj(self):
        """
        """
        return go.Figure(data=[dict(type='sankey', orientation='h', 
                           node = dict(pad = 10, thickness=10,
                                       label=self.sankey_dict['node_data']['node_labels'],
                                       customdata = self.sankey_dict['node_data']['label'],
                                       hovertemplate= '%{customdata}',
                                       color = self.sankey_dict['node_data'][self.starter_geneset + '_hex']),
                                    
                           link = dict(source = self.sankey_dict['data']['source'],
                                       target = self.sankey_dict['data']['target'],
                                       color = self.sankey_dict['data'][self.starter_gene+'_hex'],
                                       value = self.sankey_dict['data']['value']))])
    
    def _create_exp_color_bar(self, fig):
        """
        Parameters
        ----------
        """
        fig.add_trace(go.Scatter(x=[None],
                         y=[None],
                         mode='markers',
                         visible=True,
                         marker=self.sankey_dict['exp_colorbar'][self.starter_gene]))
        
    def _create_enrichment_color_bar(self, fig):
        """
        Parameters
        ----------
        """
        fig.add_trace(go.Scatter(x=[None],
                         y=[None],
                         mode='markers',
                         visible=True,
                         marker={'colorscale':self.sankey_dict['node_data'][[self.starter_geneset + '_combined.score', self.starter_geneset+'_hex']].sort_values(self.starter_geneset + '_combined.score')[self.starter_geneset+'_hex'].tolist(), 
                                 'showscale':True, 
                                 'cmin':min(pbmc_gsea[pbmc_gsea['gene.set'].isin([self.starter_geneset])]['combined.score'].tolist()),
                                 'cmax':max(pbmc_gsea[pbmc_gsea['gene.set'].isin([self.starter_geneset])]['combined.score'].tolist()), 
                                 'colorbar': {'x':1.1}}))
    
        
    def _create_dropdown_menus(self, fig):
        """
        Parameters
        ----------
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
                dict(y=-0.2,
                     buttons=[dict(label='Snap',method='restyle', args=['arrangement', 'snap']),
                              dict(label='Perpendicular', method='restyle',args=['arrangement', 'perpendicular']),
                              dict(label='Freeform', method='restyle',args=['arrangement', 'freeform']),
                              dict(label='Fixed', method='restyle',args=['arrangement', 'fixed'])]),
                dict(y=0,
                      buttons=[dict(label='Small gap',method='restyle',args=['node.pad', 15]),
                               dict(label='Large gap',method='restyle',args=['node.pad', 20])]),
                 dict(y=0.2,
                      buttons=[dict(label='Horizontal', method='restyle', args=['orientation', 'h']),
                               dict(label='Vertical',method='restyle',args=['orientation', 'v'])]),
                 dict(y=0.4,
                      buttons=[dict(label='Light', method='relayout', args=['paper_bgcolor', 'white']),
                               dict(label='Dark', method='relayout', args=['paper_bgcolor', 'black'])])])
        
    def run(self):
        """
        """
        self._create_node_label()
        fig = self._create_sankey_obj()
        self._create_exp_color_bar(fig)
        self._create_enrichment_color_bar(fig)
        self._create_dropdown_menus(fig)
        return fig, self.sankey_dict