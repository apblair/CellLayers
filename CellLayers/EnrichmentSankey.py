class EnrichmentSankey:
    def __init__(self, 
                 sankey_dict,
                 genset_oi,
                 enrichment_df,
                 leading_edge
                ):
        """
        Keyword arguments:
        - sankey_dict
        - geneset_oi
        - enrichment_df
        - leading_edge
        """
        self.sankey_dict = sankey_dict
        self.geneset_oi = geneset_oi
        self.enrichment_df = enrichment_df
        self.leading_edge = leading_edge
        
    def _create_enrichment_dict(self):
        """
        """
        # Retrieve geneset(s) of interest's enrichment scores across multi-resolutions
        enrichment_dict = {} # keys are the geneset name and values are the multi-resolution enrichment scores
        for geneset in self.geneset_oi:
            enrichment_dict[geneset] = self.enrichment_df.loc[self.enrichment_df['gene.set']== geneset, ['res_cluster', 'combined.score']]
        return enrichment_dict

    def _add_enrichment(self, enrichment_dict):
        """
        Parameters
        ----------
        enrichment_dict : dict
        """
        for geneset_name, geneset_df in enrichment_dict.items():
            geneset_df = geneset_df.rename(columns = {'combined.score':geneset_name+'_combined.score'})
            self.sankey_dict['node_data'] = pd.merge(self.sankey_dict['node_data'],
                                                     geneset_df,
                                                     left_on='node_labels',
                                                     right_on='res_cluster')
            norm = matplotlib.colors.Normalize(vmin=min(self.sankey_dict['node_data'][geneset_name+'_combined.score']),
                                               vmax=max(self.sankey_dict['node_data'][geneset_name+'_combined.score']),
                                               clip=False)
            mapper = cm.ScalarMappable(norm=norm, cmap=cm.YlGn)
            sankey_dict['node_data'][geneset_name+'_hex'] = [mcolors.to_hex(mapper.to_rgba(mapper_color)) 
                                                             for mapper_color in self.sankey_dict['node_data'][geneset_name+'_combined.score'].tolist()]
        self.sankey_dict['node_data'] = pd.merge(self.sankey_dict['node_data'],
                                                 self.leading_edge,
                                                 left_on='node_labels',
                                                 right_on='res_cluster')
        del self.sankey_dict['node_data']['res_cluster_x']
        del self.sankey_dict['node_data']['res_cluster_y']

    def run(self):
        enrichment_dict = self._create_enrichment_dict()
        self._add_enrichment(enrichment_dict)
        return sankey_fig, sankey_dict