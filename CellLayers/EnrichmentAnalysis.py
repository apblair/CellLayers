from .MultiResolutionAnalysis import *

class EnrichmentAnalysis:
    """

    Parameters
    ----------
    sankey_dict: dict
        Dictionary containing the multi-resolution cluster analysis for building a Sankey Network
    geneset_oi:
    enrichment_df:
    leading_edge:

    Attributes
    ----------
    sankey_dict: dict
        Dictionary containing the multi-resolution cluster analysis for building a Sankey Network
    """
    #TODO: Add test

    def __init__(self, 
                 sankey_dict,
                 geneset_oi,
                 enrichment_df,
                 leading_edge
                ):
        self.sankey_dict = sankey_dict
        self._geneset_oi = geneset_oi
        self._enrichment_df = enrichment_df
        self._leading_edge = leading_edge
        
    def _create_enrichment_dict(self) -> dict:
        """Create a dictionary for the multi-resolution single cell enrichment analysis"""
        # Retrieve geneset(s) of interest's enrichment scores across multi-resolutions
        enrichment_dict = {} # keys are the geneset name and values are the multi-resolution enrichment scores
        for geneset in self._geneset_oi:
            enrichment_dict[geneset] = self._enrichment_df.loc[self._enrichment_df['gene.set']== geneset, ['res_cluster', 'combined.score']]
        return enrichment_dict

    def _enrichment(self, enrichment_dict):
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
            self.sankey_dict['node_data'][geneset_name+'_hex'] = [mcolors.to_hex(mapper.to_rgba(mapper_color)) 
                                                             for mapper_color in self.sankey_dict['node_data'][geneset_name+'_combined.score'].tolist()]
        self.sankey_dict['node_data'] = pd.merge(self.sankey_dict['node_data'],
                                                 self._leading_edge,
                                                 left_on='node_labels',
                                                 right_on='res_cluster')
        del self.sankey_dict['node_data']['res_cluster_x']
        del self.sankey_dict['node_data']['res_cluster_y']

    def compute(self):
        enrichment_dict = self._create_enrichment_dict()
        self._enrichment(enrichment_dict)