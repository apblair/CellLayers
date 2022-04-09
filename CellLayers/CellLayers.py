from .BuildSankey import *
from .BuildEnrichmentSankey import *
from .BuildCoExpressionSankey import *

def build_sankey(
    exp_df, 
    meta_df, 
    genes,
    coexpressed_genes,
    exp_color,
    modularity=None,
    silhouette=None,
    node_color='#F7ED32'):
    """
    Parameters
    ----------
        exp_df: DataFrame (Pandas)
            Cell barcode x gene dataframe
        meta_df: DataFrame (Pandas)
            Cell barcode x metadata attributes
        genes: List[str]
            A list of strings that are gene names. The gene names must be present in exp_df.
        coexpressed_genes: List[List[str]] (optional)
            Nested lists of strings of length three that are gene names. The gene names must be present in exp_df.
        exp_color: str
            A string that is a matplotlib colormap continuous color.
        modularity: DataFrame (Pandas) (optional)
            Cluster resolution value by modularity score
        silhouette: DataFrame (Pandas) (optional)
            Cluster resolution and community assignment by silhouette score
        node_color: str

    """

    cl = MultiResolutionAnalysis(
        exp_df, 
        meta_df, 
        genes,
        coexpressed_genes,
        exp_color,
        modularity,
        silhouette)
    cl.compute()
    sankey_fig = BuildSankey(cl.sankey_dict, node_color=node_color).run()
    return sankey_fig, cl.sankey_dict

def build_enrichment_sankey(sankey_dict, 
                            geneset_oi,
                            enrichment_df, 
                            leading_edge, 
                            genes):
    """
    
    Parameters
    ----------
    """
    sankey_dict = EnrichmentSankey(sankey_dict,
                               geneset_oi,
                               enrichment_df,
                               leading_edge).run()
    enrichment_fig, enrichment_sankey_dict = BuildEnrichmentSankey(sankey_dict,
                                              genes,
                                              geneset_oi,
                                              enrichment_df).run()
    return enrichment_fig, enrichment_sankey_dict

def build_coexpression_sankey(sankey_dict):
    """
    """
    coexpression_sankey_fig = BuildCoExpressionSankey(sankey_dict)
    return coexpression_sankey_fig