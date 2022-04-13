from .Sankey import *
from .BuildEnrichmentSankey import *
from .CoExpressionSankey import *

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
    Construct a MultiResolutionAnalysis class object with user defined parameters and build a Sankey plotly graph object.

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
    mra = MultiResolutionAnalysis(
        exp_df, 
        meta_df, 
        genes,
        coexpressed_genes,
        exp_color,
        modularity,
        silhouette)
    mra.compute()
    sankey_fig = Sankey(mra.sankey_dict, node_color=node_color).build()
    return sankey_fig, mra.sankey_dict

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
    Parameters
    ----------
    """
    coexpression_sankey_fig = CoExpressionSankey(sankey_dict).build()
    return coexpression_sankey_fig