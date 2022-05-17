from .Sankey import *
from .EnrichmentSankey import *
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
    Construct a MultiResolutionAnalysis class object with user defined parameters and build a Sankey

    Parameters
    ----------
    exp_df: DataFrame (Pandas)
        Cell barcode x gene dataframe
    meta_df: DataFrame (Pandas)
        Cell barcode x metadata attributes
    genes: List[str]
        A list of strings that are gene names. The gene names must be present in exp_df
    coexpressed_genes: List[List[str]] (optional)
        Nested lists of strings of length three that are gene names. The gene names must be present in exp_df
    exp_color: str
        A string that is a matplotlib colormap continuous color
    modularity: DataFrame (Pandas) (optional)
        Cluster resolution value by modularity score
    silhouette: DataFrame (Pandas) (optional)
        Cluster resolution and community assignment by silhouette score
    node_color: str
        String denoting the coloring of the sankey nodes
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
                            genes,
                            enrichment_df, 
                            top_genes, 
                            cmap='YlGn'
                            ):
    """
    Construct an EnrichmentAnalysis object with user defined parameters and build an enrichment sankey

    Parameters
    ----------
    sankey_dict: dict
        Dictionary containing the multi-resolution cluster analysis for building a Sankey network
    geneset_oi: list
        List of strings that are the genesets of interest
    enrichment_df: DataFrame (Pandas)
        DataFrame containing the multi-resolution enrichment analysis
    top_genes: DataFrame (Pandas)
        DataFrame containing the differentially expression genes from a multi-resolution analysis
    cmap: str (default YlGn)
        String denoting the coloring of the sankey links by the enrichment scores
    """
    ea = EnrichmentAnalysis(sankey_dict, geneset_oi, enrichment_df,top_genes, cmap)
    ea.compute()
    enrichment_fig, enrichment_sankey_dict = EnrichmentSankey(ea.sankey_dict,
                                              geneset_oi,
                                              genes,
                                              enrichment_df).build()
    return enrichment_fig, enrichment_sankey_dict

def build_coexpression_sankey(sankey_dict):
    """
    Build a co-expression sankey and ternary figure

    Parameters
    ----------
    sankey_dict: dict
        Dictionary containing the multi-resolution cluster analysis for building a Sankey network
    """
    coexpression_sankey_fig = CoExpressionSankey(sankey_dict).build()
    return coexpression_sankey_fig