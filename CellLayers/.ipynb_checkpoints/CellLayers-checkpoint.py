from .BuildSankeyNetwork import *

def run(exp_df, meta_df,
               modularity=None,
               silhouette=None,
               genes=None, exp_color=None, 
               coexpressed_genes=None, coexp_color=None, 
               tri_coexpressed_genes=None,
               edge_cutoff=None,
       node_color='#F7ED32'):
    while genes == None:
        print('Please add gene(s) for Cell Layers.')
        break
    else:
        sankey_dict = MultiResolutionAnalysis(exp_df, meta_df,
                                 modularity,
                                 silhouette,
                                 genes, exp_color,
                                 coexpressed_genes, coexp_color,
                                 tri_coexpressed_genes,
                                 edge_cutoff).compute()
        return BuildSankeyNetwork(sankey_dict, node_color=node_color).run()
    