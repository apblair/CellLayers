
from .CellLayersCompute import CellLayersCompute
from .CellLayersConfig import CellLayersConfig


def make_config(exp_df, meta_df,
               modularity=None,
               silhouette=None,
               genes=None, exp_color=None, 
               coexpressed_genes=None, coexp_color=None, 
               tri_coexpressed_genes=None,
               edge_cutoff=None):
    while genes == None:
        print('Please add gene(s) for Cell Layers.')
        break
    else:
        return CellLayersCompute(exp_df, meta_df,
                                 modularity,
                                 silhouette,
                                 genes, exp_color,
                                 coexpressed_genes, coexp_color,
                                 tri_coexpressed_genes,
                                 edge_cutoff).compute()

def make_sankey(sankey_dict, node_color='#F7ED32'):
    return CellLayersConfig(sankey_dict, node_color=node_color)