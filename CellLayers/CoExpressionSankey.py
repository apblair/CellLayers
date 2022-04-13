import numpy as np

from plotly.subplots import make_subplots
import plotly.graph_objects as go

import matplotlib
import matplotlib.colors as mcolors
class CoExpressionSankey:
    """
    Class for building a co-expression Sankey network for multi-resolution co-expression enrichment analysis  

    Parameters
    ----------
    sankey_dict: dict
        Dictionary containing the multi-resolution cluster analysis for building a Sankey network 

    Attributes
    ----------
    sankey_dict: dict
        Dictionary containing the multi-resolution cluster analysis for building a Sankey network
    """

    #TODO: Update hovertemplate; currently hardcoded
    
    def __init__(self, sankey_dict):
        self.sankey_dict = sankey_dict
        self._starting_coexpressed_genes = list(sankey_dict['coexp_color'].keys())[0]
    
    def _create_axis(self, gene, tickangle):
        """
        Create the ternary plot axis

        Parameters
        ----------
        gene: str
            String that is a gene name for labeling the ternary plot 
        tickangle: int
            Integer used to set the tick angle on the ternary plot
        """
        return {
            'title': gene,
            'titlefont': { 'size': 20 },
            'tickangle': tickangle,
            'tickfont': { 'size': 15 },
            'tickcolor': 'rgba(0,0,0,0)',
            'ticklen': 5,
            'showline': True,
            'showgrid': True}
    
    def _create_subplot(self):
        """Create a subplot for the ternary and sankey plots"""
        fig = make_subplots(rows=1, cols=2, specs=[[{"type": "ternary"},{"type": "sankey"}]])
        fig.update_xaxes(showticklabels=False) # hide all the xticks
        fig.update_yaxes(showticklabels=False) # hide all the xticks
        fig['layout']['showlegend'] = False
        fig['layout']['xaxis']['showgrid'] = False
        fig['layout']['yaxis']['showgrid'] = False
        fig.update_layout(xaxis_zeroline=False, yaxis_zeroline=False)
        fig.update_layout({
        'plot_bgcolor': 'rgba(0, 0, 0, 0)',
        'paper_bgcolor': 'rgba(0, 0, 0, 0)'})
        return fig

        
    def _create_silhouette_colorbar(self, fig):
        """
        Create a silhouette colorbar

        Parameters
        ----------
        fig: plotly.graph_objs._figure.Figure
            Plotly Figure subclass from the graph_objects class
        """
        fig.add_trace(go.Scatter(x=[None],
                   y=[None],
                   mode='markers',
                   visible=True,
                   marker={'colorscale':self.sankey_dict['silhouette_colorbar'],
                           'showscale':True, 'cmin':-1, 'cmax':1, 'colorbar': {'x':1.1}}))
    
    def _create_sankey(self, fig):
        """
        Create a Sankey Plotly graph object

        Parameters
        ----------
        fig: plotly.graph_objs._figure.Figure
            Plotly Figure subclass from the graph_objects class
        """
        self._create_silhouette_colorbar(fig)
        fig.add_trace(go.Sankey(
            node=go.sankey.Node(label= self.sankey_dict['node_labels'],
                                color = self.sankey_dict['node_data']['silhoutte_hex']),
            link=go.sankey.Link(source = self.sankey_dict['data']['source'],
                                target = self.sankey_dict['data']['target'],
                                color = [matplotlib.colors.to_hex(x) for x in self.sankey_dict['coexp_color'][self._starting_coexpressed_genes]],
                                value = self.sankey_dict['data']['value']),), row=1, col=2)
    
    def _create_ternary(self, fig):
        """
        Create a ternary Plotly graph object

        Parameters
        ----------
        fig: plotly.graph_objs._figure.Figure
            Plotly Figure subclass from the graph_objects class
        """
        fig.add_trace(go.Scatterternary({
            'mode':'markers',
            'a':np.array(self.sankey_dict['coexp_dict'][self._starting_coexpressed_genes])[:,0],
            'b':np.array(self.sankey_dict['coexp_dict'][self._starting_coexpressed_genes])[:,1],
            'c':np.array(self.sankey_dict['coexp_dict'][self._starting_coexpressed_genes])[:,2],
            'text': 'Percentile Expression',
            'hovertemplate':'MS4A1: %{a} <br> FCER1A: %{b} <br> FCGR3A: %{c}',
            'marker':{'color':self.sankey_dict['coexp_color'][self._starting_coexpressed_genes],
                      'size': np.log2(self.sankey_dict['data']['value'].tolist())}}), row=1, col=1)
        fig.update_layout(height=700, showlegend=False)
        fig.update_layout({
            'ternary':{
            'sum': 100,
            'aaxis': self._create_axis(self._starting_coexpressed_genes[0], 0),
            'baxis': self._create_axis('<br>'+self._starting_coexpressed_genes[1], 45),
            'caxis': self._create_axis('<br>'+self._starting_coexpressed_genes[2], -45)}})
        
    def build(self):
        """Build the co-expression network"""
        fig = self._create_subplot()
        self._create_sankey(fig)
        self._create_ternary(fig)
        return fig