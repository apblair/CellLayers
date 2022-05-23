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
    #TODO: Update gene, cluster evaluation, and go bpa buttons

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
            'titlefont': { 'size': 30 },
            'tickangle': tickangle,
            'tickfont': { 'size': 30 },
            'tickcolor': 'rgba(0,0,0,0)',
            'ticklen': 5,
            'showline': True,
            'showgrid': True}
    
    def _create_subplot(self):
        """Create a subplot for the ternary and sankey plots"""
        fig = make_subplots(rows=2, cols=2, specs=[
            [{"type": "ternary", "rowspan":2},{"type": "sankey","rowspan":2}],
            [{}, {}]
            ])
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
                           'showscale':True, 'cmin':-1, 'cmax':1, 'colorbar': {'x':1.1, 'tickfont':{'size':10}},
                           }
                           ))
    
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
                                value = self.sankey_dict['data']['value']),textfont=dict(color='black',size=20)), row=1, col=2)
    
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
            'hovertemplate':self._starting_coexpressed_genes[0]+': %{a} <br>'+self._starting_coexpressed_genes[1]+': %{b} <br>'+self._starting_coexpressed_genes[2]+': %{c}',
            'marker':{'color':self.sankey_dict['coexp_color'][self._starting_coexpressed_genes],
                      'size': 2 * np.log2(self.sankey_dict['data']['value'].tolist())}}), row=1, col=1)
        fig.update_layout(height=700, showlegend=False)
        fig.update_layout({
            'ternary':{
            'sum': 100,
            'aaxis': self._create_axis('<i>'+self._starting_coexpressed_genes[0]+'</i>', 0),
            'baxis': self._create_axis('<br>'+'<i>'+self._starting_coexpressed_genes[1]+'</i>', 45),
            'caxis': self._create_axis('<br>'+'<i>'+self._starting_coexpressed_genes[2]+'</i>', -45)}})
    

    def _add_functionality(self, fig):
        """
        Add drop down menus for selecting genes and configuring the layout

        Parameters
        ----------
        fig: plotly.graph_objs._figure.Figure
            Plotly Figure subclass from the graph_objects class
        """
        gene_buttons = [dict(label=','.join(str(gene) for gene in self._starting_coexpressed_genes), 
        method='update', 
        args=[{"visible":[True]}])] # temporary fix
        cluster_evaluation_buttons = [dict(label='Silhouette Scores', method='update', args=[{"visible":[True]}])] # temporary fix

        fig.update_layout(
            updatemenus=[

                dict(x=0, y=1.05, buttons=list(gene_buttons),font=dict(size=20)), # temporary fix
                dict(x=0, y=.9, buttons=list(cluster_evaluation_buttons),font=dict(size=20)), # temporary fix

                dict(x=0, y=1.25,
                buttons=[dict(label='Snap',method='restyle', args=['arrangement', 'snap']),
                        dict(label='Perpendicular', method='restyle',args=['arrangement', 'perpendicular']),
                        dict(label='Freeform', method='restyle',args=['arrangement', 'freeform']),
                        dict(label='Fixed', method='restyle',args=['arrangement', 'fixed'])], font=dict(size=20)),    
                
                dict(x=0.2, y=1.25,
                buttons=[dict(label='Light', method='relayout', args=['paper_bgcolor', 'white']),
                        dict(label='Dark', method='relayout', args=['paper_bgcolor', 'black'])], font=dict(size=20)),
                
                dict(x=0.4, y=1.25,
                buttons=[dict(label='Thin', method='restyle',args=['node.thickness', 8]),
                        dict(label='Thick',method='restyle',args=['node.thickness', 15])], font=dict(size=20)),
                
                dict(x=0.6, y=1.25,
                buttons=[dict(label='Small gap',method='restyle',args=['node.pad', 15]),
                        dict(label='Large gap',method='restyle',args=['node.pad', 20])], font=dict(size=20)),
                
                dict(x=0.8, y=1.25,
                buttons=[dict(label='Horizontal', method='restyle', args=['orientation', 'h']),
                        dict(label='Vertical',method='restyle',args=['orientation', 'v'])],font=dict(size=20))])
        
    def build(self):
        """Build the co-expression network"""
        fig = self._create_subplot()
        self._create_sankey(fig)
        self._create_ternary(fig)
        self._add_functionality(fig)
        return fig