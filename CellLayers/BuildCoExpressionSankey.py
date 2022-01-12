import numpy as np

from plotly.subplots import make_subplots
import plotly.graph_objects as go

import matplotlib
import matplotlib.colors as mcolors


class BuildCoExpressionSankey:
    
    def __init__(self, sankey_dict):
        self.sankey_dict = sankey_dict
        self.starter_tri_coexpressed_genes = list(sankey_dict['tri_coexp_color'].keys())[0]
    
    def _make_axis(self, title, tickangle):
        return {
            'title': title,
            'titlefont': { 'size': 20 },
            'tickangle': tickangle,
            'tickfont': { 'size': 15 },
            'tickcolor': 'rgba(0,0,0,0)',
            'ticklen': 5,
            'showline': True,
            'showgrid': True}
    
    def _create_silhouette_color_bar(self, fig):
        fig.add_trace(go.Scatter(x=[None],
                   y=[None],
                   mode='markers',
                   visible=True,
                   marker={'colorscale':[mcolors.to_hex(self.sankey_dict['silhouette_mapper'].to_rgba(x)) 
                                         for x in self.sankey_dict['silhouette_list']],
                           'showscale':True, 'cmin':-1, 'cmax':1, 'colorbar': {'x':1.1}}))
    
    def _create_subplot_obj(self):
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
    
    def _create_sankey(self, fig):
        fig.add_trace(go.Sankey(
            node=go.sankey.Node(label= self.sankey_dict['node_labels']),
            link=go.sankey.Link(source = self.sankey_dict['data']['source'],
                                target = self.sankey_dict['data']['target'],
                                color = [matplotlib.colors.to_hex(x) for x in self.sankey_dict['tri_coexp_color'][self.starter_tri_coexpressed_genes]],
                                value = self.sankey_dict['data']['value']),),
                      row=1, 
                      col=2)
    
    def _create_ternary(self, fig):
        fig.add_trace(go.Scatterternary({
            'mode':'markers',
            'a':np.array(self.sankey_dict['tri_coexp_dict'][self.starter_tri_coexpressed_genes])[:,0],
            'b':np.array(self.sankey_dict['tri_coexp_dict'][self.starter_tri_coexpressed_genes])[:,1],
            'c':np.array(self.sankey_dict['tri_coexp_dict'][self.starter_tri_coexpressed_genes])[:,2],
            'text': 'Percentile Expression',
            'hovertemplate':'MS4A1: %{a} <br> FCER1A: %{b} <br> FCGR3A: %{c}',
            'marker':{'color':self.sankey_dict['tri_coexp_color'][self.starter_tri_coexpressed_genes],
                      'size': np.log2(self.sankey_dict['data']['value'].tolist())}}),
                      row=1, 
                      col=1)
        
    def run(self):
        fig = self._create_subplot_obj()
        self._create_silhouette_color_bar(fig)
        self._create_sankey(fig)
        self._create_ternary(fig)
        fig.update_layout(height=700, showlegend=False)
        fig.update_layout({'ternary':{
            'sum': 100,
            'aaxis': self._make_axis(self.starter_tri_coexpressed_genes[0], 0),
            'baxis': self._make_axis('<br>'+self.starter_tri_coexpressed_genes[1], 45),
            'caxis': self._make_axis('<br>'+self.starter_tri_coexpressed_genes[2], -45)}})
        return fig