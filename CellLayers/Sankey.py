from .MultiResolutionAnalysis import *

class Sankey:
    """
    Class for building a Sankey for multi-resolution single cell analysis  

    Parameters
    ----------
    sankey_dict: dict (required)
        Dictionary containing the multi-resolution cluster analysis for building a Sankey network
    node_color: str (optional)
        String to define node color

    Attributes
    ----------
    sankey_dict: dict
        Dictionary containing the multi-resolution cluster analysis for building a Sankey network
    """
    #TODO: Update node_color accessibility to be in self.sankey_dict
    #TODO: Create a method in MultiResolutionAnalysis for coloring nodes
    #TODO: Create a method for accessing additional meta attributes

    def __init__(self, sankey_dict, node_color='#F7ED32'):
        self.sankey_dict = sankey_dict
        self._starter_gene = self.sankey_dict['genes'][0]
        self._node_color = node_color

    def _create_sankey(self) -> 'plotly.graph_objs._figure.Figure':
        """Create a Sankey Plotly graph object"""
        fig = go.Figure(data=[dict(type='sankey', orientation='h', 
                node = dict(pad = 10, thickness=10,
                            label=self.sankey_dict['node_data']['node_labels'],
                            # customdata = self.sankey_dict['node_data']['new_label'],
                            hovertemplate= '%{customdata}',
                            color = self.sankey_dict['node_data']['silhoutte_hex']),
                
                link = dict(source = self.sankey_dict['data']['source'],
                            target = self.sankey_dict['data']['target'],
                            color =  self.sankey_dict['data'][self._starter_gene+'_hex'],
                            value = self.sankey_dict['data']['value']),textfont=dict(color='black',size=15))])
        fig.update_xaxes(showticklabels=False) 
        fig.update_yaxes(showticklabels=False)
        fig['layout']['showlegend'] = False
        fig['layout']['xaxis']['showgrid'] = False
        fig['layout']['yaxis']['showgrid'] = False
        fig.update_layout(xaxis_zeroline=False, yaxis_zeroline=False)
        fig.update_layout({'plot_bgcolor': 'rgba(0, 0, 0, 0)','paper_bgcolor': 'rgba(0, 0, 0, 0)'})
        return fig
    
    def _create_gene_expression(self, fig) -> list:
        """
        Add gene expression features to the Sankey
        
        Parameters
        ----------
        fig: plotly.graph_objs._figure.Figure
            Plotly Figure subclass from the graph_objects class
        """
        gene_buttons = []
        fig.add_trace(go.Scatter(x=[None],
                                 y=[None],
                                 mode='markers',
                                 visible=True,
                                 marker={'colorscale':self.sankey_dict['node_data']['silhoutte_hex'].tolist()}))
        for i in range(len(self.sankey_dict['genes'])):
            if i == 0:
                self._create_expression_colorbar(fig, self.sankey_dict['genes'][i], True)
            else:
                self._create_expression_colorbar(fig, self.sankey_dict['genes'][i], False)
            gene_bools = np.full(len(self.sankey_dict['genes']), False).tolist()
            gene_bools[i] = True
            gene_bools[0:0] = [True]
            gene_bools.insert(1, True)
            gene_dict = dict(
                label=self.sankey_dict['genes'][i], 
                method='update', 
                args=[{"visible":gene_bools,
                'link' : dict(source = self.sankey_dict['data']['source'],
                target = self.sankey_dict['data']['target'],
                color = self.sankey_dict['data'][self.sankey_dict['genes'][i]+'_hex'], # Gene expression or largest dag trace parameter
                value = self.sankey_dict['data']['value'])}])
            gene_buttons.append(gene_dict)
            
        return gene_buttons
    
    def _create_expression_colorbar(self, fig, gene, visible=False):
        """"
        Create an expression color bar for a gene

        Parameters
        ----------
        fig: plotly.graph_objs._figure.Figure
            Plotly Figure subclass from the graph_objects class
        gene: str
            String that is a gene name
        visible: boolean
            Boolean argument for showing expression color bar
        """
        fig.add_trace(go.Scatter(x=[None],
                                 y=[None],
                                 mode='markers',
                                 visible=visible,
                                 marker=self.sankey_dict['exp_colorbar'][gene]))
    
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
    
    def _add_functionality(self, fig):
        """
        Add drop down menus for selecting genes and configuring the layout

        Parameters
        ----------
        fig: plotly.graph_objs._figure.Figure
            Plotly Figure subclass from the graph_objects class
        """
        gene_buttons = self._create_gene_expression(fig)
        fig.update_layout(
            updatemenus=[
                dict(y=0.9, buttons=list(gene_buttons), font=dict(size=15)),
                
                dict(y=0.50, buttons=[dict(label='Silhouette Scores', 
                                           method='update', args=[{"visible":True}])],font=dict(size=15)),
                
                dict(x=0, y=1.30,
                buttons=[dict(label='Snap',method='restyle', args=['arrangement', 'snap']),
                        dict(label='Perpendicular', method='restyle',args=['arrangement', 'perpendicular']),
                        dict(label='Freeform', method='restyle',args=['arrangement', 'freeform']),
                        dict(label='Fixed', method='restyle',args=['arrangement', 'fixed'])],font=dict(size=15)),    
                
                dict(x=0.2, y=1.30,
                buttons=[dict(label='Light', method='relayout', args=['paper_bgcolor', 'white']),
                        dict(label='Dark', method='relayout', args=['paper_bgcolor', 'black'])],font=dict(size=15)),
                
                dict(x=0.4, y=1.30,
                buttons=[dict(label='Thin', method='restyle',args=['node.thickness', 8]),
                        dict(label='Thick',method='restyle',args=['node.thickness', 15])],font=dict(size=15)),
                
                dict(x=0.6, y=1.30,
                buttons=[dict(label='Small gap',method='restyle',args=['node.pad', 15]),
                        dict(label='Large gap',method='restyle',args=['node.pad', 20])],font=dict(size=15)),
                
                dict(x=0.8, y=1.30,
                buttons=[dict(label='Horizontal', method='restyle', args=['orientation', 'h']),
                        dict(label='Vertical',method='restyle',args=['orientation', 'v'])],font=dict(size=15))])
                        
        self._create_silhouette_colorbar(fig)

        
    def build(self):
        """Build the Sankey network"""
        fig = self._create_sankey()
        self._add_functionality(fig)       
        return fig