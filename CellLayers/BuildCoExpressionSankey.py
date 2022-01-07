import numpy as np

from plotly.subplots import make_subplots
import plotly.graph_objects as go

import matplotlib
import matplotlib.colors as mcolors


import numpy as np

from plotly.subplots import make_subplots
import plotly.graph_objects as go

import matplotlib
import matplotlib.colors as mcolors

def makeAxis(title, tickangle):
    return {
        'title': title,
        'titlefont': { 'size': 20 },
        'tickangle': tickangle,
        'tickfont': { 'size': 15 },
        'tickcolor': 'rgba(0,0,0,0)',
        'ticklen': 5,
        'showline': True,
        'showgrid': True}

fig = make_subplots(rows=1, cols=2, specs=[[{"type": "ternary"},{"type": "sankey"}]])

fig.add_trace(go.Scatter(x=[None],
                   y=[None],
                   mode='markers',
                   visible=True,
                   marker={'colorscale':[mcolors.to_hex(sankey_dict['silhouette_mapper'].to_rgba(x)) 
                                         for x in sankey_dict['silhouette_list']],
                           'showscale':True, 'cmin':-1, 'cmax':1, 'colorbar': {'x':1.1}}))


fig.update_xaxes(showticklabels=False) # hide all the xticks
fig.update_yaxes(showticklabels=False) # hide all the xticks
fig['layout']['showlegend'] = False
fig['layout']['xaxis']['showgrid'] = False
fig['layout']['yaxis']['showgrid'] = False
fig.update_layout(xaxis_zeroline=False, yaxis_zeroline=False)
fig.update_layout({
'plot_bgcolor': 'rgba(0, 0, 0, 0)',
'paper_bgcolor': 'rgba(0, 0, 0, 0)'
})


fig.add_trace(go.Sankey(
    node=go.sankey.Node(label=sankey_dict['node_labels']),
    link=go.sankey.Link(source = sankey_dict['data']['source'],
                        target = sankey_dict['data']['target'],
                        color = [matplotlib.colors.to_hex(x) for x in sankey_dict['tri_coexp_color'][('MS4A1', 'FCER1A', 'FCGR3A')]],
                        value = sankey_dict['data']['value']
                       ),),row=1, col=2)

fig.add_trace(go.Scatterternary({
    'mode':'markers',
    'a':np.array(sankey_dict['tri_coexp_dict'][('MS4A1', 'FCER1A', 'FCGR3A')])[:,0],
    'b':np.array(sankey_dict['tri_coexp_dict'][('MS4A1', 'FCER1A', 'FCGR3A')])[:,1],
    'c':np.array(sankey_dict['tri_coexp_dict'][('MS4A1', 'FCER1A', 'FCGR3A')])[:,2],
    'text': 'Percentile Expression',
    'hovertemplate':'MS4A1: %{a} <br> FCER1A: %{b} <br> FCGR3A: %{c}',
    'marker':{'color':sankey_dict['tri_coexp_color'][tuple(['MS4A1', 'FCER1A', 'FCGR3A'])],
              'size': np.log2(sankey_dict['data']['value'].tolist())}}),
              row=1, col=1)

fig.update_layout(height=700, showlegend=False)

fig.update_layout({'ternary':{
    'sum': 100,
    'aaxis': makeAxis('MS4A1', 0),
    'baxis': makeAxis('<br>FCER1A', 45),
    'caxis': makeAxis('<br>FCGR3A', -45)}})


fig.show()