# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'CellLayers'
# copyright = '2021, Graziella'
author = 'Andrew Blair'

release = '0.1'
version = '0.1.0'

# -- General configuration

extensions =  ['nbsphinx',
    'sphinx.ext.mathjax',
    'sphinx_copybutton',
    'sphinx.ext.githubpages']
# extensions = [
#     'sphinx.ext.duration',
#     'sphinx.ext.doctest',
#     'sphinx.ext.autodoc',
#     'sphinx.ext.autosummary',
#     'sphinx.ext.intersphinx',
#     'sphinx.ext.mathjax',
#     'sphinx.ext.githubpages',
#     'nbsphinx'
# ]

exclude_patterns = ['**.ipynb_checkpoints']

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'
html_scaled_image_link = True #False
html_static_path = ['_static']
# html_style = '_static/style.css'



import sphinx_rtd_theme
html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
html_style = 'css/custom.css'

html_context = { 
'css_files': [
    'https://media.readthedocs.org/css/sphinx_rtd_theme.css',
    'https://media.readthedocs.org/css/readthedocs-doc-embed.css',
    '_static/css/style.css',
],  
}   

  
# -- Options for EPUB output
epub_show_urls = 'footnote'

nbsphinx_execute_arguments = [
    "--InlineBackend.figure_formats={'svg', 'pdf'}",
    "--InlineBackend.rc={'figure.dpi': 200}",
]



import plotly.io as pio
pio.renderers.default = 'sphinx_gallery'