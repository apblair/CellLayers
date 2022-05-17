from setuptools import setup
required = ["matplotlib==3.5.1", "numpy==1.21.2", "pandas==1.4.1", "plotly==5.7.0", "chart-studio==1.1.0", "pandas==1.4.1", "scipy==1.7.3", "pytest"]

setup(name='CellLayers',
      description='Cell Layers is an interactive Sankey tool for the quantitative investigation of gene expression, coexpression, biological processes, and cluster integrity across clustering resolutions.',
      url='https://github.com/apblair/CellLayers',
      author='Andrew Blair',
      author_email='andrew.blair@ucsf.edu',
      license='MIT',
      packages=['CellLayers'],
      install_requires=required,
      packages=setup.find_packages(where="src"),
      zip_safe=False)