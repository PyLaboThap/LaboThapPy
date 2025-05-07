import os
import sys

# Add the path to the library directory
sys.path.insert(0, os.path.abspath('../../labothappy'))




# -- Project information -----------------------------------------------------
project = 'LaboThapPy'
copyright = '2024, Basile Chaudoir, Elise Neven'
author = 'Elise Neven, Basile Chaudoir'
release = '1.0'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',   # To generate documentation from docstrings
    'sphinx.ext.napoleon',  # For Google/NumPy style docstrings
    'sphinx.ext.viewcode',  # Add links to highlighted source code
    'sphinx.ext.todo',      # Support for TODO notes
    'sphinx.ext.mathjax',  # For mathematical equations
    'nbsphinx',  # For including Jupyter Notebooks
]

# Configure nbsphinx
nbsphinx_execute = 'auto'  # Options are 'auto', 'always', 'never'
nbsphinx_allow_errors = True

autodoc_member_order = 'bysource'
autodoc_typehints = 'description'

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**/.ipynb_checkpoints']

# -- Options for HTML output -------------------------------------------------
# The theme to use for HTML and HTML Help pages.
html_theme = 'furo'
html_static_path = ['_static']

def setup(app):
    app.add_css_file('custom.css')

html_title = "LaboThapPy Documentation"

# -- Options for mathematical equations ------------------------------------------------

# MathJax settings
mathjax_path = 'https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js'



