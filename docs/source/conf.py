# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os,sys

sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('../../pide'))

project = 'pide'
copyright = '2024, Sinan Ozaydin'
author = 'Sinan Ozaydin'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'myst_parser', 
    'sphinx.ext.mathjax',
    'nbsphinx',           # For rendering Jupyter Notebooks
]

myst_enable_extensions = [
    "deflist",           # Optional: Add definition lists support
    "linkify",           # Optional: Auto-convert URLs to links
]

# You may also want to customize the MathJax options
mathjax_config = {
    'tex': {
        'inlineMath': [['$', '$'], ['\\(', '\\)']],  # Inline math delimiters
        'displayMath': [['$$', '$$'], ['\\[', '\\]']],  # Block math delimiters
    },
}

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

nbsphinx_execute = 'never'  # Disables automatic execution of notebooks
