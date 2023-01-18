import os
import sys

sys.path.insert(0, os.path.abspath(".."))  # NOQA
import audioflux

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'AudioFlux'
copyright = '2022, immusician'
author = 'immusician'
version = release = audioflux.__version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "sphinx.ext.autosectionlabel",
    "numpydoc",
    "matplotlib.sphinxext.plot_directive",
]

plot_include_source = True
plot_html_show_formats = False
plot_html_show_source_link = False
plot_formats = [("png", 100)]

numpydoc_use_plots = True
numpydoc_show_class_members = True
numpydoc_class_members_toctree = False

# mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"
mathjax_path = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.3/MathJax.js?config=TeX-AMS-MML_HTMLorMML"

default_role = "autolink"

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

language = 'en'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'
html_favicon = '../image/icon.png'
html_static_path = ['_static']
html_css_files = [
    'css/custom.css',
]

autodoc_member_order = 'bysource'

autosectionlabel_prefix_document = True
