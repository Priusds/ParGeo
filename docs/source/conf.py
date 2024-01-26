# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import sys, os
sys.path.insert(0, os.path.abspath('../bubbles'))

project = 'Bubbles'
copyright = '2024, Till'
author = 'Till'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
extensions = [
  'sphinx.ext.autodoc',
  'sphinx.ext.autosummary',
  'sphinx.ext.doctest',
  'sphinx.ext.duration',
  'sphinx.ext.napoleon',
]




# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']
