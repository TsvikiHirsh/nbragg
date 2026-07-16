# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'nbragg'
copyright = '2026, Tsviki Y. Hirsh'
author = 'Tsviki Y. Hirsh'

try:
    from importlib.metadata import version as _pkg_version
    release = _pkg_version('nbragg')
except Exception:
    release = '1.0.0'
version = '.'.join(release.split('.')[:2])

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx_autodoc_typehints',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',  # Support for math equations
]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
templates_path = ['_templates']
pygments_style = 'sphinx'
autosummary_generate = True  # Generate autosummary pages

# Add your project source directory to sys.path
import os
import sys
sys.path.insert(0, os.path.abspath('../..'))
