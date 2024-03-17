import sys
from pathlib import Path

root_path = Path(__file__).parent.parent.parent
sys.path.insert(0, str(root_path.absolute()))     # Sphinxがライブラリを探索できるようにパスを追加


# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'my_sphinx'
copyright = '2024, naoki'
author = 'naoki'
release = '2024'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon','sphinx.ext.mathjax', 'sphinx.ext.todo', 'sphinx.ext.githubpages']

templates_path = ['_templates']
exclude_patterns = []

language = 'ja'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_title = "    "
html_static_path = ['_static']
html_logo = "https://github.com/naonao07/my_sphinx/_static/logo.png"
#html_logo = "logo.png"


