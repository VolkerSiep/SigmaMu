# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

# -- Project information -----------------------------------------------------
import sys
from os.path import abspath
from simu import __version__ as release

sys.path.insert(0, abspath('.'))

project = 'SiMu'
copyright = '2021-2024, Volker Siepmann'
author = 'Volker Siepmann'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.todo',
    'sphinx.ext.autosummary',
    'sphinx_copybutton',
    'sphinxcontrib.bibtex',
    'sphinx.ext.doctest',
    'custom_directives'
]

# nitpicky = True  # to check that all internal references work

# autoclass_content = 'both'
autodoc_member_order = 'bysource'
autodoc_typehints = 'signature'
# autodoc_imported_members = True

# copy-button config
copybutton_prompt_text = r'>>> |\.\.\. '
copybutton_prompt_is_regexp = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns: list[str] = []

# bibtex file
bibtex_bibfiles = ['bibliographies.bib']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_theme = 'nature'

# other ok themes:
# 'sphinx_rtd_theme' (search hangs),
# 'pyramid' (headings not standing out),
# 'classic' (a bit boring)

# manual configuration of side-bars to have the global toc folded and not
# messed up with all autodoc entries
html_sidebars = {
    '**': ['globaltoc.html', 'sourcelink.html', 'searchbox.html'],
    'using/windows': ['windowssidebar.html', 'searchbox.html'],
}

html_css_files = ["custom.css"]

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_show_sourcelink = False
html_favicon = "_static/goose.png"
todo_include_todos = True

mathjax3_config = {
    "tex": {
        "macros": {
            "emul": r'{\stackrel{E}{\cdot}}',
            "vec": [r'{\mathbf{\boldsymbol{#1}}}', 1],
            "standard": r'{\circ\hspace{-1.45ex}-}'
        }
    }
}

rst_prolog = r"""
.. |m3| replace:: m\ :sup:`3`
.. |degC| replace:: °C
"""
