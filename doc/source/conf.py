# more info on configuring this file:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys
from os.path import abspath
from simu import __version__ as _release

# -- Path setup --------------------------------------------------------------
sys.path.insert(0, abspath('.'))

# -- Project information -----------------------------------------------------
project = 'SigmaMu'
copyright = '2021-2025, Volker Siepmann'
author = 'Volker Siepmann'
release = _release


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx_licenseinfo',
    'sphinx.ext.todo',
    'sphinx.ext.autosummary',
    'sphinx_copybutton',
    'sphinxcontrib.bibtex',
    'sphinx.ext.doctest',
    'custom_directives',
    'sphinx.ext.intersphinx',
    # 'sphinx_autodoc_typehints'
]

nitpicky = True  # to check that all internal references work


intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
}

nitpick_ignore = [
    ("py:class", "simu.core.utilities.types.__V"),
    ("py:data", "simu.core.utilities.structures._V"),
    ("py:data", "simu.core.utilities.qstructures._V"),
    ("py:class", "casadi.SX"),
    ("py:class", "casadi.casadi.SX"),  # no clue why it tries to find this!?
]


# autoclass_content = 'both'
autodoc_member_order = 'bysource'
autodoc_class_signature = 'separated'
autodoc_typehints = 'signature'

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

# break signature into multiple lines when length is greater as given value
maximum_signature_line_length = 90  #


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_theme = 'nature'

# other ok themes:
# 'sphinx_rtd_theme' (search hangs),
# 'pyramid' (headings not standing out),
# 'classic' (a bit boring)

# manual configuration of sidebars to have the global toc folded and not
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
rst_epilog = r"""
.. _pytest: https://docs.pytest.org/
.. _CasADi: https://web.casadi.org
.. _NumPy: https://numpy.org/
.. _SciPy: https://scipy.org/
.. _Pint: https://pint.readthedocs.io
.. _PyYAML: https://pyyaml.org/
.. _PyPardiso: https://pypi.org/project/pypardiso/
.. _matplotlib: https://matplotlib.org/
.. _Sphinx: https://www.sphinx-doc.org
.. _sphinxcontrib-bibtex: https://github.com/mcmtroffaes/sphinxcontrib-bibtex
.. _sphinx-licenseinfo: https://sphinx-licenseinfo.readthedocs.io/
.. _sphinx_copybutton: https://github.com/executablebooks/sphinx-copybutton
.. _pytest-doctestplus: https://github.com/scientific-python/pytest-doctestplus
.. _Matplotlib: https://matplotlib.org/
"""