# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# autodoc imports the *installed* ttcrpy: rgrid and tmesh are compiled
# extensions, so the in-tree src/ directory holds only .pyx sources and must not
# be put on sys.path -- doing so shadows the installed package and every
# automodule directive then silently produces an empty page.
import sphinx_rtd_theme
from importlib.metadata import version as get_version


# -- Project information -----------------------------------------------------

project = 'ttcrpy'
copyright = '2026, Bernard Giroux'
author = 'Bernard Giroux'

# The full version, including alpha/beta/rc tags
release = get_version('ttcrpy')

# The short X.Y version.  The epub3 builder refuses an empty `version`, and
# .readthedocs.yaml builds epub with fail_on_warning, so leaving this unset
# fails the build even though the HTML is fine.  Taking the first two
# components keeps the dev suffix setuptools_scm appends out of it:
# "1.4.3.dev50+g533efaf14" -> "1.4".
version = '.'.join(release.split('.')[:2])


# -- General configuration ---------------------------------------------------

master_doc = 'index'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx_rtd_theme'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Napoleon settings
#
# Write array shapes in a parameter's type as "np.ndarray with shape (n, 3)",
# never "of shape" or ", shape".  Sphinx splits a type field on the delimiters
# in PyXrefMixin._delimiters_re -- brackets, commas, "|", and (since Sphinx 9)
# both "or" and "of" -- and cross-references every token in between.  The other
# two spellings therefore leave a bare "shape" token, which resolves against
# the four Grid*.shape properties: ambiguous when referenced from tmesh, which
# is an error under the fail_on_warning setting in .readthedocs.yaml, and
# silently wrong from rgrid, where it linked an array shape to the property.
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = False
napoleon_use_ivar = True
napoleon_use_param = True
napoleon_use_rtype = True
