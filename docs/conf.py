# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/stable/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.insert(0, os.path.abspath(os.pardir))


# -- Project information -----------------------------------------------------

project = "OpenFF Fragmenter"
copyright = "2018, Chaya D. Stern"
author = "Chaya D. Stern"

# The short X.Y version
version = ""
# The full version, including alpha/beta/rc tags
release = ""


# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "nbsphinx",
    "nbsphinx_link",
    "sphinxcontrib.bibtex",
    "sphinxcontrib.autodoc_pydantic",
    "openff_sphinx_theme",
]

source_suffix = ".rst"

master_doc = "index"

language = None

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "default"

# Autodoc settings
autosummary_generate = True

autodoc_default_options = {
    "member-order": "bysource",
}

autodoc_mock_imports = [
    "openff.toolkit",
]

# Napoleon settings
napoleon_numpy_docstring = True
napoleon_use_rtype = False

# autodoc_pydantic settings
autodoc_pydantic_show_config = False
autodoc_pydantic_model_show_config = False
autodoc_pydantic_show_validators = False
autodoc_pydantic_model_show_validators = False

# nbsphinx settings
nbsphinx_execute = "never"

# sphinx bibtext settings
bibtex_bibfiles = ["index.bib"]

# Set up the intershinx mappings.
intersphinx_mapping = {
    "python": ("https://docs.python.org/", None),
    "numpy": ("https://docs.scipy.org/doc/numpy/", None),
    "openff.toolkit": (
        "https://open-forcefield-toolkit.readthedocs.io/en/latest/",
        None,
    ),
}

# Set up mathjax.
mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js"

# -- Options for HTML output -------------------------------------------------

html_theme = "openff_sphinx_theme"
html_sidebars = {"**": ["globaltoc.html", "localtoc.html", "searchbox.html"]}

html_theme_options = {
    # Repository integration
    # Set the repo url for the link to appear
    "repo_url": "https://github.com/openforcefield/openff-fragmenter",
    # The name of the repo. If must be set if repo_url is set
    "repo_name": "openff-fragmenter",
    # Must be one of github, gitlab or bitbucket
    "repo_type": "github",
    # Colour for sidebar captions and other accents. One of
    # openff-blue, openff-toolkit-blue, openff-dataset-yellow,
    # openff-evaluator-orange, aquamarine, lilac, amaranth, grape,
    # violet, pink, pale-green, green, crimson, eggplant, turquoise,
    # or a tuple of three ints in the range [0, 255] corresponding to
    # a position in RGB space.
    "color_accent": "aquamarine",
}

html_static_path = ["_static"]

if importlib.util.find_spec("notfound"):
    extensions.append("notfound.extension")

    notfound_urls_prefix = (
        "https://docs.openforcefield.org/projects/fragmenter/en/stable/"
    )
    notfound_context = {
        "title": "404: File Not Found",
        "body": f"""
    <h1>404: File Not Found</h1>
    <p>
        Sorry, we couldn't find that page. This often happens as a result of
        following an outdated link. Please check the
        <a href="{notfound_urls_prefix}">latest stable version</a>
        of the docs, unless you're sure you want an earlier version, and
        try using the search box or the navigation menu on the left.
    </p>
    <p>
    </p>
    """,
    }

# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = "fragmenterdoc"

# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    "papersize": "letterpaper",
    "pointsize": "10pt",
    "preamble": "",
    "figure_align": "htbp",
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, "fragmenter.tex", "OpenFF Fragmenter Documentation", author, "manual"),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, "openff-fragmenter", "OpenFF Fragmenter Documentation", [author], 1)
]

# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        master_doc,
        "openff-fragmenter",
        "OpenFF Fragmenter Documentation",
        author,
        "openff-fragmenter",
        "Fragment molecules for quantum mechanics torsion scans.",
        "Miscellaneous",
    ),
]
