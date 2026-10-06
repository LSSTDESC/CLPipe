project = "CLPipe"
author = "LSST DESC CLPipe Contributors"
release = "0.0.1"

import os
import sys

sys.path.insert(0, os.path.abspath(".."))

root_doc = "index"
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
    "sphinx_copybutton",
]

autosummary_generate = True
autodoc_member_order = "bysource"
autodoc_default_options = {"members": True, "show-inheritance": True}
# Heavy DESC dependencies are not installed in the standalone docs workflow.
autodoc_mock_imports = [
    "ceci",
    "crow",
    "firecrown",
    "pyccl",
    "sacc",
    "tjpcov",
    "cosmosis",
    "yaml",
    "numpy",
]
exclude_patterns = [
	"_build",
	"venv/**",
	"Thumbs.db",
	".DS_Store",
]

html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "prev_next_buttons_location": None,
    "collapse_navigation": False,
    "titles_only": True,
}
