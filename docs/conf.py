project = "CLPipe"
author = "LSST DESC CLPipe Contributors"
release = "1.0.0"

import os
import shutil
import sys

sys.path.insert(0, os.path.abspath(".."))

# The tutorials live in ../tutorials. Copy their pages next to the docs so that
# Sphinx renders them, without keeping a second version in git.
_here = os.path.dirname(os.path.abspath(__file__))
shutil.rmtree(os.path.join(_here, "tutorials"), ignore_errors=True)
shutil.copytree(
    os.path.join(_here, "..", "tutorials"),
    os.path.join(_here, "tutorials"),
    ignore=shutil.ignore_patterns("data", "outputs", ".ipynb_checkpoints", "*.log"),
)

root_doc = "index"
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
    "sphinx_copybutton",
    "myst_nb",
]

# The notebooks are stored with their outputs, so they are not run here.
nb_execution_mode = "off"
myst_enable_extensions = ["dollarmath"]

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
