# Build the ttcr documentation.
#
# Two independent doc sets, built by different tools into separate trees:
#
#   C++     Doxygen, from the headers reachable from the three command-line
#           programs (ttcr3d, ttcr2d, ttcr2ds)   ->  ./html
#   Python  Sphinx, from the ttcrpy docstrings   ->  ../_build/html
#
#   make              # C++ only (default)
#   make cpp          # C++ only
#   make python       # Python only
#   make both         # both
#   make open         # build C++ and open it (macOS)
#   make open-python  # build Python and open it (macOS)
#   make input        # regenerate the Doxyfile INPUT list from the compiler
#   make check-input  # report drift in that list without rewriting it
#   make clean / clean-python / clean-all
#   make help

DOXYGEN ?= doxygen

# Sphinx lives one level up and owns its own Makefile and _build tree.
SPHINXDIR ?= ..

.PHONY: all cpp python both open open-python input check-input \
        clean clean-python clean-all help

all: cpp

# ---------------------------------------------------------------------------
# C++ (Doxygen)
# ---------------------------------------------------------------------------

cpp: html/index.html

# The prerequisite list is every header, not just those in INPUT: editing any of
# them should trigger a rebuild, and a header newly pulled into the include
# closure should not go unnoticed.
html/index.html: Doxyfile $(wildcard ../../ttcr/*.h)
	$(DOXYGEN) Doxyfile

open: cpp
	open html/index.html

# ---------------------------------------------------------------------------
# Python (Sphinx)
# ---------------------------------------------------------------------------
# Delegated to the Sphinx Makefile in the parent directory, so exactly one place
# knows how to build the Python docs. Requires ttcrpy to be importable
# (pip install -e .) — autodoc imports the package rather than parsing it.

python:
	$(MAKE) -C $(SPHINXDIR) html

open-python: python
	open $(SPHINXDIR)/_build/html/index.html

both: cpp python

# ---------------------------------------------------------------------------
# Maintenance
# ---------------------------------------------------------------------------
# Recompute the Doxyfile INPUT list — the transitive include closure of the
# three programs, as reported by the compiler. Override these when the third-
# party include paths differ:
#
#   make input INCFLAGS="-I../../ttcr -I/opt/local/include/eigen3"
#   make input VTKFLAGS="-DVTK -I/opt/local/include/vtk-9.6"

CXX      ?= c++
SRCDIR   ?= ../../ttcr
PROGRAMS ?= ttcr3d.cpp ttcr2d.cpp ttcr2ds.cpp
INCFLAGS ?= -I$(SRCDIR) -I../../eigen-5.0.0 -I../../boost_1_91_0
VTKFLAGS ?= -DVTK

# The --opt=value form is required, not --opt value: these values begin with a
# dash (-DVTK, -I...), which argparse would otherwise read as another option.
input:
	@python3 update_input.py --cxx="$(CXX)" --srcdir="$(SRCDIR)" \
	    --programs="$(PROGRAMS)" --incflags="$(INCFLAGS)" --vtkflags="$(VTKFLAGS)"

check-input:
	@python3 update_input.py --check --cxx="$(CXX)" --srcdir="$(SRCDIR)" \
	    --programs="$(PROGRAMS)" --incflags="$(INCFLAGS)" --vtkflags="$(VTKFLAGS)"

# ---------------------------------------------------------------------------
# Cleaning
# ---------------------------------------------------------------------------

clean:
	rm -rf html doxygen-warnings.log

clean-python:
	$(MAKE) -C $(SPHINXDIR) clean

clean-all: clean clean-python

help:
	@echo "C++ (Doxygen)   -> ./html"
	@echo "  make cpp           build"
	@echo "  make open          build and open"
	@echo "  make clean         remove output"
	@echo ""
	@echo "Python (Sphinx) -> $(SPHINXDIR)/_build/html"
	@echo "  make python        build"
	@echo "  make open-python   build and open"
	@echo "  make clean-python  remove output"
	@echo ""
	@echo "  make both          build both"
	@echo "  make clean-all     remove both"
	@echo "  make input         regenerate the Doxyfile INPUT list"
	@echo "  make check-input   report drift in that list"
