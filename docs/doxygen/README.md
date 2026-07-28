# ttcr documentation

Two independent documentation sets, built by different tools into separate
output trees. The Makefile here drives both, but nothing is merged: the C++
reference and the Python reference stay entirely separate.

| | tool | source | output |
|---|---|---|---|
| **C++** | Doxygen | `../../ttcr/*.h` | `docs/doxygen/html` |
| **Python** | Sphinx | `ttcrpy` docstrings | `docs/_build/html` |

## Requirements

- **C++** — only [Doxygen](https://www.doxygen.nl/). This configuration needs
  neither Graphviz/`dot` nor LaTeX.
- **Python** — Sphinx with `sphinx_rtd_theme`, and `ttcrpy` importable
  (`pip install -e .` from the repository root). Sphinx `autodoc` imports the
  package rather than parsing it, so an unbuilt or uninstalled `ttcrpy` makes
  the Python build fail.

```bash
# macOS
brew install doxygen
# Debian/Ubuntu
sudo apt-get install doxygen
```

## Building

```bash
cd docs/doxygen

make              # C++ only (the default)
make cpp          # C++ only
make python       # Python only
make both         # both

make open         # build the C++ docs and open them (macOS)
make open-python  # build the Python docs and open them (macOS)

make clean        # remove the C++ output
make clean-python # remove the Python output
make clean-all    # remove both

make help         # this list
```

Doxygen problems are collected in `doxygen-warnings.log`; the build is currently
warning-free.

## Scope of the C++ reference

`INPUT` in the `Doxyfile` is **not** the whole `ttcr/` directory. It is an
explicit list of the 77 headers reachable from the three command-line programs
— `ttcr3d.cpp`, `ttcr2d.cpp` and `ttcr2ds.cpp` — that is, their transitive
include closure.

Three headers in `ttcr/` are therefore absent, none of them reachable from those
programs:

| file | why |
|---|---|
| `Interface.h` | unused throughout the project |
| `msh2vtk_io.h` | belongs to the separate `msh2vtk` utility |
| `structs_msh2vtk.h` | likewise |

### Keeping the list honest

Do not edit `INPUT` by hand. Regenerate it — the compiler is asked for each
program's include closure (`-MM`) and the block is rewritten in place:

```bash
make input        # rewrite the INPUT list
make check-input  # report drift, change nothing, exit non-zero if stale
```

`make check-input` is the one to run in CI: it fails if a header has entered or
left the closure without the `Doxyfile` being updated.

Both targets take the usual overrides when the third-party include paths differ:

```bash
make input INCFLAGS="-I../../ttcr -I/opt/local/include/eigen3 -I/opt/local/include"
make input VTKFLAGS="-DVTK -I/opt/local/include/vtk-9.6"
```

`VTKFLAGS` defaults to `-DVTK` with no extra include path. VTK matters here
because `VTUReader.h` is reachable *only* from the VTK code paths — probe
without it and that header silently drops out of the closure. If the VTK probe
fails, `update_input.py` says so on stderr and falls back to the non-VTK
closure, so a stale list is reported rather than quietly produced.

## Notes

- `EXTRACT_ALL` is on, so the reference stays complete even for headers that are
  not yet fully annotated.
- Optional feature sections gated behind `VTK` / `_OPENMP` are omitted from the
  rendered output unless those macros are added to `PREDEFINED` in the
  `Doxyfile`. This is independent of `VTKFLAGS` above, which only affects which
  files are *listed* in `INPUT`.
- `ctpl_stl.h` (vendored thread pool) is part of the closure and so is
  documented. `nanoflann.hpp` is also a dependency but is excluded by
  `FILE_PATTERNS = *.h`.
- The Python build currently emits some pre-existing "duplicate object
  description" warnings, from the `Grid2d_d`/`Grid3d_d` dtype variants
  documenting identically-named attributes.
