#!/usr/bin/env python3
"""Regenerate the Doxyfile INPUT list.

The C++ reference is scoped to the headers actually reachable from the three
command-line programs (ttcr3d, ttcr2d, ttcr2ds). Rather than curating that list
by hand, this asks the compiler for each program's transitive include closure
(``-MM``) and rewrites the INPUT block in place.

Run it through the Makefile::

    make input
    make input CXX=g++ INCFLAGS="-I../../ttcr -I/opt/include"

VTK is enabled when probing, so VTUReader.h — reachable only from the VTK code
paths — is not silently dropped. If the VTK probe fails the program is retried
without it, and any header that only VTK pulls in is reported as missing.
"""

import argparse
import os
import re
import shlex
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DOXYFILE = os.path.join(HERE, "Doxyfile")

# Matches a dependency line entry pointing into the ttcr source directory.
HEADER_RE = re.compile(r"(?:^|/)ttcr/([A-Za-z0-9_]+\.h)$")


# Where a VTK installation's headers usually live. Only consulted when -DVTK was
# requested but the given flags do not make the VTK headers findable.
VTK_GLOBS = [
    "/opt/local/include/vtk-*",      # MacPorts
    "/opt/homebrew/include/vtk-*",   # Homebrew (Apple silicon)
    "/usr/local/include/vtk-*",      # Homebrew (Intel) / manual
    "/usr/include/vtk-*",            # Debian/Ubuntu
    "/usr/include/vtk",
]


def find_vtk_include():
    """Best-effort search for a VTK include directory. Newest version wins."""
    import glob
    hits = []
    for pat in VTK_GLOBS:
        hits.extend(d for d in glob.glob(pat) if os.path.isdir(d))
    if not hits:
        return None
    return sorted(hits)[-1]


def _preprocess(cxx, path, flags):
    res = subprocess.run([cxx, "-std=c++17", *flags, "-MM", path],
                         capture_output=True, text=True)
    return (res.stdout, None) if res.returncode == 0 else (None, res.stderr)


def closure(cxx, srcdir, programs, incflags, vtkflags):
    """Return the set of ttcr headers the given programs include, transitively.

    When VTK is requested the probe must succeed: VTUReader.h is reachable only
    through the VTK code paths, so falling back to a non-VTK preprocess would
    quietly produce a list that is missing a header. That is reported as an
    error rather than papered over.
    """
    if vtkflags:
        probe = os.path.join(srcdir, programs[0])
        out, _ = _preprocess(cxx, probe, vtkflags + incflags)
        if out is None:
            extra = find_vtk_include()
            if extra:
                vtkflags = vtkflags + ["-I" + extra]
                out, err = _preprocess(cxx, probe, vtkflags + incflags)
                if out is None:
                    sys.exit(f"error: VTK requested but its headers are not usable\n"
                             f"       tried -I{extra}\n"
                             f"       set VTKFLAGS to the right path, or "
                             f"VTKFLAGS= to scope the list to a non-VTK build\n"
                             f"{err.strip()[:400]}")
                print(f"  found VTK headers in {extra}", file=sys.stderr)
            else:
                sys.exit("error: VTK requested (-DVTK) but no VTK include "
                         "directory was found.\n"
                         "       Set VTKFLAGS=\"-DVTK -I/path/to/vtk-X.Y\", or "
                         "VTKFLAGS= to scope\n"
                         "       the list to a non-VTK build (drops VTUReader.h).")

    flags = vtkflags + incflags
    found = set()
    for prog in programs:
        path = os.path.join(srcdir, prog)
        if not os.path.isfile(path):
            sys.exit(f"error: no such program source: {path}")
        out, err = _preprocess(cxx, path, flags)
        if out is None:
            sys.exit(f"error: {cxx} could not preprocess {prog}; "
                     f"check CXX and INCFLAGS\n{err.strip()[:400]}")
        for tok in out.replace("\\", " ").split():
            m = HEADER_RE.search(tok)
            if m:
                found.add(m.group(1))
    return found


def render(headers):
    lead, cont = "INPUT                  = ", " " * 25
    lines = []
    for i, h in enumerate(sorted(headers, key=str.lower)):
        prefix = lead if i == 0 else cont
        suffix = "" if i == len(headers) - 1 else " \\"
        lines.append(f"{prefix}../../ttcr/{h}{suffix}")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cxx", default=os.environ.get("CXX", "c++"))
    ap.add_argument("--srcdir", default=os.path.join(HERE, "..", "..", "ttcr"))
    ap.add_argument("--programs", default="ttcr3d.cpp ttcr2d.cpp ttcr2ds.cpp")
    ap.add_argument("--incflags", default="")
    ap.add_argument("--vtkflags", default="")
    ap.add_argument("--check", action="store_true",
                    help="report drift without rewriting the Doxyfile")
    args = ap.parse_args()

    srcdir = os.path.normpath(args.srcdir)
    incflags = shlex.split(args.incflags) or ["-I" + srcdir]
    vtkflags = shlex.split(args.vtkflags)
    programs = args.programs.split()

    headers = closure(args.cxx, srcdir, programs, incflags, vtkflags)
    if not headers:
        sys.exit("error: empty include closure; check CXX and INCFLAGS")

    src = open(DOXYFILE).read()
    block = render(headers)
    new, n = re.subn(r"^INPUT                  = .*?(?=\nINPUT_ENCODING)",
                     lambda _: block, src, count=1, flags=re.S | re.M)
    if n != 1:
        sys.exit("error: could not locate the INPUT block in Doxyfile")

    present = set(re.findall(r"\.\./\.\./ttcr/([A-Za-z0-9_]+\.h)", src))
    added, removed = sorted(headers - present), sorted(present - headers)

    if args.check:
        if added or removed:
            for f in added:
                print(f"  + {f} (in the closure, missing from Doxyfile)")
            for f in removed:
                print(f"  - {f} (in Doxyfile, not in the closure)")
            sys.exit(1)
        print(f"Doxyfile INPUT is up to date ({len(headers)} headers)")
        return

    if new == src:
        print(f"Doxyfile INPUT already up to date ({len(headers)} headers)")
        return
    open(DOXYFILE, "w").write(new)
    for f in added:
        print(f"  + {f}")
    for f in removed:
        print(f"  - {f}")
    print(f"Doxyfile INPUT updated: {len(headers)} headers")


if __name__ == "__main__":
    main()
