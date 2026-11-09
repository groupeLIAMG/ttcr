"""Build a copy of ttcrpy with ThreadSanitizer, for the interpreter running this.

Usage: python tests/tsan/build.py BUILD_DIR

The checkout is left alone: the sources are copied to BUILD_DIR and setup.py is
patched there.  Only the extensions are instrumented, and for arm64 only.
Builds for several interpreters can share a BUILD_DIR.
"""
import os
import pathlib
import platform
import shutil
import subprocess
import sys

REPO = pathlib.Path(__file__).resolve().parents[2]
IGNORED = shutil.ignore_patterns('.git', 'build', 'dist', 'tests', '__pycache__',
                                 '*.so', 'rgrid.cpp', 'tmesh.cpp', '_version.py')
PATCHES = (
    ("'-O3', '-arch', 'x86_64',",
     "'-O1', '-g', '-fsanitize=thread', '-fno-omit-frame-pointer',"),
    ("'-undefined', 'dynamic_lookup', '-arch', 'x86_64',",
     "'-undefined', 'dynamic_lookup', '-fsanitize=thread',"),
)


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    if platform.system() != 'Darwin' or platform.machine() != 'arm64':
        sys.exit('only written for macOS on arm64')

    build_dir = pathlib.Path(sys.argv[1]).resolve()
    if build_dir.is_relative_to(REPO):
        sys.exit('BUILD_DIR must be outside the checkout')
    shutil.copytree(REPO, build_dir, ignore=IGNORED, dirs_exist_ok=True)

    setup = (REPO / 'setup.py').read_text()
    for old, new in PATCHES:
        if setup.count(old) != 1:
            sys.exit('setup.py has changed, cannot patch: %s' % old)
        setup = setup.replace(old, new)
    (build_dir / 'setup.py').write_text(setup)

    env = dict(os.environ, ARCHFLAGS='-arch arm64',
               SETUPTOOLS_SCM_PRETEND_VERSION='0.0.0')
    env.pop('PYTHONPATH', None)
    subprocess.run([sys.executable, 'setup.py', 'build_ext', '--inplace', '--force',
                    '--parallel', '4'], cwd=build_dir, env=env, check=True)


if __name__ == '__main__':
    main()
