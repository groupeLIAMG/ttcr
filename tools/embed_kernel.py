#!/usr/bin/env python
"""Embed an OpenCL .cl file into a C++ header as a raw string literal.

The Unix makefiles do this with printf and cat.  Makefile.win cannot: in a
Visual Studio Developer Command Prompt the shell is cmd.exe, which has neither,
and the header text contains quotes and '#' that cmd's echo will not pass
through intact.  This does the same job portably.

The output matches setup.py's _embed_kernel and the makefile recipes byte for
byte, apart from the provenance comment on the first line.

    python tools/embed_kernel.py <in.cl> <out.h> <symbol> <delimiter>
"""

import io
import os
import sys


def embed(cl_path, hdr_path, sym_name, delim):
    with io.open(cl_path, 'r', encoding='utf-8') as f:
        src = f.read()

    # A raw string ends at the first )delim" -- if the kernel contains that
    # sequence the generated header would not compile, so refuse rather than
    # emit something broken.
    if (')' + delim + '"') in src:
        sys.exit('%s: kernel source contains the raw-string delimiter %s; '
                 'choose another' % (cl_path, delim))

    guard = 'ttcr_' + os.path.basename(hdr_path).replace('.', '_')
    content = (
        '// Auto-generated from %s. Do not edit.\n' % os.path.basename(cl_path) +
        '#ifndef %s\n' % guard +
        '#define %s\n' % guard +
        'namespace ttcr {\n' +
        'static const char* const %s = R"%s(\n%s)%s";\n' % (sym_name, delim, src, delim) +
        '}\n'
        '#endif\n'
    )
    with io.open(hdr_path, 'w', encoding='utf-8', newline='\n') as f:
        f.write(content)


def main(argv):
    if len(argv) != 5:
        sys.exit('usage: embed_kernel.py <in.cl> <out.h> <symbol> <delimiter>')
    embed(argv[1], argv[2], argv[3], argv[4])


if __name__ == '__main__':
    main(sys.argv)
