# -*- coding: utf-8 -*-
"""VTK is optional: ttcrpy imports and raytraces without it.

Deliberately imports nothing from VTK itself.  The check runs in a child
interpreter where 'import vtk' fails, since this process may already hold it.
"""
import subprocess
import sys
import textwrap
import unittest

SCRIPT = textwrap.dedent('''
    import sys
    sys.modules['vtk'] = None           # any "import vtk" now raises
    import numpy as np
    import ttcrpy.rgrid as rg
    import ttcrpy.tmesh as tm           # imports too

    x = np.arange(6) * 0.1
    g = rg.Grid3d(x, x.copy(), x.copy(), 1, cell_slowness=True, method='FSM')
    g.set_slowness(np.full(125, 0.5))
    tt = g.raytrace(np.array([[0.12, 0.13, 0.14]]),
                    np.array([[0.41, 0.38, 0.33]]))
    assert tt.shape == (1,) and tt[0] > 0, tt

    try:
        g.to_vtk({'s': np.full(125, 0.5)}, 'nowhere')
    except ImportError as e:
        assert 'ttcrpy[vtk]' in str(e), str(e)
    else:
        raise SystemExit('to_vtk ran without VTK')
    print('OK')
''')


class TestVtkIsOptional(unittest.TestCase):

    def test_import_and_raytrace_without_vtk(self):
        r = subprocess.run([sys.executable, '-c', SCRIPT],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn('OK', r.stdout)


if __name__ == '__main__':
    unittest.main()
