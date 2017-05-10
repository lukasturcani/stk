"""
Tests functions which use MOPAC.

These tests are only run when the --macromodel py.test option is used.
MacroModel is 3rd party software, it does not come with MMEA.

"""

import pytest
import sys
import os
from os.path import join
import numpy as np
from .. import mopac_opt, Molecule

mopac = pytest.mark.skipif(
    all('mopac' not in x for x in sys.argv),
    reason="only run when explicitly asked")

# Possible installation directories of MacroModel. Your computer's
# must be present in order for this test to run successfully.
dirs = [r'C:\Program Files\mopac\MOPAC2016.exe',
        '/opt/mopac/MOPAC2016.exe']
mopac_path = next((x for x in dirs if os.path.exists(x)), None)

c1 = Molecule.load(join('data', 'mopac', 'small_mol.json'))
c2 = Molecule.load(join('data', 'mopac', 'small_mol2.json'))
outdir = 'mopac_tests_output'
try:
    os.mkdir(outdir)
except:
    ...


@mopac
def test_mopac_opt():
    if outdir not in os.getcwd():
        os.chdir(outdir)

    mopac_opt(c1, mopac_path)


@mopac
def test_mopac_ip():
    if outdir not in os.getcwd():
        os.chdir(outdir)
    assert np.allclose(
        c2.energy.mopac_ip(mopac_path), 5.3975, atol=1e-2)
