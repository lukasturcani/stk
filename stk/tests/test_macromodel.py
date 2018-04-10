"""
Tests functions which use MacroModel.

These tests are only run when the --macromodel py.test option is used.
MacroModel is 3rd party software, it does not come with MMEA.

"""

import pytest
import sys
import os
from os.path import join
import numpy as np
from tempfile import TemporaryDirectory
from .. import macromodel_opt, macromodel_cage_opt, Molecule

macromodel = pytest.mark.skipif(
    all('macromodel' not in x for x in sys.argv),
    reason="only run when explicitly asked")

# Possible installation directories of MacroModel. Your computer's
# must be present in order for this test to run successfully.
dirs = [r'C:\Program Files\Schrodinger2016-3',
        '/home/lukas/program_files/schrodinger2016-4',
        '/opt/schrodinger2016-4',
        '/home/lukas/opt/schrodinger2016-4',
        '/opt/schrodinger2017-2']
mm_path = next((x for x in dirs if os.path.exists(x)), None)

c1 = Molecule.load(join('data', 'macromodel', 'cage.json'))
c2 = Molecule.load(join('data', 'macromodel', 'small_mol.json'))
outdir = 'macromodel_tests_output'
try:
    os.mkdir(outdir)
except:
    ...


@macromodel
def test_macromodel_opt():
    if outdir not in os.getcwd():
        os.chdir(outdir)

    macromodel_opt(c1, mm_path,
    {'md' : True, 'gradient' : 1, 'restricted' : False},
    {'gradient' : 1, 'sim_time' : 20, 'eq_time' : 2, 'confs' : 2})


@macromodel
def test_macromodel_cage_opt():
    if outdir not in os.getcwd():
        os.chdir(outdir)

    macromodel_cage_opt(c1, mm_path,
    {'md' : True, 'gradient' : 1, 'restricted' : False},
    {'gradient' : 1, 'sim_time' : 20, 'eq_time' : 2, 'confs' : 2})


@macromodel
def test_macromodel_eng():
    if outdir not in os.getcwd():
        os.chdir(outdir)
    assert np.allclose(
        c2.energy.macromodel(16, mm_path), 23.48, atol=1e-2)
