"""
Tests functions which use MOPAC.

These tests are only run when the --mopac_path pytest option is used.
MOPAC is 3rd party software, it does not come with ``stk``.

"""

import pytest
import sys
import os
import stk

mopac = pytest.mark.skipif(
    all('mopac' not in x for x in sys.argv),
    reason="Only run when explicitly asked.")

outdir = 'mopac_tests_output'
if not os.path.exists(outdir):
    os.mkdir(outdir)


@mopac
def test_mopac_opt(tmp_amine2, mopac_path):
    # Give conformer a distinct geometry.
    tmp_amine2.set_position_from_matrix(
        pos_mat=tmp_amine2.mol.GetConformer().GetPositions().T*4,
        conformer=1)

    tmp_amine2.write(os.path.join(outdir, 'before_opt.mol'))
    mopac = stk.MOPAC(mopac_path)
    mopac.optimize(tmp_amine2)
    tmp_amine2.write(os.path.join(outdir, 'after_opt.mol'))


@mopac
def test_mopac_ip(amine2, mopac_path):
    if outdir not in os.getcwd():
        os.chdir(outdir)
    assert abs(amine2.energy.mopac_ip(mopac_path)-5.3975) < 1e02
