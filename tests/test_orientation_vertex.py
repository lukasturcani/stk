from dataclasses import replace

import numpy as np
import rdkit.Chem.AllChem as rdkit
import stk


def test_place() -> None:
    bb = stk.BuildingBlock.from_smiles("BrC(CCC)CCCBr", stk.bromo())
    br1, br2 = bb.position_matrix[[0, 8]]
    bb = replace(bb, rotation_anchor=stk.RotationAnchor(br2 - br1, None))
    assert bb.rotation_anchor is not None
    vertex = stk.OrientationVertex(np.array([10, 20, 30]), np.array([0, 1, 0]))
    matrix = vertex.place(
        matrix=bb.position_matrix,
        position_anchor=bb.position_anchor,
        rotation_anchor_axis=bb.rotation_anchor.axis,
    )
    br1, br2 = matrix[[0, 8]]
    assert np.allclose(stk.normalize_vector(br2 - br1), [0, 1, 0], atol=1e-6)
