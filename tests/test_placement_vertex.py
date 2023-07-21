import numpy as np
import stk


def test_place() -> None:
    bb = stk.BuildingBlock.from_smiles("CC")
    vertex = stk.PlacementVertex(np.array([10, 20, 30]))
    matrix = vertex.place(bb.position_matrix, bb.position_anchor)
    assert np.allclose(stk.get_centroid(matrix), [10, 20, 30])
