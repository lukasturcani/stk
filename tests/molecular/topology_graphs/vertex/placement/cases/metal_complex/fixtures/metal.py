import pytest
import numpy as np
import stk
from rdkit.Chem import AllChem as rdkit

from ....case_data import CaseData

vertices = stk.metal_complex.vertices

single_atom = rdkit.MolFromSmiles('[Fe]')
single_atom.AddConformer(rdkit.Conformer(single_atom.GetNumAtoms()))


@pytest.fixture(
    params=(
        CaseData(
            vertex=vertices._MetalVertex(
                id=0,
                position=(1, 2, 3),
            ),
            edges=(),
            building_block=stk.BuildingBlock.init_from_rdkit_mol(
                single_atom
            ),
            position=np.array([1, 2, 3], dtype=np.float64),
            alignment_tests={},
            functional_group_edges={},
        ),
    ),
)
def metal(request):
    return request.param
