import pytest
import stk
import numpy as np

from .case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            dejsonizer=stk.ConstructedMoleculeDejsonizer(),
            molecule_json={
                'a': (
                    (35, 0),
                    (6, 2),
                    (6, 2),
                    (6, 2),
                    (6, 2),
                    (35, 0),
                ),
                'b': (
                    (0, 1, 1, (0, 0, 0)),
                    (1, 2, 1, (0, 0, 0)),
                    (3, 4, 1, (0, 0, 0)),
                    (4, 5, 1, (0, 0, 0)),
                    (2, 3, 1, (0, 0, 0)),
                ),
                'InChI': 'InChI=1S/C4Br2/c5-3-1-2-4-6/q+8',
                'InChIKey': 'CXAFVTYJXQJZSL-UHFFFAOYSA-N',
            },
            constructed_molecule_json={
                'InChI': 'InChI=1S/C4Br2/c5-3-1-2-4-6/q+8',
                'InChIKey': 'CXAFVTYJXQJZSL-UHFFFAOYSA-N',
                'BB': (
                ),
                'nBB': (2, ),
                'aI': (
                    (0, 0),
                    (0, 0),
                    (0, 0),
                    (0, 1),
                    (0, 1),
                    (0, 1),
                ),
                'bI': (
                    (0, 0),
                    (0, 0),
                    (0, 1),
                    (0, 1),
                    (None, None),
                ),
                'ConstructedMoleculeKey': str((
                    'CXAFVTYJXQJZSL-UHFFFAOYSA-N',
                    ('UWAHASCVLDBPQQ-UHFFFAOYSA-N', ),
                    (2, ),
                )),
            },
            position_matrix=np.array([
                [1, 0, 0],
                [2, 0, 0],
                [3, 0, 0],
                [4, 0, 0],
                [5, 0, 0],
                [6, 0, 0],
            ], dtype=np.float64),
            building_blocks=(
                stk.BuildingBlock(
                    smiles='Br[C+2][C+2]Br',
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='Br[C+2][C+2]Br',
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    repeating_unit='A',
                    num_repeating_units=2,
                ),
            ).with_position_matrix(np.array([
                [1, 0, 0],
                [2, 0, 0],
                [3, 0, 0],
                [4, 0, 0],
                [5, 0, 0],
                [6, 0, 0],
            ], dtype=np.float64)),
        ),
    ),
)
def case_data(request):
    return request.param
