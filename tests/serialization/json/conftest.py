import stk
import pytest

from .case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            jsonizer=stk.MoleculeJsonizer(
                key_makers=(
                    stk.Inchi(),
                    stk.InchiKey(),
                )
            ),
            molecule=stk.BuildingBlock('Br[C+2][C+2]Br'),
            json={
                'a': (
                    (35, 0),
                    (6, 2),
                    (6, 2),
                    (35, 0),
                ),
                'b': (
                    (0, 1, 1, (0, 0, 0)),
                    (1, 2, 1, (0, 0, 0)),
                    (2, 3, 1, (0, 0, 0)),
                ),
                'InChI': 'InChI=1S/C2Br2/c3-1-2-4/q+4',
                'InChIKey': 'UWAHASCVLDBPQQ-UHFFFAOYSA-N',
            },
        ),
        CaseData(
            jsonizer=stk.ConstructedMoleculeJsonizer(
                key_makers=(
                    stk.Inchi(),
                    stk.InchiKey(),
                    stk.ConstructedMoleculeKeyMaker(),
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
            ),
            json={
                'InChI': 'InChI=1S/C4Br2/c5-3-1-2-4-6/q+8',
                'InChIKey': 'CXAFVTYJXQJZSL-UHFFFAOYSA-N',
                'BB': (
                    {
                        'InChI': 'InChI=1S/C2Br2/c3-1-2-4/q+4',
                        'InChIKey': 'UWAHASCVLDBPQQ-UHFFFAOYSA-N',
                    },
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
        ),
    ),
)
def case_data(request):
    return request.param
