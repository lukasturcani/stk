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
                'atoms': (
                    {'atomic_number': 35, 'charge': 0},
                    {'atomic_number': 6, 'charge': 2},
                    {'atomic_number': 6, 'charge': 2},
                    {'atomic_number': 35, 'charge': 0},
                ),
                'bonds': (
                    {
                        'atom1': 0,
                        'atom2': 1,
                        'order': 1,
                        'periodicity': (0, 0, 0),
                    },
                    {
                        'atom1': 1,
                        'atom2': 2,
                        'order': 1,
                        'periodicity': (0, 0, 0),
                    },
                    {
                        'atom1': 2,
                        'atom2': 3,
                        'order': 1,
                        'periodicity': (0, 0, 0),
                    },
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
                'building_blocks': (
                    {
                        'InChI': 'InChI=1S/C2Br2/c3-1-2-4/q+4',
                        'InChIKey': 'UWAHASCVLDBPQQ-UHFFFAOYSA-N',
                    },
                ),
                'num_building_blocks': {0: 2},
                'atom_infos': (
                    {'building_block': 0, 'building_block_id': 0},
                    {'building_block': 0, 'building_block_id': 0},
                    {'building_block': 0, 'building_block_id': 0},
                    {'building_block': 0, 'building_block_id': 1},
                    {'building_block': 0, 'building_block_id': 1},
                    {'building_block': 0, 'building_block_id': 1},
                ),
                'bond_infos': (
                    {'building_block': 0, 'building_block_id': 0},
                    {'building_block': 0, 'building_block_id': 0},
                    {'building_block': 0, 'building_block_id': 1},
                    {'building_block': 0, 'building_block_id': 1},
                    {
                        'building_block': None,
                        'building_block_id': None,
                    },
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
