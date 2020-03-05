import pytest
from pytest_lazyfixture import lazy_fixture
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            building_blocks=(
                stk.BuildingBlock('BrC#CBr', [stk.BromoFactory()]),
                stk.BuildingBlock(
                    smiles='Br[C+]=NC#CBr',
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            topology_graph=stk.polymer.Linear('AB', 2),
            num_new_atoms=0,
            num_new_bonds=3,
        ),
    ),
)
def linear(request):
    return request.param


@pytest.fixture(
    params=(
        CaseData(
            building_blocks=(
                stk.BuildingBlock(
                    smiles='Br[C+]=NC#CBr',
                    functional_groups=[stk.BromoFactory()],
                ),
                stk.BuildingBlock(
                    smiles=(
                        'Br[C+]1[C+2][C+](Br)[C+](F)[C+](Br)[C+2]1'
                    ),
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            topology_graph=stk.cof.Honeycomb((2, 2, 1)),
            num_new_atoms=...,
            num_new_bonds=...
        ),
        CaseData(
            building_blocks=(
                stk.BuildingBlock(
                    smiles='Br[C+]=NC#CBr',
                    functional_groups=[stk.BromoFactory()],
                ),
                stk.BuildingBlock(
                    smiles=(
                        'Br[C+]1[C+2][C+](Br)[C+](F)[C+](Br)[C+2]1'
                    ),
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            topology_graph=stk.cof.Honeycomb(
                lattice_size=(2, 2, 1),
                vertex_alignments={0: 1, 1: 1, 2: 1, 3: 1, 4: 1},
            ),
            num_new_atoms=...,
            num_new_bonds=...,
        ),
        CaseData(
            building_blocks=(
                stk.BuildingBlock(
                    smiles='Br[C+]=NC#CBr',
                    functional_groups=[stk.BromoFactory()],
                ),
                stk.BuildingBlock(
                    smiles=(
                        'Br[C+]1[C+2][C+](Br)[C+](F)[C+](Br)[C+2]1'
                    ),
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            topology_graph=stk.cof.Honeycomb(
                lattice_size=(2, 2, 1),
                vertex_alignments={0: 2, 1: 2},
            ),
            num_new_atoms=...,
            num_new_bonds=...,
        ),
        CaseData(
            building_blocks=(
                stk.BuildingBlock(
                    smiles=(
                        'Br[C+]1[C+2][C+](Br)[C+](F)[C+](Br)[C+2]1'
                    ),
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            topology_graph=stk.cof.LinkerlessHoneycomb((2, 2, 1)),
            num_new_atoms=...,
            num_new_bonds=...,
        ),
        CaseData(
            building_blocks=(
                stk.BuildingBlock(
                    smiles='Br[C+]=NC#CBr',
                    functional_groups=[stk.BromoFactory()],
                ),
                stk.BuildingBlock(
                    smiles=(
                        'Br[C+]1[C+](F)[C+](I)[C+](I)[C+](Br)'
                        'C1Br'
                    ),
                    functional_groups=[
                        stk.BromoFactory(),
                        stk.IodoFactory(),
                        stk.FluoroFactory(),
                    ],
                ),
            ),
            topology_graph=stk.cof.Hexagonal(
                lattice_size=(2, 2, 1),
                vertex_alignments={0: 5},
            ),
            num_new_atoms=...,
            num_new_bonds=...,
        ),
        CaseData(
            building_blocks=(
                stk.BuildingBlock(
                    smiles='Br[C+]=NC#CBr',
                    functional_groups=[stk.BromoFactory()],
                ),
                stk.BuildingBlock(
                    smiles=(
                        'Br[C+]1[C+](Br)[C+](F)[C+](Br)[C+](Br)'
                        '[C+2]1'
                    ),
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            topology_graph=stk.cof.Kagome((2, 2, 1)),
            num_new_atoms=...,
            num_new_bonds=...,
        ),
        CaseData(
            building_blocks=(
                stk.BuildingBlock(
                    smiles='Br[C+]=NC#CBr',
                    functional_groups=[stk.BromoFactory()],
                ),
                stk.BuildingBlock(
                    smiles='BrC1=C(Br)C(F)(Br)[C+]1Br',
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            topology_graph=stk.cof.Square((2, 2, 1)),
            num_new_atoms=...,
            num_new_bonds=...,
        ),
    ),
)
def cof(request):
    return request.param


@pytest.fixture(
    params=(
        lazy_fixture('linear'),
        lazy_fixture('cof'),
    ),
)
def case_data(request):
    return request.param
