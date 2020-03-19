import stk
import itertools as it
import pytest


from .case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            building_block=stk.BuildingBlock('Br[C+2][C+2]Br'),
            functional_groups=(),
            core_atoms=(stk.Br(0), stk.C(1), stk.C(2), stk.Br(3)),
            placers=(stk.Br(0), stk.C(1), stk.C(2), stk.Br(3)),
        ),
        CaseData(
            building_block=stk.BuildingBlock(
                smiles='Br[C+2][C+2]Br',
                functional_groups=[stk.BromoFactory()],
            ),
            functional_groups=(
                stk.Bromo(
                    bromine=stk.Br(0),
                    atom=stk.C(1),
                    bonders=(stk.C(1), ),
                    deleters=(stk.Br(0), ),
                ),
                stk.Bromo(
                    bromine=stk.Br(3),
                    atom=stk.C(2),
                    bonders=(stk.C(2), ),
                    deleters=(stk.Br(3), ),
                ),
            ),
            core_atoms=(stk.C(1), stk.C(2)),
            placers=(stk.C(1), stk.C(2)),
        ),
    ),
)
def case_data(request):
    """
    A :class:`.CaseData` instance.

    """

    return request.param.clone()


@pytest.fixture(
    params=(
        lambda molecule:
            stk.BromoFactory().get_functional_groups(molecule),
        lambda molecule:
            stk.PrimaryAminoFactory().get_functional_groups(molecule),
        lambda molecule: it.chain(
            stk.PrimaryAminoFactory().get_functional_groups(molecule),
            stk.BromoFactory().get_functional_groups(molecule)),
    )
)
def get_functional_groups(request):
    """
    Yield the functional groups of a `molecule`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule whose functional groups should be gotten.

    Yields
    ------
    :class:`.FunctionalGroup`
        A functional group of `molecule`.

    """

    return request.param
