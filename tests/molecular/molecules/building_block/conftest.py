import stk
import itertools as it
import pytest


@pytest.fixture(
    params=[
        stk.BuildingBlock('NCCN'),
        stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()]),
        stk.BuildingBlock(
            smiles='BrCC(Br)C(Br)C(Br)C(Br)C(Br)C(Br)CBr',
            functional_groups=[stk.BromoFactory()],
        ),
        stk.BuildingBlock('N[C+][C+2]N'),
    ],
)
def building_block(request):
    """
    A :class:`.BuildingBlock` instance.

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
