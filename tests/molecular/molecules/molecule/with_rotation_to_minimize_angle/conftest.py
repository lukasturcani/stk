import pytest
import stk


@pytest.fixture(
    scope="session",
    params=(
        lambda: stk.BuildingBlock("NCCN"),
        lambda: stk.BuildingBlock(
            smiles="Brc1ccc(Br)cc1Br",
            functional_groups=[stk.BromoFactory()],
        ),
        lambda: stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    stk.BuildingBlock("BrCCBr", [stk.BromoFactory()]),
                    stk.BuildingBlock(
                        smiles="BrCNCCBr",
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                repeating_unit="AB",
                num_repeating_units=2,
            ),
        ),
    ),
)
def molecule(request) -> stk.Molecule:
    """
    A :class:`.Molecule` instance with at least 3 atoms.

    """

    return request.param()
