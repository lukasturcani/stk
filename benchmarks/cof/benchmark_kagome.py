from __future__ import annotations

import pytest
import stk


def build_kagome(
    lattice_size: tuple[int, int, int],
) -> stk.ConstructedMolecule:
    bb1 = stk.BuildingBlock("BrCCBr", [stk.BromoFactory()])
    bb2 = stk.BuildingBlock(
        smiles="BrC1C(Br)CC(Br)C(Br)C1",
        functional_groups=[stk.BromoFactory()],
    )
    cof = stk.ConstructedMolecule(
        topology_graph=stk.cof.Kagome(
            building_blocks=(bb1, bb2),
            lattice_size=lattice_size,
        ),
    )
    return cof


@pytest.fixture(
    params=(
        (1, 1, 1),
        (2, 2, 2),
        (4, 4, 4),
    ),
)
def lattice_size(request) -> tuple[int, int, int]:
    return request.param


def benchmark_kagome(
    benchmark,
    lattice_size: tuple[int, int, int],
) -> None:
    benchmark(build_kagome, lattice_size)
