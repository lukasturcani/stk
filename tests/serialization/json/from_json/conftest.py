import numpy as np
import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda: CaseData(
            dejsonizer=stk.MoleculeDejsonizer(),
            json={
                "molecule": {
                    "a": (
                        (35, 0),
                        (6, 2),
                        (6, 2),
                        (35, 0),
                    ),
                    "b": (
                        (0, 1, 1, (0, 0, 0)),
                        (1, 2, 1, (0, 0, 0)),
                        (2, 3, 1, (0, 0, 0)),
                    ),
                    "InChI": "InChI=1S/C2Br2/c3-1-2-4/q+4",
                    "InChIKey": "UWAHASCVLDBPQQ-UHFFFAOYSA-N",
                },
                "matrix": {
                    "InChI": "InChI=1S/C2Br2/c3-1-2-4/q+4",
                    "InChIKey": "UWAHASCVLDBPQQ-UHFFFAOYSA-N",
                    "m": [
                        [0.0, 0.0, 0.0],
                        [1.0, 0.0, 0.0],
                        [2.0, 0.0, 0.0],
                        [3.0, 0.0, 0.0],
                    ],
                },
            },
            molecule=stk.BuildingBlock(
                smiles="Br[C+2][C+2]Br",
            ).with_position_matrix(
                np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [1.0, 0.0, 0.0],
                        [2.0, 0.0, 0.0],
                        [3.0, 0.0, 0.0],
                    ]
                )
            ),
        ),
        lambda: CaseData(
            dejsonizer=stk.ConstructedMoleculeDejsonizer(),
            json={
                "molecule": {
                    "a": (
                        (35, 0),
                        (6, 2),
                        (6, 2),
                        (6, 2),
                        (6, 2),
                        (35, 0),
                    ),
                    "b": (
                        (0, 1, 1, (0, 0, 0)),
                        (1, 2, 1, (0, 0, 0)),
                        (3, 4, 1, (0, 0, 0)),
                        (4, 5, 1, (0, 0, 0)),
                        (2, 3, 1, (0, 0, 0)),
                    ),
                    "InChI": "InChI=1S/C4Br2/c5-3-1-2-4-6/q+8",
                    "InChIKey": "CXAFVTYJXQJZSL-UHFFFAOYSA-N",
                },
                "constructedMolecule": {
                    "InChI": "InChI=1S/C4Br2/c5-3-1-2-4-6/q+8",
                    "InChIKey": "CXAFVTYJXQJZSL-UHFFFAOYSA-N",
                    "BB": (
                        {
                            "InChI": "InChI=1S/C2Br2/c3-1-2-4/q+4",
                            "InChIKey": "UWAHASCVLDBPQQ-UHFFFAOYSA-N",
                        },
                    ),
                    "nBB": (2,),
                    "aI": (
                        (0, 0, 0),
                        (0, 0, 1),
                        (0, 0, 2),
                        (0, 1, 1),
                        (0, 1, 2),
                        (0, 1, 3),
                    ),
                    "bI": (
                        (0, 0),
                        (0, 0),
                        (0, 1),
                        (0, 1),
                        (None, None),
                    ),
                },
                "matrix": {
                    "m": [
                        [1.0, 0.0, 0.0],
                        [2.0, 0.0, 0.0],
                        [3.0, 0.0, 0.0],
                        [4.0, 0.0, 0.0],
                        [5.0, 0.0, 0.0],
                        [6.0, 0.0, 0.0],
                    ],
                    "InChI": "InChI=1S/C4Br2/c5-3-1-2-4-6/q+8",
                    "InChIKey": "CXAFVTYJXQJZSL-UHFFFAOYSA-N",
                },
                "buildingBlocks": (
                    {
                        "molecule": {
                            "a": (
                                (35, 0),
                                (6, 2),
                                (6, 2),
                                (35, 0),
                            ),
                            "b": (
                                (0, 1, 1, (0, 0, 0)),
                                (1, 2, 1, (0, 0, 0)),
                                (2, 3, 1, (0, 0, 0)),
                            ),
                            "InChI": "InChI=1S/C2Br2/c3-1-2-4/q+4",
                            "InChIKey": "UWAHASCVLDBPQQ-UHFFFAOYSA-N",
                        },
                        "matrix": {
                            "InChI": "InChI=1S/C2Br2/c3-1-2-4/q+4",
                            "InChIKey": "UWAHASCVLDBPQQ-UHFFFAOYSA-N",
                            "m": [
                                [0.0, 0.0, 0.0],
                                [1.0, 0.0, 0.0],
                                [2.0, 0.0, 0.0],
                                [3.0, 0.0, 0.0],
                            ],
                        },
                    },
                ),
            },
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles="Br[C+2][C+2]Br",
                            functional_groups=[stk.BromoFactory()],
                        ).with_position_matrix(
                            np.array(
                                [
                                    [0.0, 0.0, 0.0],
                                    [1.0, 0.0, 0.0],
                                    [2.0, 0.0, 0.0],
                                    [3.0, 0.0, 0.0],
                                ],
                                dtype=np.float64,
                            )
                        ),
                    ),
                    repeating_unit="A",
                    num_repeating_units=2,
                ),
            ).with_position_matrix(
                np.array(
                    [
                        [1, 0, 0],
                        [2, 0, 0],
                        [3, 0, 0],
                        [4, 0, 0],
                        [5, 0, 0],
                        [6, 0, 0],
                    ],
                    dtype=np.float64,
                )
            ),
        ),
    ),
)
def case_data(request) -> CaseData:
    """
    A :class:`.CaseData` instance.

    """

    return request.param()
