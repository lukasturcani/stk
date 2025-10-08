import pytest

import stk

from .case_data import CaseData


def _get_cof_case() -> CaseData:
    bb1 = stk.BuildingBlock("BrCCBr", [stk.BromoFactory()])
    bb2 = stk.BuildingBlock("BrCC(CBr)CBr", [stk.BromoFactory()])
    topology_graph = stk.cof.PeriodicHoneycomb(
        building_blocks=(bb1, bb2),
        lattice_size=(1, 1, 1),
    )
    construction_result = topology_graph.construct()
    cof = stk.ConstructedMolecule.init_from_construction_result(
        construction_result=construction_result,
    )
    return CaseData(
        molecule=cof,
        writer=stk.TurbomoleWriter(),
        string=(
            "$periodic 3\n$cell angs\n   36.393   36.392  105.054  90.00  90.0"
            "0  60.00\n$coord angs\n 19.448 11.228 52.5268 C\n 18.2028 10.5111"
            " 52.9458 C\n 16.948 11.2277 52.5268 C\n 18.193 9.0604 52.5268 C\n"
            " 20.3376 10.6938 52.9171 H\n 19.4166 12.2371 52.9845 H\n 18.1974 "
            "10.5145 54.0612 H\n 16.857 12.1765 53.0679 H\n 16.0939 10.5852 52"
            ".8465 H\n 18.1737 8.915 51.4353 H\n 17.2817 8.5734 52.9302 H\n 35"
            ".141 20.2881 52.5268 C\n 36.3862 21.0049 52.9458 C\n 37.641 20.28"
            "83 52.5268 C\n 36.396 22.4557 52.5268 C\n 34.2514 20.8222 52.9171"
            " H\n 35.1724 19.2789 52.9845 H\n 36.3916 21.0015 54.0612 H\n 37.7"
            "32 19.3396 53.0679 H\n 38.4951 20.9308 52.8465 H\n 36.4153 22.601"
            " 51.4353 H\n 37.3073 22.9426 52.9302 H\n 26.63 15.3744 52.5268 C"
            "\n 27.959 16.1417 52.5268 C\n 26.6792 14.4269 51.9881 H\n 26.2573"
            " 15.2518 53.5624 H\n 28.1979 16.498 51.5049 H\n 28.7707 15.5371 "
            "52.9674 H\n 18.1963 0.7673 52.5268 C\n 18.1963 -0.7673 52.5268 C"
            "\n 17.2418 1.1984 52.2215 H\n 18.538 1.1514 53.5077 H\n 18.1186 -"
            "1.1524 51.4904 H\n 17.41 -1.1679 53.1898 H\n 9.7627 15.3744 52.52"
            "68 C\n 8.4337 16.1417 52.5268 C\n 10.2629 15.3786 51.5572 H\n 10."
            "4179 15.741 53.3409 H\n 7.6835 15.6125 51.9058 H\n 8.5773 17.1917"
            " 52.2182 H\n$end\n"
        ),
        periodic_info=construction_result.get_periodic_info(),  # type: ignore[attr-defined]
    )


@pytest.fixture(
    params=(
        lambda: CaseData(
            molecule=stk.BuildingBlock("BrCCBr", [stk.BromoFactory()]),
            writer=stk.TurbomoleWriter(),
            string=(
                "$coord angs\n -1.4302 1.5499 0.3629 Br\n -0.7383 -0.2592 0.11"
                "57 C\n 0.7182 -0.1112 -0.3446 C\n 1.6324 0.852 1.0901 Br\n -1"
                ".3636 -0.7969 -0.5986 H\n -0.7775 -0.7041 1.1291 H\n 0.7784 0"
                ".5649 -1.2208 H\n 1.1805 -1.0955 -0.5339 H\n$end\n"
            ),
            periodic_info=None,
        ),
        _get_cof_case,
    ),
)
def case_data(request) -> CaseData:
    return request.param()
