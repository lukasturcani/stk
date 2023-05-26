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
            "$periodic 3\n$cell angs\n   36.432   36.431  105.166 "
            " 90.00  90.00  60.00\n$coord angs\n 19.4675 11.2392 5"
            "2.5831 C\n 18.2205 10.5214 53.0065 C\n 16.9674 11.238"
            "9 52.5831 C\n 18.2126 9.0717 52.5831 C\n 20.3559 10.6"
            "937 52.9602 H\n 19.4494 12.2444 53.05 H\n 18.2209 10."
            "5188 54.1253 H\n 16.8933 12.1964 53.1311 H\n 16.1039 "
            "10.6163 52.92 H\n 18.198 8.9414 51.4898 H\n 17.295 8."
            "586 52.9733 H\n 35.18 20.3106 52.5831 C\n 36.427 21.0"
            "285 53.0065 C\n 37.6801 20.3109 52.5831 C\n 36.4349 2"
            "2.4781 52.5831 C\n 34.2916 20.8561 52.9602 H\n 35.198"
            "1 19.3054 53.05 H\n 36.4266 21.031 54.1253 H\n 37.754"
            "2 19.3534 53.1311 H\n 38.5436 20.9335 52.92 H\n 36.44"
            "95 22.6084 51.4898 H\n 37.3525 22.9638 52.9733 H\n 26"
            ".6592 15.3912 52.5831 C\n 27.9883 16.1586 52.5831 C\n"
            " 25.8627 15.9486 53.0853 H\n 26.3889 15.0919 51.5516 "
            "H\n 28.3774 16.2543 53.6163 H\n 27.887 17.1318 52.075"
            "5 H\n 18.2158 0.7674 52.5831 C\n 18.2158 -0.7674 52.5"
            "831 C\n 18.8755 1.1784 53.3532 H\n 18.4472 1.1511 51."
            "5703 H\n 17.762 -1.1522 53.518 H\n 19.2283 -1.1662 52"
            ".4069 H\n 9.7725 15.3912 52.5831 C\n 8.4433 16.1586 5"
            "2.5831 C\n 9.6258 14.315 52.7155 H\n 10.357 15.6362 5"
            "1.6749 H\n 7.9667 16.1026 53.582 H\n 7.7801 15.8075 5"
            "1.7756 H\n$end\n"
        ),
        periodic_info=construction_result.get_periodic_info(),
    )


@pytest.fixture(
    params=(
        lambda: CaseData(
            molecule=stk.BuildingBlock("BrCCBr", [stk.BromoFactory()]),
            writer=stk.TurbomoleWriter(),
            string=(
                "$coord angs\n -1.4238 1.5615 0.3223 Br\n -0.7405 -0.2"
                "573 0.128 C\n 0.7148 -0.1157 -0.3383 C\n 1.6267 0.889"
                "6 1.0687 Br\n -1.3518 -0.8075 -0.5939 H\n -0.7769 -0."
                "6964 1.144 H\n 0.7695 0.528 -1.2387 H\n 1.1821 -1.102"
                "2 -0.4922 H\n$end\n"
            ),
            periodic_info=None,
        ),
        _get_cof_case,
    ),
)
def case_data(request) -> CaseData:
    return request.param()
