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
            "$periodic 3\n$cell angs\n   36.348   36.347  104.924  90.00  "
            "90.00  60.00\n$coord angs\n 19.4236 11.2139 52.4621 C\n 18.1939 "
            "10.4901 52.8934 C\n 16.929 11.2242 52.4621 C\n 18.1692 9.0391 52"
            ".4621 C\n 20.3081 10.6461 52.8138 H\n 19.4563 12.2214 52.9174 H"
            "\n 18.1856 10.4831 54.0041 H\n 16.871 12.2049 52.9835 H\n 16.0593"
            " 10.6426 52.8347 H\n 18.0564 8.8975 51.3812 H\n 17.2567 8.589 52."
            "9526 H\n 35.0982 20.2633 52.4621 C\n 36.3279 20.9871 52.8934 C\n "
            "37.5928 20.253 52.4621 C\n 36.3526 22.4381 52.4621 C\n 34.2136 20"
            ".8312 52.8138 H\n 35.0654 19.2558 52.9174 H\n 36.3362 20.9941 54."
            "0041 H\n 37.6508 19.2723 52.9835 H\n 38.4624 20.8346 52.8347 H\n "
            "36.4653 22.5797 51.3812 H\n 37.265 22.8883 52.9526 H\n 26.6019 15"
            ".3582 52.4621 C\n 27.9198 16.1191 52.4621 C\n 26.6425 14.396 51.9"
            "319 H\n 26.2373 15.2066 53.5136 H\n 28.1283 16.4845 51.4331 H\n 2"
            "8.7543 15.5207 52.8863 H\n 18.1739 0.7609 52.4621 C\n 18.1739 -0."
            "7609 52.4621 C\n 17.2162 1.2069 52.158 H\n 18.4819 1.1525 53.4688"
            " H\n 18.1266 -1.1242 51.4125 H\n 17.3715 -1.1844 53.1033 H\n 9.74"
            "59 15.3582 52.4621 C\n 8.428 16.1191 52.4621 C\n 10.2575 15.3524 "
            "51.4891 H\n 10.4298 15.7596 53.2576 H\n 7.6812 15.5521 51.8647 H"
            "\n 8.5459 17.1704 52.1226 H\n$end\n"
        ),
        periodic_info=construction_result.get_periodic_info(),  # type: ignore[attr-defined]
    )


@pytest.fixture(
    params=(
        lambda: CaseData(
            molecule=stk.BuildingBlock("BrCCBr", [stk.BromoFactory()]),
            writer=stk.TurbomoleWriter(),
            string=(
                "$coord angs\n -1.3795 1.5853 0.2769 Br\n -0.7294 -0.2312 "
                "0.1386 C\n 0.71 -0.1421 -0.3475 C\n 1.6335 0.8945 1.0205 "
                "Br\n -1.3861 -0.8197 -0.5178 H\n -0.7646 -0.6353 1.1859 "
                "H\n 0.7466 0.4912 -1.2604 H\n 1.1695 -1.1427 -0.4961 H\n"
                "$end\n"
            ),
            periodic_info=None,
        ),
        _get_cof_case,
    ),
)
def case_data(request) -> CaseData:
    return request.param()
