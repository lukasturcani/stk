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
        writer=stk.PdbWriter(),
        string=(
            "CRYST1   36.348   36.347  104.924  90.00  90.00  60.00 P 1\nHETAT"
            "M    1 C1   UNL     1      19.424  11.214  52.462  1.00  0.00    "
            "       C 0\nHETATM    2 C2   UNL     1      18.194  10.490  52.89"
            "3  1.00  0.00           C 0\nHETATM    3 C3   UNL     1      16.9"
            "29  11.224  52.462  1.00  0.00           C 0\nHETATM    4 C4   UN"
            "L     1      18.169   9.039  52.462  1.00  0.00           C 0\nHE"
            "TATM    5 H1   UNL     1      20.308  10.646  52.814  1.00  0.00 "
            "          H 0\nHETATM    6 H2   UNL     1      19.456  12.221  52"
            ".917  1.00  0.00           H 0\nHETATM    7 H3   UNL     1      1"
            "8.186  10.483  54.004  1.00  0.00           H 0\nHETATM    8 H4  "
            " UNL     1      16.871  12.205  52.984  1.00  0.00           H 0"
            "\nHETATM    9 H5   UNL     1      16.059  10.643  52.835  1.00  0"
            ".00           H 0\nHETATM   10 H6   UNL     1      18.056   8.898"
            "  51.381  1.00  0.00           H 0\nHETATM   11 H7   UNL     1   "
            "   17.257   8.589  52.953  1.00  0.00           H 0\nHETATM   12 "
            "C5   UNL     1      35.098  20.263  52.462  1.00  0.00           "
            "C 0\nHETATM   13 C6   UNL     1      36.328  20.987  52.893  1.00"
            "  0.00           C 0\nHETATM   14 C7   UNL     1      37.593  20."
            "253  52.462  1.00  0.00           C 0\nHETATM   15 C8   UNL     1"
            "      36.353  22.438  52.462  1.00  0.00           C 0\nHETATM   "
            "16 H8   UNL     1      34.214  20.831  52.814  1.00  0.00        "
            "   H 0\nHETATM   17 H9   UNL     1      35.065  19.256  52.917  1"
            ".00  0.00           H 0\nHETATM   18 H10  UNL     1      36.336  "
            "20.994  54.004  1.00  0.00           H 0\nHETATM   19 H11  UNL   "
            "  1      37.651  19.272  52.984  1.00  0.00           H 0\nHETATM"
            "   20 H12  UNL     1      38.462  20.835  52.835  1.00  0.00     "
            "      H 0\nHETATM   21 H13  UNL     1      36.465  22.580  51.381"
            "  1.00  0.00           H 0\nHETATM   22 H14  UNL     1      37.26"
            "5  22.888  52.953  1.00  0.00           H 0\nHETATM   23 C9   UNL"
            "     1      26.602  15.358  52.462  1.00  0.00           C 0\nHET"
            "ATM   24 C10  UNL     1      27.920  16.119  52.462  1.00  0.00  "
            "         C 0\nHETATM   25 H15  UNL     1      26.642  14.396  51."
            "932  1.00  0.00           H 0\nHETATM   26 H16  UNL     1      26"
            ".237  15.207  53.514  1.00  0.00           H 0\nHETATM   27 H17  "
            "UNL     1      28.128  16.485  51.433  1.00  0.00           H 0\n"
            "HETATM   28 H18  UNL     1      28.754  15.521  52.886  1.00  0.0"
            "0           H 0\nHETATM   29 C11  UNL     1      18.174   0.761  "
            "52.462  1.00  0.00           C 0\nHETATM   30 C12  UNL     1     "
            " 18.174  -0.761  52.462  1.00  0.00           C 0\nHETATM   31 H1"
            "9  UNL     1      17.216   1.207  52.158  1.00  0.00           H "
            "0\nHETATM   32 H20  UNL     1      18.482   1.152  53.469  1.00  "
            "0.00           H 0\nHETATM   33 H21  UNL     1      18.127  -1.12"
            "4  51.412  1.00  0.00           H 0\nHETATM   34 H22  UNL     1  "
            "    17.372  -1.184  53.103  1.00  0.00           H 0\nHETATM   35"
            " C13  UNL     1       9.746  15.358  52.462  1.00  0.00          "
            " C 0\nHETATM   36 C14  UNL     1       8.428  16.119  52.462  1.0"
            "0  0.00           C 0\nHETATM   37 H23  UNL     1      10.258  15"
            ".352  51.489  1.00  0.00           H 0\nHETATM   38 H24  UNL     "
            "1      10.430  15.760  53.258  1.00  0.00           H 0\nHETATM  "
            " 39 H25  UNL     1       7.681  15.552  51.865  1.00  0.00       "
            "    H 0\nHETATM   40 H26  UNL     1       8.546  17.170  52.123  "
            "1.00  0.00           H 0\nCONECT    1    2               \nCONECT"
            "    2    3               \nCONECT    2    4               \nCONEC"
            "T    1    5               \nCONECT    1    6               \nCONE"
            "CT    2    7               \nCONECT    3    8               \nCON"
            "ECT    3    9               \nCONECT    4   10               \nCO"
            "NECT    4   11               \nCONECT   12   13               \nC"
            "ONECT   13   14               \nCONECT   13   15               \n"
            "CONECT   12   16               \nCONECT   12   17               "
            "\nCONECT   13   18               \nCONECT   14   19              "
            " \nCONECT   14   20               \nCONECT   15   21             "
            "  \nCONECT   15   22               \nCONECT   23   24            "
            "   \nCONECT   23   25               \nCONECT   23   26           "
            "    \nCONECT   24   27               \nCONECT   24   28          "
            "     \nCONECT   29   30               \nCONECT   29   31         "
            "      \nCONECT   29   32               \nCONECT   30   33        "
            "       \nCONECT   30   34               \nCONECT   35   36       "
            "        \nCONECT   35   37               \nCONECT   35   38      "
            "         \nCONECT   36   39               \nCONECT   36   40     "
            "          \nCONECT    1   23               \nCONECT   12   24    "
            "           \nCONECT    4   29               \nCONECT   15   30   "
            "            \nCONECT    3   35               \nCONECT   14   36  "
            "             \nEND\n"
        ),
        periodic_info=construction_result.get_periodic_info(),  # type: ignore[attr-defined]
    )


@pytest.fixture(
    scope="session",
    params=(
        lambda: CaseData(
            molecule=stk.BuildingBlock("BrCCBr", [stk.BromoFactory()]),
            writer=stk.PdbWriter(),
            string=(
                "HETATM    1 Br1  UNL     1      -1.379   1.585   0.277  1.00 "
                " 0.00          Br 0\nHETATM    2 C1   UNL     1      -0.729  "
                "-0.231   0.139  1.00  0.00           C 0\nHETATM    3 C2   UN"
                "L     1       0.710  -0.142  -0.348  1.00  0.00           C 0"
                "\nHETATM    4 Br2  UNL     1       1.634   0.894   1.020  1.0"
                "0  0.00          Br 0\nHETATM    5 H1   UNL     1      -1.386"
                "  -0.820  -0.518  1.00  0.00           H 0\nHETATM    6 H2   "
                "UNL     1      -0.765  -0.635   1.186  1.00  0.00           H"
                " 0\nHETATM    7 H3   UNL     1       0.747   0.491  -1.260  1"
                ".00  0.00           H 0\nHETATM    8 H4   UNL     1       1.1"
                "69  -1.143  -0.496  1.00  0.00           H 0\nCONECT    1    "
                "2               \nCONECT    2    3               \nCONECT    "
                "3    4               \nCONECT    2    5               \nCONEC"
                "T    2    6               \nCONECT    3    7               \n"
                "CONECT    3    8               \nEND\n"
            ),
            periodic_info=None,
        ),
        _get_cof_case,
    ),
)
def case_data(request) -> CaseData:
    return request.param()
