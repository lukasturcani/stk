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
            "CRYST1   36.393   36.392  105.054  90.00  90.00  60.00 P 1\n"
            "HETATM    1 C1   UNL     1      19.448  11.228  52.527  1.00  "
            "0.00           C 0\nHETATM    2 C2   UNL     1      18.203  10.51"
            "1  52.946  1.00  0.00           C 0\nHETATM    3 C3   UNL     1  "
            "    16.948  11.228  52.527  1.00  0.00           C 0\nHETATM    4"
            " C4   UNL     1      18.193   9.060  52.527  1.00  0.00          "
            " C 0\nHETATM    5 H1   UNL     1      20.338  10.694  52.917  1.0"
            "0  0.00           H 0\nHETATM    6 H2   UNL     1      19.417  12"
            ".237  52.984  1.00  0.00           H 0\nHETATM    7 H3   UNL     "
            "1      18.197  10.515  54.061  1.00  0.00           H 0\nHETATM  "
            "  8 H4   UNL     1      16.857  12.176  53.068  1.00  0.00       "
            "    H 0\nHETATM    9 H5   UNL     1      16.094  10.585  52.847  "
            "1.00  0.00           H 0\nHETATM   10 H6   UNL     1      18.174 "
            "  8.915  51.435  1.00  0.00           H 0\nHETATM   11 H7   UNL  "
            "   1      17.282   8.573  52.930  1.00  0.00           H 0\nHETAT"
            "M   12 C5   UNL     1      35.141  20.288  52.527  1.00  0.00    "
            "       C 0\nHETATM   13 C6   UNL     1      36.386  21.005  52.94"
            "6  1.00  0.00           C 0\nHETATM   14 C7   UNL     1      37.6"
            "41  20.288  52.527  1.00  0.00           C 0\nHETATM   15 C8   UN"
            "L     1      36.396  22.456  52.527  1.00  0.00           C 0\nHE"
            "TATM   16 H8   UNL     1      34.251  20.822  52.917  1.00  0.00 "
            "          H 0\nHETATM   17 H9   UNL     1      35.172  19.279  52"
            ".984  1.00  0.00           H 0\nHETATM   18 H10  UNL     1      3"
            "6.392  21.002  54.061  1.00  0.00           H 0\nHETATM   19 H11 "
            " UNL     1      37.732  19.340  53.068  1.00  0.00           H 0"
            "\nHETATM   20 H12  UNL     1      38.495  20.931  52.847  1.00  0"
            ".00           H 0\nHETATM   21 H13  UNL     1      36.415  22.601"
            "  51.435  1.00  0.00           H 0\nHETATM   22 H14  UNL     1   "
            "   37.307  22.943  52.930  1.00  0.00           H 0\nHETATM   23 "
            "C9   UNL     1      26.630  15.374  52.527  1.00  0.00           "
            "C 0\nHETATM   24 C10  UNL     1      27.959  16.142  52.527  1.00"
            "  0.00           C 0\nHETATM   25 H15  UNL     1      26.679  14."
            "427  51.988  1.00  0.00           H 0\nHETATM   26 H16  UNL     1"
            "      26.257  15.252  53.562  1.00  0.00           H 0\nHETATM   "
            "27 H17  UNL     1      28.198  16.498  51.505  1.00  0.00        "
            "   H 0\nHETATM   28 H18  UNL     1      28.771  15.537  52.967  1"
            ".00  0.00           H 0\nHETATM   29 C11  UNL     1      18.196  "
            " 0.767  52.527  1.00  0.00           C 0\nHETATM   30 C12  UNL   "
            "  1      18.196  -0.767  52.527  1.00  0.00           C 0\nHETATM"
            "   31 H19  UNL     1      17.242   1.198  52.221  1.00  0.00     "
            "      H 0\nHETATM   32 H20  UNL     1      18.538   1.151  53.508"
            "  1.00  0.00           H 0\nHETATM   33 H21  UNL     1      18.11"
            "9  -1.152  51.490  1.00  0.00           H 0\nHETATM   34 H22  UNL"
            "     1      17.410  -1.168  53.190  1.00  0.00           H 0\nHET"
            "ATM   35 C13  UNL     1       9.763  15.374  52.527  1.00  0.00  "
            "         C 0\nHETATM   36 C14  UNL     1       8.434  16.142  52."
            "527  1.00  0.00           C 0\nHETATM   37 H23  UNL     1      10"
            ".263  15.379  51.557  1.00  0.00           H 0\nHETATM   38 H24  "
            "UNL     1      10.418  15.741  53.341  1.00  0.00           H 0\n"
            "HETATM   39 H25  UNL     1       7.684  15.612  51.906  1.00  0.0"
            "0           H 0\nHETATM   40 H26  UNL     1       8.577  17.192  "
            "52.218  1.00  0.00           H 0\nCONECT    1    2               "
            "\nCONECT    2    3               \nCONECT    2    4              "
            " \nCONECT    1    5               \nCONECT    1    6             "
            "  \nCONECT    2    7               \nCONECT    3    8            "
            "   \nCONECT    3    9               \nCONECT    4   10           "
            "    \nCONECT    4   11               \nCONECT   12   13          "
            "     \nCONECT   13   14               \nCONECT   13   15         "
            "      \nCONECT   12   16               \nCONECT   12   17        "
            "       \nCONECT   13   18               \nCONECT   14   19       "
            "        \nCONECT   14   20               \nCONECT   15   21      "
            "         \nCONECT   15   22               \nCONECT   23   24     "
            "          \nCONECT   23   25               \nCONECT   23   26    "
            "           \nCONECT   24   27               \nCONECT   24   28   "
            "            \nCONECT   29   30               \nCONECT   29   31  "
            "             \nCONECT   29   32               \nCONECT   30   33 "
            "              \nCONECT   30   34               \nCONECT   35   36"
            "               \nCONECT   35   37               \nCONECT   35   3"
            "8               \nCONECT   36   39               \nCONECT   36   "
            "40               \nCONECT    1   23               \nCONECT   12  "
            " 24               \nCONECT    4   29               \nCONECT   15 "
            "  30               \nCONECT    3   35               \nCONECT   14"
            "   36               \nEND\n"
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
                "HETATM    1 Br1  UNL     1      -1.430   1.550   0.363  1.00 "
                " 0.00          Br 0\nHETATM    2 C1   UNL     1      -0.738  "
                "-0.259   0.116  1.00  0.00           C 0\nHETATM    3 C2   UN"
                "L     1       0.718  -0.111  -0.345  1.00  0.00           C 0"
                "\nHETATM    4 Br2  UNL     1       1.632   0.852   1.090  1.0"
                "0  0.00          Br 0\nHETATM    5 H1   UNL     1      -1.364"
                "  -0.797  -0.599  1.00  0.00           H 0\nHETATM    6 H2   "
                "UNL     1      -0.778  -0.704   1.129  1.00  0.00           H"
                " 0\nHETATM    7 H3   UNL     1       0.778   0.565  -1.221  1"
                ".00  0.00           H 0\nHETATM    8 H4   UNL     1       1.1"
                "81  -1.096  -0.534  1.00  0.00           H 0\nCONECT    1    "
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
