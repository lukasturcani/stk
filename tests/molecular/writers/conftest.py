import pytest
import stk

from .case_data import CaseData

bb1 = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])
bb2 = stk.BuildingBlock('BrCC(CBr)CBr', [stk.BromoFactory()])
topology_graph = stk.cof.Honeycomb(
    building_blocks=(bb1, bb2),
    lattice_size=(1, 1, 1),
    periodic=True
)
cof = stk.ConstructedMolecule(topology_graph)


@pytest.fixture(
    params=(
        CaseData(
            molecule=bb1,
            writer=stk.TurbomoleWriter(),
            string=(
                '$coord angs\n -1.4238 1.5615 0.3223 Br\n -0.7405 -0.2'
                '573 0.128 C\n 0.7148 -0.1157 -0.3383 C\n 1.6267 0.889'
                '6 1.0687 Br\n -1.3518 -0.8075 -0.5939 H\n -0.7769 -0.'
                '6964 1.144 H\n 0.7695 0.528 -1.2387 H\n 1.1821 -1.102'
                '2 -0.4922 H\n$end\n'
            ),
            periodic_cell=None,
        ),
        CaseData(
            molecule=bb1,
            writer=stk.PdbWriter(),
            string=(
                'HETATM    1 Br1  UNL     1      -1.424   1.561   0.32'
                '2  1.00  0.00          Br 0\nHETATM    2 C1   UNL    '
                ' 1      -0.741  -0.257   0.128  1.00  0.00           '
                'C 0\nHETATM    3 C2   UNL     1       0.715  -0.116  '
                '-0.338  1.00  0.00           C 0\nHETATM    4 Br2  UN'
                'L     1       1.627   0.890   1.069  1.00  0.00      '
                '    Br 0\nHETATM    5 H1   UNL     1      -1.352  -0.'
                '807  -0.594  1.00  0.00           H 0\nHETATM    6 H2'
                '   UNL     1      -0.777  -0.696   1.144  1.00  0.00 '
                '          H 0\nHETATM    7 H3   UNL     1       0.769'
                '   0.528  -1.239  1.00  0.00           H 0\nHETATM   '
                ' 8 H4   UNL     1       1.182  -1.102  -0.492  1.00  '
                '0.00           H 0\nCONECT    1    2               \n'
                'CONECT    2    3               \nCONECT    3    4    '
                '           \nCONECT    2    5               \nCONECT '
                '   2    6               \nCONECT    3    7           '
                '    \nCONECT    3    8               \nEND\n'
            ),
            periodic_cell=None,
        ),

        CaseData(
            molecule=cof,
            writer=stk.TurbomoleWriter(),
            string=(
                '$periodic 3\n$cell angs\n   36.432   36.431  105.166 '
                ' 90.00  90.00  60.00\n$coord angs\n 19.4675 11.2392 5'
                '2.5831 C\n 18.2205 10.5214 53.0065 C\n 16.9674 11.238'
                '9 52.5831 C\n 18.2126 9.0717 52.5831 C\n 20.3559 10.6'
                '937 52.9602 H\n 19.4494 12.2444 53.05 H\n 18.2209 10.'
                '5188 54.1253 H\n 16.8933 12.1964 53.1311 H\n 16.1039 '
                '10.6163 52.92 H\n 18.198 8.9414 51.4898 H\n 17.295 8.'
                '586 52.9733 H\n 35.18 20.3106 52.5831 C\n 36.427 21.0'
                '285 53.0065 C\n 37.6801 20.3109 52.5831 C\n 36.4349 2'
                '2.4781 52.5831 C\n 34.2916 20.8561 52.9602 H\n 35.198'
                '1 19.3054 53.05 H\n 36.4266 21.031 54.1253 H\n 37.754'
                '2 19.3534 53.1311 H\n 38.5436 20.9335 52.92 H\n 36.44'
                '95 22.6084 51.4898 H\n 37.3525 22.9638 52.9733 H\n 26'
                '.6592 15.3912 52.5831 C\n 27.9883 16.1586 52.5831 C\n'
                ' 26.5553 14.749 53.4629 H\n 25.811 16.0929 52.4614 H'
                '\n 28.838 15.4564 52.469 H\n 28.081 16.7958 53.4778 H'
                '\n 18.2158 0.7674 52.5831 C\n 18.2158 -0.7674 52.5831'
                ' C\n 17.7115 1.1784 53.4629 H\n 19.2476 1.1511 52.461'
                '4 H\n 17.1829 -1.1522 52.469 H\n 18.7214 -1.1662 53.4'
                '778 H\n 9.7725 15.3912 52.5831 C\n 8.4433 16.1586 52.'
                '5831 C\n 10.3806 15.6224 53.4629 H\n 9.5889 14.3059 5'
                '2.4614 H\n 8.6266 17.2456 52.469 H\n 7.8451 15.9202 5'
                '3.4778 H\n$end\n'
            ),
            periodic_cell=topology_graph.get_periodic_cell(),
        ),
        CaseData(
            molecule=cof,
            writer=stk.PdbWriter(),
            string=(
                'CRYST1   36.432   36.431  105.166  90.00  90.00  60.0'
                '0 P 1\nHETATM    1 C1   UNL     1      19.467  11.239'
                '  52.583  1.00  0.00           C 0\nHETATM    2 C2   '
                'UNL     1      18.220  10.521  53.006  1.00  0.00    '
                '       C 0\nHETATM    3 C3   UNL     1      16.967  1'
                '1.239  52.583  1.00  0.00           C 0\nHETATM    4 '
                'C4   UNL     1      18.213   9.072  52.583  1.00  0.0'
                '0           C 0\nHETATM    5 H1   UNL     1      20.3'
                '56  10.694  52.960  1.00  0.00           H 0\nHETATM '
                '   6 H2   UNL     1      19.449  12.244  53.050  1.00'
                '  0.00           H 0\nHETATM    7 H3   UNL     1     '
                ' 18.221  10.519  54.125  1.00  0.00           H 0\nHE'
                'TATM    8 H4   UNL     1      16.893  12.196  53.131 '
                ' 1.00  0.00           H 0\nHETATM    9 H5   UNL     1'
                '      16.104  10.616  52.920  1.00  0.00           H '
                '0\nHETATM   10 H6   UNL     1      18.198   8.941  51'
                '.490  1.00  0.00           H 0\nHETATM   11 H7   UNL '
                '    1      17.295   8.586  52.973  1.00  0.00        '
                '   H 0\nHETATM   12 C5   UNL     1      35.180  20.31'
                '1  52.583  1.00  0.00           C 0\nHETATM   13 C6  '
                ' UNL     1      36.427  21.028  53.006  1.00  0.00   '
                '        C 0\nHETATM   14 C7   UNL     1      37.680  '
                '20.311  52.583  1.00  0.00           C 0\nHETATM   15'
                ' C8   UNL     1      36.435  22.478  52.583  1.00  0.'
                '00           C 0\nHETATM   16 H8   UNL     1      34.'
                '292  20.856  52.960  1.00  0.00           H 0\nHETATM'
                '   17 H9   UNL     1      35.198  19.305  53.050  1.0'
                '0  0.00           H 0\nHETATM   18 H10  UNL     1    '
                '  36.427  21.031  54.125  1.00  0.00           H 0\nH'
                'ETATM   19 H11  UNL     1      37.754  19.353  53.131'
                '  1.00  0.00           H 0\nHETATM   20 H12  UNL     '
                '1      38.544  20.934  52.920  1.00  0.00           H'
                ' 0\nHETATM   21 H13  UNL     1      36.450  22.608  5'
                '1.490  1.00  0.00           H 0\nHETATM   22 H14  UNL'
                '     1      37.353  22.964  52.973  1.00  0.00       '
                '    H 0\nHETATM   23 C9   UNL     1      26.659  15.3'
                '91  52.583  1.00  0.00           C 0\nHETATM   24 C10'
                '  UNL     1      27.988  16.159  52.583  1.00  0.00  '
                '         C 0\nHETATM   25 H15  UNL     1      26.555 '
                ' 14.749  53.463  1.00  0.00           H 0\nHETATM   2'
                '6 H16  UNL     1      25.811  16.093  52.461  1.00  0'
                '.00           H 0\nHETATM   27 H17  UNL     1      28'
                '.838  15.456  52.469  1.00  0.00           H 0\nHETAT'
                'M   28 H18  UNL     1      28.081  16.796  53.478  1.'
                '00  0.00           H 0\nHETATM   29 C11  UNL     1   '
                '   18.216   0.767  52.583  1.00  0.00           C 0\n'
                'HETATM   30 C12  UNL     1      18.216  -0.767  52.58'
                '3  1.00  0.00           C 0\nHETATM   31 H19  UNL    '
                ' 1      17.712   1.178  53.463  1.00  0.00           '
                'H 0\nHETATM   32 H20  UNL     1      19.248   1.151  '
                '52.461  1.00  0.00           H 0\nHETATM   33 H21  UN'
                'L     1      17.183  -1.152  52.469  1.00  0.00      '
                '     H 0\nHETATM   34 H22  UNL     1      18.721  -1.'
                '166  53.478  1.00  0.00           H 0\nHETATM   35 C1'
                '3  UNL     1       9.772  15.391  52.583  1.00  0.00 '
                '          C 0\nHETATM   36 C14  UNL     1       8.443'
                '  16.159  52.583  1.00  0.00           C 0\nHETATM   '
                '37 H23  UNL     1      10.381  15.622  53.463  1.00  '
                '0.00           H 0\nHETATM   38 H24  UNL     1       '
                '9.589  14.306  52.461  1.00  0.00           H 0\nHETA'
                'TM   39 H25  UNL     1       8.627  17.246  52.469  1'
                '.00  0.00           H 0\nHETATM   40 H26  UNL     1  '
                '     7.845  15.920  53.478  1.00  0.00           H 0'
                '\nCONECT    1    2               \nCONECT    2    3  '
                '             \nCONECT    2    4               \nCONEC'
                'T    1    5               \nCONECT    1    6         '
                '      \nCONECT    2    7               \nCONECT    3 '
                '   8               \nCONECT    3    9               '
                '\nCONECT    4   10               \nCONECT    4   11  '
                '             \nCONECT   12   13               \nCONEC'
                'T   13   14               \nCONECT   13   15         '
                '      \nCONECT   12   16               \nCONECT   12 '
                '  17               \nCONECT   13   18               '
                '\nCONECT   14   19               \nCONECT   14   20  '
                '             \nCONECT   15   21               \nCONEC'
                'T   15   22               \nCONECT   23   24         '
                '      \nCONECT   23   25               \nCONECT   23 '
                '  26               \nCONECT   24   27               '
                '\nCONECT   24   28               \nCONECT   29   30  '
                '             \nCONECT   29   31               \nCONEC'
                'T   29   32               \nCONECT   30   33         '
                '      \nCONECT   30   34               \nCONECT   35 '
                '  36               \nCONECT   35   37               '
                '\nCONECT   35   38               \nCONECT   36   39  '
                '             \nCONECT   36   40               \nCONEC'
                'T    1   23               \nCONECT   12   24         '
                '      \nCONECT    4   29               \nCONECT   15 '
                '  30               \nCONECT    3   35               '
                '\nCONECT   14   36               \nEND\n'
            ),
            periodic_cell=topology_graph.get_periodic_cell(),
        ),
    ),
)
def case_data(request):

    return request.param
