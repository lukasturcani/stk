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
    ),
)
def case_data(request):

    return request.param
