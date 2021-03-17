import pytest
import stk

from ...case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cof.PeriodicKagome(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='BrC1=C(Br)[C+]=N1',
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles=(
                                'Br[C+]1C2(Br)[C+]=N[C+]2[C+](Br)[C+]('
                                'Br)[C+2]1'
                            ),
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    lattice_size=(2, 2, 1),
                ),
            ),
            smiles=(
                '[C+]1=NC2=C1[C+]1[C+2][C+]3C4=C(N=[C+]4)C45[C+]=N[C+]'
                '4[C+]4C6=C([C+]=N6)[C+]6[C+]7[C+2][C+]8C9=C(N=[C+]9)['
                'C+]9[C+2][C+]%10C%11=C([C+]=N%11)[C+]%11[C+]%12[C+2]['
                'C+]%13C%14=C(N=[C+]%14)C%14%15[C+]=N[C+]%14[C+]%14C%1'
                '6=C([C+]=N%16)[C+]%16[C+]%17[C+2][C+]%18C%19=C(N=[C+]'
                '%19)[C+]%19[C+2][C+]%20C%21=C([C+]=N%21)C3%21[C+]=N[C'
                '+]%21[C+]1C1=C(N=[C+]1)C13[C+]=N[C+]1[C+]1C%21=C(N=[C'
                '+]%21)[C+]%20[C+]%20N=[C+]C%19%20C%19=C([C+]=N%19)[C+'
                ']%19[C+2][C+](C%20=C(N=[C+]%20)C8%20[C+]=N[C+]6%20)[C'
                '+](C6=C(N=[C+]6)[C+]4[C+2][C+]5C4=C([C+]=N4)[C+]4[C+2'
                '][C+]2[C+](C2=C(N=[C+]2)[C+]%10[C+]2N=[C+]C92C2=C([C+'
                ']=N2)[C+]2[C+2][C+](C5=C(N=[C+]5)C%185[C+]=N[C+]%165)'
                '[C+](C5=C(N=[C+]5)[C+]%14[C+2][C+]%15C5=C([C+]=N5)[C+'
                ']3[C+2][C+]1C1=C([C+]=N1)C%131[C+]=N[C+]%111)[C+]1N=['
                'C+]C21C1=C7N=[C+]1)[C+]1N=[C+]C41C1=C%12[C+]=N1)[C+]1'
                'N=[C+]C%191C1=C%17N=[C+]1'
            ),
        ),
    ),
)
def cof_periodic_kagome(request):
    return request.param
