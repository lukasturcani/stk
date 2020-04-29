import pytest
import stk
from rdkit.Chem import AllChem as rdkit

from .case_data import CaseData

dative_molecule = rdkit.MolFromSmiles('[Fe+2]<-N')
dative_molecule.AddConformer(rdkit.Conformer(
    dative_molecule.GetNumAtoms())
)


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.BuildingBlock('NCCN'),
            bonds=(
                stk.Bond(stk.N(0), stk.C(1), 1),
                stk.Bond(stk.C(1), stk.C(2), 1),
                stk.Bond(stk.C(2), stk.N(3), 1),
                stk.Bond(stk.N(0), stk.H(4), 1),
                stk.Bond(stk.N(0), stk.H(5), 1),
                stk.Bond(stk.C(1), stk.H(6), 1),
                stk.Bond(stk.C(1), stk.H(7), 1),
                stk.Bond(stk.C(2), stk.H(8), 1),
                stk.Bond(stk.C(2), stk.H(9), 1),
                stk.Bond(stk.N(3), stk.H(10), 1),
                stk.Bond(stk.N(3), stk.H(11), 1),
            ),
        ),
        CaseData(
            molecule=stk.BuildingBlock.init_from_rdkit_mol(
                dative_molecule
            ),
            bonds=(
                stk.Bond(
                    stk.Fe(0, charge=2),
                    stk.N(1),
                    1,
                    is_dative=True
                ),
            ),
        ),
    ),
)
def case_data(request):
    return request.param
