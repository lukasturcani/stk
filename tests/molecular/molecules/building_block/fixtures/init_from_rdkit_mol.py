import pytest
import rdkit.Chem.AllChem as rdkit
import stk

from ..case_data import CaseData


rdkit_molecule = rdkit.MolFromSmiles('Br[C+2][C+2]Br')
rdkit.EmbedMolecule(rdkit_molecule, rdkit.ETKDGv2())


@pytest.fixture(
    params=(
        CaseData(
            building_block=stk.BuildingBlock.init_from_rdkit_mol(
                molecule=rdkit_molecule,
            ),
            functional_groups=(),
            core_atom_ids=(0, 1, 2, 3),
            placer_ids=(0, 1, 2, 3),
        ),
        CaseData(
            building_block=stk.BuildingBlock.init_from_rdkit_mol(
                molecule=rdkit_molecule,
                functional_groups=[stk.BromoFactory()],
            ),
            functional_groups=(
                stk.Bromo(
                    bromine=stk.Br(0),
                    atom=stk.C(1),
                    bonders=(stk.C(1), ),
                    deleters=(stk.Br(0), ),
                ),
                stk.Bromo(
                    bromine=stk.Br(3),
                    atom=stk.C(2),
                    bonders=(stk.C(2), ),
                    deleters=(stk.Br(3), ),
                ),
            ),
            core_atom_ids=(1, 2),
            placer_ids=(1, 2),
        ),
        CaseData(
            building_block=stk.BuildingBlock.init_from_rdkit_mol(
                molecule=rdkit_molecule,
                placer_ids=(1, 2),
            ),
            functional_groups=(),
            core_atom_ids=(0, 1, 2, 3),
            placer_ids=(1, 2),
        ),
        CaseData(
            building_block=stk.BuildingBlock.init_from_rdkit_mol(
                molecule=rdkit_molecule,
                functional_groups=[stk.BromoFactory()],
                placer_ids=(0, 3),
            ),
            functional_groups=(
                stk.Bromo(
                    bromine=stk.Br(0),
                    atom=stk.C(1),
                    bonders=(stk.C(1), ),
                    deleters=(stk.Br(0), ),
                ),
                stk.Bromo(
                    bromine=stk.Br(3),
                    atom=stk.C(2),
                    bonders=(stk.C(2), ),
                    deleters=(stk.Br(3), ),
                ),
            ),
            core_atom_ids=(1, 2),
            placer_ids=(0, 3),
        ),
        CaseData(
            building_block=stk.BuildingBlock.init_from_rdkit_mol(
                molecule=rdkit_molecule,
                functional_groups=[stk.IodoFactory()],
            ),
            functional_groups=(),
            core_atom_ids=(0, 1, 2, 3),
            placer_ids=(0, 1, 2, 3),
        ),
    ),
)
def init_from_rdkit_mol(request):
    return request.param
