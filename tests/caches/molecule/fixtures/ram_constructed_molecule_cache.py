import pytest
import rdkit.Chem.AllChem as rdkit
import stk

from ..case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            cache=stk.RamConstructedMoleculeCache(),
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='BrCCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    repeating_unit='A',
                    num_repeating_units=2,
                ),
            ),
            key=rdkit.MolToInchiKey(rdkit.MolFromSmiles('BrCCCCBr')),
        ),
    ),
)
def ram_constructed_molecule_cache(request):
    return request.param
