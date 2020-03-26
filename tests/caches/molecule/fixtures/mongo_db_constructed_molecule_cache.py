import pytest
import stk
import rdkit.Chem.AllChem as rdkit

from .utilities import MockMongoClient
from ..case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            cache=stk.MongoDbConstructedMoleculeCache(
                mongo_client=MockMongoClient(),
                lru_cache_size=0,
            ),
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
            key={
                'InChIKey':
                    rdkit.MolToInchiKey(rdkit.MolFromSmiles(
                        SMILES='BrCCCCBr'
                    )),
            },
        ),
        CaseData(
            cache=stk.MongoDbConstructedMoleculeCache(
                mongo_client=MockMongoClient(),
                jsonizer=stk.ConstructedMoleculeJsonizer(
                    key_makers=(
                        stk.MoleculeKeyMaker(
                            smiles='SMILES',
                            get_key=lambda molecule: rdkit.MolToSmiles(
                                mol=molecule.to_rdkit_mol(),
                            )
                        ),
                    ),
                ),
                lru_cache_size=0,
            ),
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='Br[C+2][C+2]Br',
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    repeating_unit='A',
                    num_repeating_units=2,
                ),
            ),
            key={
                'SMILES':
                    rdkit.MolToInchiKey(rdkit.MolFromSmiles(
                        SMILES='Br[C+2][C+2][C+2][C+2]Br'
                    )),
            },
        ),
    ),
)
def mongo_db_constructed_molecule_cache(request):
    return request.param
