import pytest
import stk
import rdkit.Chem.AllChem as rdkit
import pymongo

from ..case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            database=stk.ConstructedMoleculeMongoDb(
                mongo_client=pymongo.MongoClient(),
                database='_stk_test_database_for_testing',
                put_lru_cache_size=0,
                get_lru_cache_size=0,
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
            database=stk.ConstructedMoleculeMongoDb(
                mongo_client=pymongo.MongoClient(),
                database='_stk_test_database_for_testing',
                put_lru_cache_size=128,
                get_lru_cache_size=128,
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
            database=stk.ConstructedMoleculeMongoDb(
                mongo_client=pymongo.MongoClient(),
                database='_stk_test_database_for_testing',
                jsonizer=stk.ConstructedMoleculeJsonizer(
                    key_makers=(
                        stk.MoleculeKeyMaker(
                            key_name='SMILES',
                            get_key=lambda molecule: rdkit.MolToSmiles(
                                mol=molecule.to_rdkit_mol(),
                            )
                        ),
                    ),
                ),
                put_lru_cache_size=0,
                get_lru_cache_size=0,
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
            key={'SMILES': 'Br[C+2][C+2][C+2][C+2]Br'},
        ),
        CaseData(
            database=stk.ConstructedMoleculeMongoDb(
                mongo_client=pymongo.MongoClient(),
                database='_stk_test_database_for_testing',
                jsonizer=stk.ConstructedMoleculeJsonizer(
                    key_makers=(
                        stk.MoleculeKeyMaker(
                            key_name='SMILES',
                            get_key=lambda molecule: rdkit.MolToSmiles(
                                mol=molecule.to_rdkit_mol(),
                            )
                        ),
                    ),
                ),
                put_lru_cache_size=128,
                get_lru_cache_size=128,
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
            key={'SMILES': 'Br[C+2][C+2][C+2][C+2]Br'},
        ),
    ),
)
def constructed_molecule_mongo_db(request):
    return request.param
