from dataclasses import dataclass
from typing import Callable

import pymongo
import pytest
import rdkit.Chem.AllChem as rdkit
import stk

from ..case_data import CaseData


@dataclass(frozen=True)
class CaseDataData:
    """
    Data used to create a :class:`.CaseData` instance.

    Attributes
    ----------
    get_database : :class:`callable`
        Creates the database to test. Takes a
        :class:`pymongo.MongoClient` as input and returns a
        :class:`.MoleculeMongoDb` instance.

    molecule : :class:`.Molecule`
        The molecule to put and get from the database being tested.

    key : :class:`object`
        The key used to retrieve the :attr:`.molecule` from the
        database.

    """

    get_database: Callable[[pymongo.MongoClient], stk.ValueMongoDb]
    molecule: stk.Molecule
    key: object


@pytest.fixture(
    params=(
        lambda: CaseDataData(
            get_database=lambda mongo_client: (
                stk.ConstructedMoleculeMongoDb(
                    mongo_client=mongo_client,
                    database="_stk_test_database_for_testing",
                    put_lru_cache_size=0,
                    get_lru_cache_size=0,
                )
            ),
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles="BrCCBr",
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    repeating_unit="A",
                    num_repeating_units=2,
                ),
            ),
            key={
                "InChIKey": rdkit.MolToInchiKey(
                    rdkit.MolFromSmiles(SMILES="BrCCCCBr")
                ),
            },
        ),
        lambda: CaseDataData(
            get_database=lambda mongo_client: (
                stk.ConstructedMoleculeMongoDb(
                    mongo_client=mongo_client,
                    database="_stk_test_database_for_testing",
                    put_lru_cache_size=128,
                    get_lru_cache_size=128,
                )
            ),
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles="BrCCBr",
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    repeating_unit="A",
                    num_repeating_units=2,
                ),
            ),
            key={
                "InChIKey": rdkit.MolToInchiKey(
                    rdkit.MolFromSmiles(SMILES="BrCCCCBr")
                ),
            },
        ),
        lambda: CaseDataData(
            get_database=lambda mongo_client: (
                stk.ConstructedMoleculeMongoDb(
                    mongo_client=mongo_client,
                    database="_stk_test_database_for_testing",
                    jsonizer=stk.ConstructedMoleculeJsonizer(
                        key_makers=(
                            stk.MoleculeKeyMaker(
                                key_name="SMILES",
                                get_key=lambda molecule: (
                                    rdkit.MolToSmiles(
                                        mol=molecule.to_rdkit_mol(),
                                    )
                                ),
                            ),
                        ),
                    ),
                    put_lru_cache_size=0,
                    get_lru_cache_size=0,
                )
            ),
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles="Br[C+2][C+2]Br",
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    repeating_unit="A",
                    num_repeating_units=2,
                ),
            ),
            key={"SMILES": "Br[C+2][C+2][C+2][C+2]Br"},
        ),
        lambda: CaseDataData(
            get_database=lambda mongo_client: (
                stk.ConstructedMoleculeMongoDb(
                    mongo_client=mongo_client,
                    database="_stk_test_database_for_testing",
                    jsonizer=stk.ConstructedMoleculeJsonizer(
                        key_makers=(
                            stk.MoleculeKeyMaker(
                                key_name="SMILES",
                                get_key=lambda molecule: (
                                    rdkit.MolToSmiles(
                                        mol=molecule.to_rdkit_mol(),
                                    )
                                ),
                            ),
                        ),
                    ),
                    put_lru_cache_size=128,
                    get_lru_cache_size=128,
                )
            ),
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles="Br[C+2][C+2]Br",
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    repeating_unit="A",
                    num_repeating_units=2,
                ),
            ),
            key={"SMILES": "Br[C+2][C+2][C+2][C+2]Br"},
        ),
    ),
)
def constructed_molecule_mongo_db(
    request,
    mongo_client: pymongo.MongoClient,
) -> CaseData:
    data = request.param()
    return CaseData(
        database=data.get_database(mongo_client),
        molecule=data.molecule,
        key=data.key,
    )
