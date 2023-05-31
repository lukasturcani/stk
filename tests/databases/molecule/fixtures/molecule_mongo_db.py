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
    Data used to create a :class:`.CaseDataData` instance.

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
            get_database=lambda mongo_client: stk.MoleculeMongoDb(
                mongo_client=mongo_client,
                database="_stk_test_database_for_testing",
                put_lru_cache_size=0,
                get_lru_cache_size=0,
            ),
            molecule=stk.BuildingBlock("BrCCBr"),
            key={
                "InChIKey": rdkit.MolToInchiKey(rdkit.MolFromSmiles("BrCCBr"))
            },
        ),
        lambda: CaseDataData(
            get_database=lambda mongo_client: stk.MoleculeMongoDb(
                mongo_client=mongo_client,
                database="_stk_test_database_for_testing",
                put_lru_cache_size=0,
                get_lru_cache_size=0,
                jsonizer=stk.MoleculeJsonizer(
                    key_makers=(
                        stk.MoleculeKeyMaker(
                            key_name="SMILES",
                            get_key=lambda molecule: rdkit.MolToSmiles(
                                molecule.to_rdkit_mol()
                            ),
                        ),
                    ),
                ),
            ),
            molecule=stk.BuildingBlock("BrBr"),
            key={"SMILES": "BrBr"},
        ),
        lambda: CaseDataData(
            get_database=lambda mongo_client: stk.MoleculeMongoDb(
                mongo_client=mongo_client,
                database="_stk_test_database_for_testing",
                put_lru_cache_size=128,
                get_lru_cache_size=128,
                jsonizer=stk.MoleculeJsonizer(
                    key_makers=(
                        stk.MoleculeKeyMaker(
                            key_name="SMILES",
                            get_key=lambda molecule: rdkit.MolToSmiles(
                                molecule.to_rdkit_mol()
                            ),
                        ),
                    ),
                ),
            ),
            molecule=stk.BuildingBlock("BrBr"),
            key={"SMILES": "BrBr"},
        ),
    ),
)
def molecule_mongo_db(
    request,
    mongo_client: pymongo.MongoClient,
) -> CaseData:
    data = request.param()
    return CaseData(
        database=data.get_database(mongo_client),
        molecule=data.molecule,
        key=data.key,
    )
