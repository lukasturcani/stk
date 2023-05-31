from dataclasses import dataclass
from typing import Callable

import pymongo
import pytest
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
        :class:`.ValueMongoDb` instance.

    molecule : :class:`.Molecule`
        The molecule to test.

    value : :class:`object`
        The value to put into the database.

    """

    get_database: Callable[[pymongo.MongoClient], stk.ValueMongoDb]
    molecule: stk.Molecule
    value: object


@pytest.fixture(
    params=(
        lambda: CaseDataData(
            get_database=lambda mongo_client: stk.ValueMongoDb(
                mongo_client=mongo_client,
                collection="values",
                database="_stk_test_database_for_testing",
                put_lru_cache_size=0,
                get_lru_cache_size=0,
            ),
            molecule=stk.BuildingBlock("BrCCBr"),
            value=12,
        ),
        lambda: CaseDataData(
            get_database=lambda mongo_client: stk.ValueMongoDb(
                mongo_client=mongo_client,
                collection="values",
                database="_stk_test_database_for_testing",
                put_lru_cache_size=128,
                get_lru_cache_size=128,
            ),
            molecule=stk.BuildingBlock("BrCCBr"),
            value=12,
        ),
    ),
)
def mongo_db(
    request,
    mongo_client: pymongo.MongoClient,
) -> CaseData:
    data = request.param()
    return CaseData(
        database=data.get_database(mongo_client),
        molecule=data.molecule,
        value=data.value,
    )
