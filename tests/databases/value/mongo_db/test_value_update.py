import stk

from ...utilities import (
    DatabaseEntry,
    DatabaseState,
    assert_database_state,
)
from .utilities import get_database_state


def test_update_1(mongo_client):
    """
    Test that existing entries are updated.

    """

    collection = "_test_update_1"
    database_name = "_test_update_1"
    mongo_client.drop_database(database_name)

    database = stk.ValueMongoDb(
        mongo_client=mongo_client,
        collection=collection,
        database=database_name,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
    )

    molecule = stk.BuildingBlock("CCC")

    database.put(molecule, 12)
    assert_database_state(
        state1=get_database_state(database),
        state2=DatabaseState(
            {
                DatabaseEntry(
                    InChIKey=stk.InchiKey().get_key(molecule),
                    v=12,
                ): 1,
            }
        ),
    )

    database.put(molecule, 43)
    assert_database_state(
        state1=get_database_state(database),
        state2=DatabaseState(
            {
                DatabaseEntry(
                    InChIKey=stk.InchiKey().get_key(molecule),
                    v=43,
                ): 1,
            }
        ),
    )


def test_update_2(mongo_client):
    """
    Test that existing entries are updated.

    In this test, you first create two separate entries, using
    different molecule keys. You then update both at the same time,
    with a database which uses both molecule keys.

    """

    collection = "_test_update_2"
    database_name = "_test_update_2"
    mongo_client.drop_database(database_name)

    database1 = stk.ValueMongoDb(
        mongo_client=mongo_client,
        collection=collection,
        database=database_name,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        key_makers=(stk.InchiKey(),),
    )
    database2 = stk.ValueMongoDb(
        mongo_client=mongo_client,
        collection=collection,
        database=database_name,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        key_makers=(stk.Smiles(),),
    )
    database3 = stk.ValueMongoDb(
        mongo_client=mongo_client,
        collection=collection,
        database=database_name,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        key_makers=(
            stk.InchiKey(),
            stk.Smiles(),
        ),
    )

    molecule = stk.BuildingBlock("CCC")

    database1.put(molecule, 12)
    assert_database_state(
        state1=get_database_state(database1),
        state2=DatabaseState(
            {
                DatabaseEntry(
                    InChIKey=stk.InchiKey().get_key(molecule),
                    v=12,
                ): 1,
            }
        ),
    )

    # Should add another entry, as a different key maker is used.
    database2.put(molecule, 32)
    assert_database_state(
        state1=get_database_state(database1),
        state2=DatabaseState(
            {
                DatabaseEntry(
                    InChIKey=stk.InchiKey().get_key(molecule),
                    v=12,
                ): 1,
                DatabaseEntry(
                    SMILES=stk.Smiles().get_key(molecule),
                    v=32,
                ): 1,
            }
        ),
    )

    # Should update both entries as both key makers are used.
    database3.put(molecule, 56)
    assert_database_state(
        state1=get_database_state(database1),
        state2=DatabaseState(
            {
                DatabaseEntry(
                    InChIKey=stk.InchiKey().get_key(molecule),
                    SMILES=stk.Smiles().get_key(molecule),
                    v=56,
                ): 2,
            }
        ),
    )


def test_update_3(mongo_client):
    """
    Test that existing entries are updated.

    In this test, you first create one entry with two keys. Then
    update the entry with databases, each using 1 different key.
    No duplicate entries should be made in the database this way.

    """

    collection = "_test_update_3"
    database_name = "_test_update_3"
    mongo_client.drop_database(database_name)

    database1 = stk.ValueMongoDb(
        mongo_client=mongo_client,
        collection=collection,
        database=database_name,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        key_makers=(
            stk.InchiKey(),
            stk.Smiles(),
        ),
    )
    database2 = stk.ValueMongoDb(
        mongo_client=mongo_client,
        collection=collection,
        database=database_name,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        key_makers=(stk.InchiKey(),),
    )
    database3 = stk.ValueMongoDb(
        mongo_client=mongo_client,
        collection=collection,
        database=database_name,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        key_makers=(stk.Smiles(),),
    )

    molecule = stk.BuildingBlock("CCC")

    database1.put(molecule, 12)
    assert_database_state(
        state1=get_database_state(database1),
        state2=DatabaseState(
            {
                DatabaseEntry(
                    InChIKey=stk.InchiKey().get_key(molecule),
                    SMILES=stk.Smiles().get_key(molecule),
                    v=12,
                ): 1
            }
        ),
    )

    # Should update the entry.
    database2.put(molecule, 32)
    assert_database_state(
        state1=get_database_state(database1),
        state2=DatabaseState(
            {
                DatabaseEntry(
                    InChIKey=stk.InchiKey().get_key(molecule),
                    SMILES=stk.Smiles().get_key(molecule),
                    v=32,
                ): 1,
            }
        ),
    )

    # Should also update the entry.
    database3.put(molecule, 62)
    assert_database_state(
        state1=get_database_state(database1),
        state2=DatabaseState(
            {
                DatabaseEntry(
                    InChIKey=stk.InchiKey().get_key(molecule),
                    SMILES=stk.Smiles().get_key(molecule),
                    v=62,
                ): 1,
            }
        ),
    )
