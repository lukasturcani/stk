import numpy as np
import stk

from ...utilities import (
    DatabaseEntry,
    DatabaseState,
    assert_database_state,
)
from .utilities import get_database_state


def to_hashable(json):
    json = dict(json)
    return {
        "m": tuple(tuple(row) for row in json.pop("m")),
        **json,
    }


def test_update_1(mongo_client):
    """
    Test that existing entries are updated.

    """

    database_name = "_test_update_1"
    mongo_client.drop_database(database_name)

    database = stk.MoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
    )
    jsonizer = stk.MoleculeJsonizer()

    molecule = stk.BuildingBlock("CCC").with_canonical_atom_ordering()
    json = jsonizer.to_json(molecule)

    database.put(molecule)
    assert_database_state(
        state1=get_database_state(database),
        state2=DatabaseState(
            {
                DatabaseEntry(**json["molecule"]): 1,
                DatabaseEntry(**to_hashable(json["matrix"])): 1,
            }
        ),
    )

    molecule2 = molecule.with_position_matrix(
        position_matrix=np.zeros((molecule.get_num_atoms(), 3)),
    )
    json2 = jsonizer.to_json(molecule2)

    database.put(molecule2)
    assert_database_state(
        state1=get_database_state(database),
        state2=DatabaseState(
            {
                # Molecule JSON should be unchanged.
                DatabaseEntry(**json["molecule"]): 1,
                # Position matrix should be updated.
                DatabaseEntry(**to_hashable(json2["matrix"])): 1,
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

    database_name = "_test_update_2"
    mongo_client.drop_database(database_name)

    jsonizer1 = stk.MoleculeJsonizer()
    database1 = stk.MoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        jsonizer=jsonizer1,
    )

    jsonizer2 = stk.MoleculeJsonizer(
        key_makers=(stk.Smiles(),),
    )
    database2 = stk.MoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        jsonizer=jsonizer2,
    )

    jsonizer3 = stk.MoleculeJsonizer(
        key_makers=(
            stk.InchiKey(),
            stk.Smiles(),
        ),
    )
    database3 = stk.MoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        jsonizer=jsonizer3,
    )

    molecule1 = stk.BuildingBlock("CCC").with_canonical_atom_ordering()
    json1 = jsonizer1.to_json(molecule1)

    database1.put(molecule1)
    assert_database_state(
        state1=get_database_state(database1),
        state2=DatabaseState(
            {
                DatabaseEntry(**json1["molecule"]): 1,
                DatabaseEntry(**to_hashable(json1["matrix"])): 1,
            }
        ),
    )

    # Should add another entry, as a different key maker is used.
    molecule2 = molecule1.with_position_matrix(
        position_matrix=np.zeros((molecule1.get_num_atoms(), 3)),
    )
    json2 = jsonizer2.to_json(molecule2)

    database2.put(molecule2)
    assert_database_state(
        state1=get_database_state(database1),
        state2=DatabaseState(
            {
                DatabaseEntry(**json1["molecule"]): 1,
                DatabaseEntry(**to_hashable(json1["matrix"])): 1,
                DatabaseEntry(**json2["molecule"]): 1,
                DatabaseEntry(**to_hashable(json2["matrix"])): 1,
            }
        ),
    )

    # Should update both entries as both key makers are used.
    molecule3 = molecule1.with_position_matrix(
        position_matrix=np.ones((molecule1.get_num_atoms(), 3)),
    )
    json3 = jsonizer3.to_json(molecule3)

    database3.put(molecule3)
    assert_database_state(
        state1=get_database_state(database1),
        state2=DatabaseState(
            {
                DatabaseEntry(**json3["molecule"]): 2,
                DatabaseEntry(**to_hashable(json3["matrix"])): 2,
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

    database_name = "_test_update_3"
    mongo_client.drop_database(database_name)

    jsonizer1 = stk.MoleculeJsonizer(
        key_makers=(
            stk.InchiKey(),
            stk.Smiles(),
        ),
    )
    database1 = stk.MoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        jsonizer=jsonizer1,
    )

    jsonizer2 = stk.MoleculeJsonizer()
    database2 = stk.MoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        jsonizer=jsonizer2,
    )

    jsonizer3 = stk.MoleculeJsonizer(
        key_makers=(stk.Smiles(),),
    )
    database3 = stk.MoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        jsonizer=jsonizer3,
    )

    molecule1 = stk.BuildingBlock("CCC").with_canonical_atom_ordering()
    json1 = jsonizer1.to_json(molecule1)

    database1.put(molecule1)
    assert_database_state(
        state1=get_database_state(database1),
        state2=DatabaseState(
            {
                DatabaseEntry(**json1["molecule"]): 1,
                DatabaseEntry(**to_hashable(json1["matrix"])): 1,
            }
        ),
    )

    molecule2 = molecule1.with_position_matrix(
        position_matrix=np.zeros((molecule1.get_num_atoms(), 3)),
    )
    json2 = jsonizer2.to_json(molecule2)
    json2["matrix"] = dict(json1["matrix"])
    json2["matrix"]["m"] = jsonizer2.to_json(molecule2)["matrix"]["m"]

    database2.put(molecule2)
    assert_database_state(
        state1=get_database_state(database1),
        state2=DatabaseState(
            {
                DatabaseEntry(**json1["molecule"]): 1,
                DatabaseEntry(**to_hashable(json2["matrix"])): 1,
            }
        ),
    )

    molecule3 = molecule1.with_position_matrix(
        position_matrix=np.zeros((molecule1.get_num_atoms(), 3)),
    )
    json3 = jsonizer3.to_json(molecule3)
    json3["matrix"] = dict(json1["matrix"])
    json3["matrix"]["m"] = jsonizer3.to_json(molecule3)["matrix"]["m"]

    database3.put(molecule3)
    assert_database_state(
        state1=get_database_state(database1),
        state2=DatabaseState(
            {
                DatabaseEntry(**json1["molecule"]): 1,
                DatabaseEntry(**to_hashable(json3["matrix"])): 1,
            }
        ),
    )
