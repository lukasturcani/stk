import numpy as np
import stk

from ...utilities import (
    DatabaseEntry,
    DatabaseState,
    HashableDict,
    assert_database_state,
)
from .utilities import get_database_state


def to_hashable_matrix(json):
    json = dict(json)
    return {
        "m": tuple(tuple(row) for row in json.pop("m")),
        **json,
    }


def to_hashable_constructed_molecule(json):
    json = dict(json)
    return {
        "BB": tuple(map(HashableDict, json.pop("BB"))),
        "aI": tuple(tuple(v) for v in json.pop("aI")),
        "bI": tuple(tuple(v) for v in json.pop("bI")),
        "nBB": tuple(json.pop("nBB")),
        **json,
    }


def test_update_1(mongo_client):
    """
    Test that existing entries are updated.

    """

    database_name = "_test_update_1"
    mongo_client.drop_database(database_name)

    database = stk.ConstructedMoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
    )
    jsonizer = stk.ConstructedMoleculeJsonizer()

    molecule = stk.BuildingBlock(
        smiles="BrCCBr",
        functional_groups=[stk.BromoFactory()],
    ).with_canonical_atom_ordering()

    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            # Use it as a building block twice, to make sure it is
            # not repeatedly added to the molecules database.
            building_blocks=(molecule, molecule),
            repeating_unit="AB",
            num_repeating_units=2,
        ),
    ).with_canonical_atom_ordering()
    json = jsonizer.to_json(polymer)

    database.put(polymer)
    assert_database_state(
        state1=get_database_state(database),
        state2=DatabaseState(
            {
                DatabaseEntry(**json["molecule"]): 1,
                DatabaseEntry(**to_hashable_matrix(json["matrix"])): 1,
                DatabaseEntry(**json["buildingBlocks"][0]["molecule"]): 1,
                DatabaseEntry(
                    **to_hashable_matrix(
                        json=json["buildingBlocks"][0]["matrix"],
                    )
                ): 1,
                DatabaseEntry(
                    **to_hashable_constructed_molecule(
                        json=json["constructedMolecule"],
                    )
                ): 1,
            }
        ),
    )

    polymer2 = polymer.with_position_matrix(
        position_matrix=np.zeros((polymer.get_num_atoms(), 3)),
    )
    json2 = jsonizer.to_json(polymer2)

    database.put(polymer2)
    assert_database_state(
        state1=get_database_state(database),
        state2=DatabaseState(
            {
                DatabaseEntry(**json["molecule"]): 1,
                DatabaseEntry(**to_hashable_matrix(json2["matrix"])): 1,
                DatabaseEntry(**json["buildingBlocks"][0]["molecule"]): 1,
                DatabaseEntry(
                    **to_hashable_matrix(
                        json=json["buildingBlocks"][0]["matrix"],
                    )
                ): 1,
                DatabaseEntry(
                    **to_hashable_constructed_molecule(
                        json=json["constructedMolecule"],
                    )
                ): 1,
            }
        ),
    )


def test_update_2(mongo_client):
    """
    Test that existing entries are updated.

    In this test, your first create two separate entries, using
    different molecule keys. You then update both at the same time,
    with a database which uses both molecule keys.

    """

    database_name = "_test_update_2"
    mongo_client.drop_database(database_name)

    jsonizer1 = stk.ConstructedMoleculeJsonizer(
        key_makers=(stk.InchiKey(),),
    )
    database1 = stk.ConstructedMoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        jsonizer=jsonizer1,
    )

    jsonizer2 = stk.ConstructedMoleculeJsonizer(
        key_makers=(stk.Smiles(),),
    )
    database2 = stk.ConstructedMoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        jsonizer=jsonizer2,
    )

    jsonizer3 = stk.ConstructedMoleculeJsonizer(
        key_makers=(
            stk.InchiKey(),
            stk.Smiles(),
        ),
    )
    database3 = stk.ConstructedMoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        jsonizer=jsonizer3,
    )

    molecule = stk.BuildingBlock(
        smiles="BrCCCBr",
        functional_groups=[stk.BromoFactory()],
    ).with_canonical_atom_ordering()

    polymer1 = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            # Use it as a building block twice, to make sure it is
            # not repeatedly added to the molecules database.
            building_blocks=(molecule, molecule),
            repeating_unit="AB",
            num_repeating_units=2,
        ),
    ).with_canonical_atom_ordering()
    json1 = jsonizer1.to_json(polymer1)

    database1.put(polymer1)
    assert_database_state(
        state1=get_database_state(database1),
        state2=DatabaseState(
            {
                DatabaseEntry(**json1["molecule"]): 1,
                DatabaseEntry(**to_hashable_matrix(json1["matrix"])): 1,
                DatabaseEntry(**json1["buildingBlocks"][0]["molecule"]): 1,
                DatabaseEntry(
                    **to_hashable_matrix(
                        json=json1["buildingBlocks"][0]["matrix"],
                    )
                ): 1,
                DatabaseEntry(
                    **to_hashable_constructed_molecule(
                        json=json1["constructedMolecule"],
                    )
                ): 1,
            }
        ),
    )

    # Should add another entry, as a different key maker is used.
    polymer2 = polymer1.with_position_matrix(
        position_matrix=np.zeros((polymer1.get_num_atoms(), 3)),
    )
    json2 = jsonizer2.to_json(polymer2)

    database2.put(polymer2)
    assert_database_state(
        state1=get_database_state(database1),
        state2=DatabaseState(
            {
                DatabaseEntry(**json1["molecule"]): 1,
                DatabaseEntry(**to_hashable_matrix(json1["matrix"])): 1,
                DatabaseEntry(**json1["buildingBlocks"][0]["molecule"]): 1,
                DatabaseEntry(
                    **to_hashable_matrix(
                        json=json1["buildingBlocks"][0]["matrix"],
                    )
                ): 1,
                DatabaseEntry(
                    **to_hashable_constructed_molecule(
                        json=json1["constructedMolecule"],
                    )
                ): 1,
                DatabaseEntry(**json2["molecule"]): 1,
                DatabaseEntry(**to_hashable_matrix(json2["matrix"])): 1,
                DatabaseEntry(**json2["buildingBlocks"][0]["molecule"]): 1,
                DatabaseEntry(
                    **to_hashable_matrix(
                        json=json2["buildingBlocks"][0]["matrix"],
                    )
                ): 1,
                DatabaseEntry(
                    **to_hashable_constructed_molecule(
                        json=json2["constructedMolecule"],
                    )
                ): 1,
            }
        ),
    )

    # Should update both entries.
    polymer3 = polymer1.with_position_matrix(
        position_matrix=np.zeros((polymer1.get_num_atoms(), 3)),
    )
    json3 = jsonizer3.to_json(polymer3)

    database3.put(polymer3)
    assert_database_state(
        state1=get_database_state(database1),
        state2=DatabaseState(
            {
                DatabaseEntry(**json3["molecule"]): 2,
                DatabaseEntry(**to_hashable_matrix(json3["matrix"])): 2,
                DatabaseEntry(**json3["buildingBlocks"][0]["molecule"]): 2,
                DatabaseEntry(
                    **to_hashable_matrix(
                        json=json3["buildingBlocks"][0]["matrix"],
                    )
                ): 2,
                DatabaseEntry(
                    **to_hashable_constructed_molecule(
                        json=json3["constructedMolecule"],
                    )
                ): 2,
                DatabaseEntry(**json3["molecule"]): 2,
                DatabaseEntry(**to_hashable_matrix(json3["matrix"])): 2,
                DatabaseEntry(**json3["buildingBlocks"][0]["molecule"]): 2,
                DatabaseEntry(
                    **to_hashable_matrix(
                        json=json3["buildingBlocks"][0]["matrix"],
                    )
                ): 2,
                DatabaseEntry(
                    **to_hashable_constructed_molecule(
                        json=json3["constructedMolecule"],
                    )
                ): 2,
            }
        ),
    )


def test_update_3(mongo_client):
    """
    Test that existing entries are updated.

    In this test, your first create one entry with two keys. Then
    update the entry with databases, each using 1 different key.
    No duplicate entries should be made in the database this way.

    """

    database_name = "_test_update_3"
    mongo_client.drop_database(database_name)

    jsonizer1 = stk.ConstructedMoleculeJsonizer(
        key_makers=(
            stk.InchiKey(),
            stk.Smiles(),
        ),
    )
    database1 = stk.ConstructedMoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        jsonizer=jsonizer1,
    )

    jsonizer2 = stk.ConstructedMoleculeJsonizer(
        key_makers=(stk.InchiKey(),),
    )
    database2 = stk.ConstructedMoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        jsonizer=jsonizer2,
    )

    jsonizer3 = stk.ConstructedMoleculeJsonizer(
        key_makers=(stk.Smiles(),),
    )
    database3 = stk.ConstructedMoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        jsonizer=jsonizer3,
    )

    molecule = stk.BuildingBlock(
        smiles="BrCCCBr",
        functional_groups=[stk.BromoFactory()],
    ).with_canonical_atom_ordering()

    polymer1 = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            # Use it as a building block twice, to make sure it is
            # not repeatedly added to the molecules database.
            building_blocks=(molecule, molecule),
            repeating_unit="AB",
            num_repeating_units=2,
        ),
    ).with_canonical_atom_ordering()
    json1 = jsonizer1.to_json(polymer1)

    database1.put(polymer1)
    assert_database_state(
        state1=get_database_state(database1),
        state2=DatabaseState(
            {
                DatabaseEntry(**json1["molecule"]): 1,
                DatabaseEntry(**to_hashable_matrix(json1["matrix"])): 1,
                DatabaseEntry(**json1["buildingBlocks"][0]["molecule"]): 1,
                DatabaseEntry(
                    **to_hashable_matrix(
                        json=json1["buildingBlocks"][0]["matrix"],
                    )
                ): 1,
                DatabaseEntry(
                    **to_hashable_constructed_molecule(
                        json=json1["constructedMolecule"],
                    )
                ): 1,
            }
        ),
    )

    # Should update the entry.
    polymer2 = polymer1.with_position_matrix(
        position_matrix=np.zeros((polymer1.get_num_atoms(), 3)),
    )
    json2 = jsonizer2.to_json(polymer2)
    json2["matrix"] = dict(json1["matrix"])
    json2["matrix"]["m"] = jsonizer2.to_json(polymer2)["matrix"]["m"]

    database2.put(polymer2)
    assert_database_state(
        state1=get_database_state(database1),
        state2=DatabaseState(
            {
                DatabaseEntry(**json1["molecule"]): 1,
                DatabaseEntry(**to_hashable_matrix(json2["matrix"])): 1,
                DatabaseEntry(**json1["buildingBlocks"][0]["molecule"]): 1,
                DatabaseEntry(
                    **to_hashable_matrix(
                        json=json1["buildingBlocks"][0]["matrix"],
                    )
                ): 1,
                DatabaseEntry(
                    **to_hashable_constructed_molecule(
                        json=json1["constructedMolecule"],
                    )
                ): 1,
            }
        ),
    )

    # Should also update the entry.
    polymer3 = polymer1.with_position_matrix(
        position_matrix=np.zeros((polymer1.get_num_atoms(), 3)),
    )
    json3 = jsonizer3.to_json(polymer3)
    json3["matrix"] = dict(json1["matrix"])
    json3["matrix"]["m"] = jsonizer3.to_json(polymer3)["matrix"]["m"]

    database3.put(polymer3)
    assert_database_state(
        state1=get_database_state(database1),
        state2=DatabaseState(
            {
                DatabaseEntry(**json1["molecule"]): 1,
                DatabaseEntry(**to_hashable_matrix(json3["matrix"])): 1,
                DatabaseEntry(**json1["buildingBlocks"][0]["molecule"]): 1,
                DatabaseEntry(
                    **to_hashable_matrix(
                        json=json1["buildingBlocks"][0]["matrix"],
                    )
                ): 1,
                DatabaseEntry(
                    **to_hashable_constructed_molecule(
                        json=json1["constructedMolecule"],
                    )
                ): 1,
            }
        ),
    )
