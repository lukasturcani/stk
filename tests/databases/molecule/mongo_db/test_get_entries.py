import stk
import pymongo

from tests.utilities import is_equivalent_molecule


def test_get_entries():
    """
    Test iteration over all entries.

    """

    database_name = '_test_get_entries_molecule'
    client = pymongo.MongoClient()
    client.drop_database(database_name)

    key_maker = stk.Inchi()
    jsonizer = stk.MoleculeJsonizer(key_makers=(key_maker, ))

    database = stk.MoleculeMongoDb(
        mongo_client=client,
        database=database_name,
        jsonizer=jsonizer,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
    )

    molecule_list = [
        stk.BuildingBlock('CCC').with_canonical_atom_ordering(),
        stk.BuildingBlock('BrCCCBr').with_canonical_atom_ordering(),
        stk.BuildingBlock('NCCN').with_canonical_atom_ordering(),
    ]
    molecule_key_dict = {
        key_maker.get_key(i): i for i in molecule_list
    }

    for molecule in molecule_list:
        database.put(molecule)

    for i, retrieved in enumerate(database.get_entries()):
        key = key_maker.get_key(retrieved)
        molecule = molecule_key_dict[key]
        is_equivalent_molecule(
            molecule1=molecule.with_canonical_atom_ordering(),
            molecule2=retrieved.with_canonical_atom_ordering(),
        )

    # Check number of molecules.
    assert i+1 == len(molecule_list)
