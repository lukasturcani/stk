import stk
import pymongo

from tests.utilities import is_equivalent_molecule


def test_get_all():
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

    molecules = (
        stk.BuildingBlock('CCC').with_canonical_atom_ordering(),
        stk.BuildingBlock('BrCCCBr').with_canonical_atom_ordering(),
        stk.BuildingBlock('NCCN').with_canonical_atom_ordering(),
    )
    molecules_by_key = {
        key_maker.get_key(molecule): molecule
        for molecule in molecules
    }

    for molecule in molecules:
        database.put(molecule)

    for i, retrieved in enumerate(database.get_all()):
        key = key_maker.get_key(retrieved)
        molecule = molecules_by_key[key]
        is_equivalent_molecule(
            molecule1=molecule.with_canonical_atom_ordering(),
            molecule2=retrieved.with_canonical_atom_ordering(),
        )

    # Check number of molecules.
    assert i+1 == len(molecules)
