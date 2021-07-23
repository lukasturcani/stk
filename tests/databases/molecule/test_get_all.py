import stk

from tests.utilities import is_equivalent_molecule


def test_get_all(mongo_client):
    """
    Test iteration over all entries.

    """

    database_name = '_test_get_entries_molecule'
    mongo_client.drop_database(database_name)

    key_maker = stk.Inchi()
    jsonizer = stk.MoleculeJsonizer(key_makers=(key_maker, ))
    database1 = stk.MoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        jsonizer=jsonizer,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        indices=(key_maker.get_key_name(), ),
    )

    key_maker = stk.Smiles()
    jsonizer = stk.MoleculeJsonizer(key_makers=(key_maker, ))
    database2 = stk.MoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        jsonizer=jsonizer,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        indices=(key_maker.get_key_name(), ),
    )

    all_molecules = (
        stk.BuildingBlock('CCC').with_canonical_atom_ordering(),
        stk.BuildingBlock('BrCCCBr').with_canonical_atom_ordering(),
        stk.BuildingBlock('NCCN').with_canonical_atom_ordering(),
        stk.BuildingBlock('CCCCCC').with_canonical_atom_ordering(),
        stk.BuildingBlock('NCCCCN').with_canonical_atom_ordering(),
    )

    # Two sets of molecules with two distinct molecules each, and one
    # shared molecule to test for overlap.
    molecules1 = (all_molecules[0], all_molecules[1], all_molecules[4])
    molecules2 = (all_molecules[2], all_molecules[3], all_molecules[4])

    for molecule in molecules1:
        database1.put(molecule)

    for molecule in molecules2:
        database2.put(molecule)

    # A new database with both indices.
    key_makers = (stk.Inchi(), stk.Smiles())
    jsonizer = stk.MoleculeJsonizer(key_makers=key_makers)
    database3 = stk.MoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        jsonizer=jsonizer,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        indices=(i.get_key_name() for i in key_makers),
    )

    molecules1_by_key = {
        key_makers[0].get_key(molecule): molecule
        for molecule in molecules1
    }
    molecules2_by_key = {
        key_makers[1].get_key(molecule): molecule
        for molecule in molecules2
    }

    for i, retrieved in enumerate(database3.get_all()):
        try:
            key = key_makers[0].get_key(retrieved)
            molecule = molecules1_by_key[key]
        except KeyError:
            key = key_makers[1].get_key(retrieved)
            molecule = molecules2_by_key[key]
        is_equivalent_molecule(
            molecule1=molecule.with_canonical_atom_ordering(),
            molecule2=retrieved.with_canonical_atom_ordering(),
        )

    # Check number of molecules.
    assert i+1 == len(all_molecules)
