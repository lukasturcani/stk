import stk

from tests.utilities import is_equivalent_molecule


def test_get_all(mongo_client):
    """
    Test iteration over all entries.

    """

    database_name = '_test_get_entries_molecule'
    mongo_client.drop_database(database_name)

    inchi = stk.Inchi()
    inchi_database = stk.MoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        jsonizer=stk.MoleculeJsonizer(
            key_makers=(inchi, ),
        ),
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        indices=(inchi.get_key_name(), ),
    )

    smiles = stk.Smiles()
    smiles_database = stk.MoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        jsonizer=stk.MoleculeJsonizer(
            key_makers=(smiles, ),
        ),
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        indices=(smiles.get_key_name(), ),
    )

    inchi_and_smiles_database = stk.MoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        jsonizer=stk.MoleculeJsonizer(
            key_makers=(smiles, inchi),
        ),
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        indices=(),
    )

    all_molecules = (
        stk.BuildingBlock('CCC'),
        stk.BuildingBlock('BrCCCBr'),
        stk.BuildingBlock('NCCN'),
        stk.BuildingBlock('CCCCC'),
        stk.BuildingBlock('NCCCCN'),
        stk.BuildingBlock('NCC(CCBr)CCN'),
        stk.BuildingBlock('NCCNCCC(CCBr)CCN'),
    )

    inchi_molecules = all_molecules[:3]
    smiles_molecules = all_molecules[3:-2]
    inchi_and_smiles_molecules = all_molecules[-2:]

    for molecule in inchi_molecules:
        inchi_database.put(molecule)

    for molecule in smiles_molecules:
        smiles_database.put(molecule)

    for molecule in inchi_and_smiles_molecules:
        inchi_and_smiles_database.put(molecule)

    smiles_to_molecule = {
        smiles.get_key(molecule): molecule
        for molecule in all_molecules
    }

    # Use an InChIKey database for get_all because none of the
    # molecules use this key, but this should not matter.
    inchi_key_database = stk.MoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        jsonizer=stk.MoleculeJsonizer(
            key_makers=(stk.InchiKey(), ),
        ),
        indices=(),
    )
    for i, retrieved in enumerate(inchi_key_database.get_all()):
        molecule = smiles_to_molecule[smiles.get_key(retrieved)]

        is_equivalent_molecule(
            molecule1=molecule.with_canonical_atom_ordering(),
            molecule2=retrieved.with_canonical_atom_ordering(),
        )

    # Check number of molecules.
    assert i+1 == len(all_molecules)
