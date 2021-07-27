import stk

from tests.utilities import is_equivalent_constructed_molecule


def test_get_all(mongo_client):
    """
    Test iteration over all molecules.

    """

    database_name = '_test_get_entries_constructed_molecule'
    mongo_client.drop_database(database_name)

    inchi = stk.Inchi()
    inchi_database = stk.ConstructedMoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        jsonizer=stk.ConstructedMoleculeJsonizer(
            key_makers=(inchi, ),
        ),
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        indices=(inchi.get_key_name(), ),
    )

    smiles = stk.Smiles()
    smiles_database = stk.ConstructedMoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        jsonizer=stk.ConstructedMoleculeJsonizer(
            key_makers=(smiles, ),
        ),
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        indices=(smiles.get_key_name(), ),
    )

    inchi_and_smiles_database = stk.ConstructedMoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        jsonizer=stk.ConstructedMoleculeJsonizer(
            key_makers=(smiles, inchi),
        ),
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        indices=(),
    )

    building_blocks = (
        stk.BuildingBlock('BrC#CBr', [stk.BromoFactory()]),
        stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
        stk.BuildingBlock('BrCCCBr', [stk.BromoFactory()]),
        stk.BuildingBlock('BrCCCCBr', [stk.BromoFactory()]),
        stk.BuildingBlock('BrCNCBr', [stk.BromoFactory()]),
        stk.BuildingBlock('BrCCNCBr', [stk.BromoFactory()]),
    )

    all_molecules = [
        stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(building_blocks[0], ),
                repeating_unit='A',
                num_repeating_units=3,
            ),
        ),
        stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    building_blocks[0], building_blocks[1]
                ),
                repeating_unit='AB',
                num_repeating_units=3,
            ),
        ),
        stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(building_blocks[2], ),
                repeating_unit='A',
                num_repeating_units=3,
            ),
        ),
        stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    building_blocks[2], building_blocks[3]
                ),
                repeating_unit='AB',
                num_repeating_units=3,
            ),
        ),
        stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(building_blocks[4], ),
                repeating_unit='A',
                num_repeating_units=3,
            ),
        ),
        stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    building_blocks[4], building_blocks[5]
                ),
                repeating_unit='AB',
                num_repeating_units=3,
            ),
        ),
    ]

    inchi_molecules = all_molecules[:2]
    smiles_molecules = all_molecules[2:-2]
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
    inchi_key_database = stk.ConstructedMoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        jsonizer=stk.ConstructedMoleculeJsonizer(
            key_makers=(stk.InchiKey(), ),
        ),
        indices=(),
    )
    for i, retrieved in enumerate(inchi_key_database.get_all()):
        expected = smiles_to_molecule[smiles.get_key(retrieved)]
        is_equivalent_constructed_molecule(
            constructed_molecule1=(
                expected.with_canonical_atom_ordering()
            ),
            constructed_molecule2=(
                retrieved.with_canonical_atom_ordering()
            ),
        )

    assert i+1 == len(all_molecules)
