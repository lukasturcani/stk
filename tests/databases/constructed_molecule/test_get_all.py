import stk

from tests.utilities import is_equivalent_constructed_molecule


def test_get_all(mongo_client):
    """
    Test iteration over all molecules.

    """

    database_name = '_test_get_entries_constructed_molecule'
    mongo_client.drop_database(database_name)

    key_maker = stk.Inchi()
    jsonizer = stk.ConstructedMoleculeJsonizer((key_maker, ))
    database1 = stk.ConstructedMoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        jsonizer=jsonizer,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        indices=(key_maker.get_key_name(), ),
    )

    key_maker = stk.Smiles()
    jsonizer = stk.ConstructedMoleculeJsonizer((key_maker, ))
    database2 = stk.ConstructedMoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        jsonizer=jsonizer,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        indices=(key_maker.get_key_name(), ),
    )

    building_blocks = (
        stk.BuildingBlock('BrCBr', [stk.BromoFactory()]),
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
                num_repeating_units=2,
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
                num_repeating_units=2,
            ),
        ),
        stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(building_blocks[4], ),
                repeating_unit='A',
                num_repeating_units=3,
            ),
        ),
    ]

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
    jsonizer = stk.ConstructedMoleculeJsonizer(key_makers)
    database3 = stk.ConstructedMoleculeMongoDb(
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
        is_equivalent_constructed_molecule(
            constructed_molecule1=(
                molecule.with_canonical_atom_ordering()
            ),
            constructed_molecule2=(
                retrieved.with_canonical_atom_ordering()
            ),
        )

    # Check number of molecules.
    assert i+1 == len(all_molecules)
