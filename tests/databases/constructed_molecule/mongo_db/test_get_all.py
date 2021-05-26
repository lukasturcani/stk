import stk

from tests.utilities import is_equivalent_constructed_molecule


def test_get_all(mongo_client):
    """
    Test iteration over all molecules.

    """

    database_name = '_test_get_entries_constructed_molecule'
    mongo_client.drop_database(database_name)

    key_maker = stk.Inchi()
    jsonizer = stk.ConstructedMoleculeJsonizer(
        key_makers=(key_maker, )
    )

    database = stk.ConstructedMoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        jsonizer=jsonizer,
        put_lru_cache_size=0,
        get_lru_cache_size=0,
    )

    molecules = [
        stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(stk.BuildingBlock(
                    smiles='BrCCCBr',
                    functional_groups=[stk.BromoFactory()]),
                ),
                repeating_unit='A',
                num_repeating_units=3,
            ),
        ),
        stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrCCBr',
                        functional_groups=[stk.BromoFactory()]
                    ),
                    stk.BuildingBlock(
                        smiles='BrCNCBr',
                        functional_groups=[stk.BromoFactory()]
                    ),
                ),
                repeating_unit='AB',
                num_repeating_units=2,
            ),
        ),
    ]
    molecules_by_key = {
        key_maker.get_key(molecule): molecule
        for molecule in molecules
    }

    for molecule in molecules:
        database.put(molecule)

    for i, retrieved in enumerate(database.get_all()):
        key = key_maker.get_key(retrieved)
        molecule = molecules_by_key[key]
        is_equivalent_constructed_molecule(
            molecule.with_canonical_atom_ordering(),
            retrieved.with_canonical_atom_ordering(),
        )

    # Check number of molecules.
    assert i+1 == len(molecules)
