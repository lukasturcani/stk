"""
Constructed Molecule JSONizer
=============================

"""

from stk.molecular import (
    MoleculeKeyMaker,
    InchiKey,
    ConstructedMoleculeKeyMaker,
)


class ConstructedMoleculeJsonizer:
    """
    Serializes :class:`.ConstructedMolecule` instances into JSON.

    See Also
    --------
    :class:`.MoleculeJsonizer`

    Examples
    --------
    It is first necessary to read the docstring of
    :class:`.MoleculeJsonizer`, since this class is almost always used
    together with it, and not as a replacement for it.

    Now, if you want get a JSON representation of a
    :class:`.ConstructedMolecule` you can do something like

    .. code-block:: python

        import stk

        # Make the molecule you want jsonize.

        polymer_topology_graph = stk.polymer.Linear(
            building_blocks=(
                stk.BuildingBlock('BrCCBr', [stk.BromoFactory()],),
            ),
            repeating_unit='A',
            num_repeating_units=3,
        )
        polymer = stk.ConstructedMolecule(polymer_topology_graph)


        # Make a JSONizer.
        jsonizer = stk.ConstructedMoleculeJsonizer()
        # Get the JSON.
        json = jsonizer.to_json(polymer)

    At this point, you might think that ``json`` holds the entire
    JSON representation of ``polymer``. This is not the case. It
    only holds the information that is exclusive to
    :class:`.ConstructedMolecule`, and only holds references to the
    information which is relevant to a :class:`.Molecule`.

    If you want get a JSON representation of the information relevant
    to a :class:`.Molecule`, like the atoms and bonds, you have
    to use a :class:`.MoleculeJsonizer` in addition to the
    :class:`.ConstructedMoleculeJsonizer`

    .. code-block:: python

        molecule_jsonizer = stk.MoleculeJsonizer()
        # Holds JSON representation of information relevant to
        # Molecule objects.
        molecule_json = molecule_jsonizer.to_json(polymer)

    The reason for this, is that it allows the building of
    JSON databases, like MongoDB, in a sensible way. This means
    information relevant to the :class:`.Molecule` class is held in a
    separate collection to the extra data carried by a
    :class:`.ConstructedMolecule`, which is held in its own
    collection.

    Here the purpose of the `key_makers` parameter becomes apparent.
    The keys made by the `key_makers` can be used to reference the
    relevant molecular data across collections. For example, instead of
    holding the building blocks of a :class:`.ConstructedMolecule`
    directly, the JSON representation merely uses `key_makers` to
    create the keys, which can be used to find them in the separate
    collection.

    """

    def __init__(
        self,
        key_makers=(
            InchiKey(),
            ConstructedMoleculeKeyMaker(),
        ),
    ):
        """
        Initializes a :class:`.ConstructedMoleculeJsonizer`.

        Parameters
        ----------
        key_makers : :class:`tuple` of \
                :class:`.MoleculeKeyMaker` and \
                :class:`.ConstructedMoleculeKeyMaker`
            Used to make the keys of molecules, which should be
            included in their JSON representations.

        """

        self._key_makers = key_makers

    def to_json(self, molecule):
        """
        Serialize `molecule`.

        Parameters
        ----------
        molecule : :class:`.ConstructedMolecule`
            The constructed molecule to serialize.

        Returns
        -------
        :class:`dict`
            A JSON representation of `molecule`.

        """

        def get_keys(building_block):
            return {
                key_maker.get_key_name(): key_maker.get_key(molecule)
                for key_maker in self._key_makers
                if isinstance(key_maker, MoleculeKeyMaker)
            }

        building_block_indices = {
            building_block: index
            for index, building_block
            in enumerate(molecule.get_building_blocks())
        }

        def atom_info_to_json(atom_info):
            return {
                'atom': atom_info.get_atom().get_id(),
                'building_block': building_block_indices[
                    atom_info.get_building_block()
                ],
                'building_block_id': atom_info.get_building_block_id(),
            }

        def bond_info_to_json(bond_data):
            bond_id, bond_info = bond_data
            return {
                'bond': bond_id,
                'building_block': building_block_indices[
                    bond_info.get_building_block()
                ],
                'building_block_id': bond_info.get_building_block_id(),
            }

        json = {
            'building_blocks': tuple(map(
                get_keys,
                molecule.get_building_blocks(),
            )),
            'atom_infos': tuple(map(
                atom_info_to_json,
                molecule.get_atom_infos(),
            )),
            'bond_infos': tuple(map(
                bond_info_to_json,
                enumerate(molecule.get_bond_infos()),
            )),
            'num_building_blocks': {
                index: molecule.get_num_building_block(building_block)
                for index, building_block
                in enumerate(molecule.get_building_blocks())
            },
        }
        for key_maker in self._key_makers:
            json[key_maker.get_key_name()] = (
                key_maker.get_key(molecule)
            )
        return json
