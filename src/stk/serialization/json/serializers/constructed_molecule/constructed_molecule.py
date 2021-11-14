"""
Constructed Molecule JSONizer
=============================

"""

from stk.molecular import InchiKey, MoleculeKeyMaker

from ..molecule import MoleculeJsonizer


class ConstructedMoleculeJsonizer:
    """
    Abstract base class for creating JSONs of constructed molecules.

    See Also
    --------
    :class:`.MoleculeJsonizer`

    Notes
    -----
    You might notice that the public methods of this abstract base
    class are implemented. These are just default implementations,
    which can be safely ignored or overridden, when implementing
    subclasses. However, the default implementation can be used
    directly, if it suits your needs.

    Examples
    --------
    *Converting a Constructed Molecule to JSON*

    You want get a JSON representation of a
    :class:`.ConstructedMolecule`

    .. testcode:: converting-a-molecule-to-json

        import stk

        # Make the molecule you want jsonize.
        polymer = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
                ),
                repeating_unit='A',
                num_repeating_units=3,
            )
        )

        # Make a JSONizer.
        jsonizer = stk.ConstructedMoleculeJsonizer()
        # Get the JSON.
        json = jsonizer.to_json(polymer)

    *Adding Additional Molecular Keys*

    Apart from atoms, bonds and the position matrix, the JSON
    representation holds additional fields, one for each
    :class:`.MoleculeKeyMaker` provided to the initializer

    .. testcode:: adding-additional-molecular-keys

        import stk

        # Make the molecule you want jsonize.
        polymer = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
                ),
                repeating_unit='A',
                num_repeating_units=3,
            )
        )

        # Make a JSONizer.
        jsonizer = stk.ConstructedMoleculeJsonizer()
        # Get the JSON.
        json = jsonizer.to_json(polymer)

    In this case, ``json`` will look something like

    .. code-block:: python

        {
            # A tuple of JSON atom representations.
            'atoms': (...),

            # A tuple of JSON bond representations.
            'bonds': (...),

            'InChI': 'The InChI of the molecule',
            'InChIKey': 'The InChIKey of the molecule',
        }

    For every :class:`.MoleculeKeyMaker` provided to `key_makers`,
    a new key will be added to the JSON representation, with its name
    given by :meth:`.MoleculeKeyMaker.get_key_name` and the value
    given by :meth:`.MoleculeKeyMaker.get_key`.

    """

    def __init__(
        self,
        key_makers=(InchiKey(), ),
    ):
        """
        Initializes a :class:`.ConstructedMoleculeJsonizer`.

        Parameters
        ----------
        key_makers : :class:`tuple` of :class:`.MoleculeKeyMaker`
            Used to make the keys of molecules, which should be
            included in their JSON representations. Keys allow
            molecular data to reference itself when split across
            multiple JSONs.

        """

        self._jsonizer = MoleculeJsonizer(key_makers=())
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
                key_maker.get_key_name():
                    key_maker.get_key(building_block)
                for key_maker in self._key_makers
                if isinstance(key_maker, MoleculeKeyMaker)
            }

        building_block_indices = {
            building_block: index
            for index, building_block
            in enumerate(molecule.get_building_blocks())
        }
        building_block_indices[None] = None

        def atom_info_to_json(atom_info):
            if atom_info.get_building_block() is None:
                return (
                    None,
                    None,
                    None,
                )

            return (
                building_block_indices[atom_info.get_building_block()],
                atom_info.get_building_block_id(),
                atom_info.get_building_block_atom().get_id(),
            )

        def bond_info_to_json(bond_info):
            return (
                building_block_indices[bond_info.get_building_block()],
                bond_info.get_building_block_id(),
            )

        molecule_json = self._jsonizer.to_json(molecule)
        constructed_molecule_json = {
            'BB': tuple(map(
                get_keys,
                molecule.get_building_blocks(),
            )),
            'aI': tuple(map(
                atom_info_to_json,
                molecule.get_atom_infos(),
            )),
            'bI': tuple(map(
                bond_info_to_json,
                molecule.get_bond_infos(),
            )),
            'nBB': tuple(map(
                molecule.get_num_building_block,
                molecule.get_building_blocks(),
            )),
        }
        for key_maker in self._key_makers:
            key_name = key_maker.get_key_name()
            key = key_maker.get_key(molecule)
            molecule_json['molecule'][key_name] = key
            molecule_json['matrix'][key_name] = key
            constructed_molecule_json[key_name] = key

        building_block_jsons = tuple(map(
            self._jsonizer.to_json,
            molecule.get_building_blocks(),
        ))

        def is_molecule_key_maker(key_maker):
            return isinstance(key_maker, MoleculeKeyMaker)

        for key_maker in filter(
            is_molecule_key_maker,
            self._key_makers,
        ):
            key_name = key_maker.get_key_name()
            for building_block, json in zip(
                molecule.get_building_blocks(),
                building_block_jsons,
            ):
                key = key_maker.get_key(building_block)
                json['molecule'][key_name] = key
                json['matrix'][key_name] = key

        return {
            'molecule': molecule_json['molecule'],
            'constructedMolecule': constructed_molecule_json,
            'matrix': molecule_json['matrix'],
            'buildingBlocks': building_block_jsons,
        }

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'{self.__class__.__name__}({self._key_makers!r})'
