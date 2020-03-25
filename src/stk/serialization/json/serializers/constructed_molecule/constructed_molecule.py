"""
Constructed Molecule JSONizer
=============================

"""

from stk.molecular import (
    MoleculeKeyMaker,
    InchiKey,
    ConstructedMoleculeKeyMaker,
)
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
    You want get a JSON representation of a
    :class:`.ConstructedMolecule`

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

    Note that the JSON representation does not include the building
    blocks of the constructed molecule. Instead the building blocks
    are only referenced via the keys provided in the `key_makers`
    parameter.

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

        self._molecule_jsonizer = MoleculeJsonizer(key_makers=())
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
            return (
                building_block_indices[atom_info.get_building_block()],
                atom_info.get_building_block_id(),
            )

        def bond_info_to_json(bond_info):
            return (
                building_block_indices[bond_info.get_building_block()],
                bond_info.get_building_block_id(),
            )

        molecule_json = self._molecule_jsonizer.to_json(molecule)
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
            molecule_json['m'][key_name] = key
            constructed_molecule_json[key_name] = key
        return {
            'm': molecule_json['m'],
            'cM': constructed_molecule_json,
            'p': molecule_json['p'],
        }

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'{self.__class__.__name__}({self._key_makers!r})'
