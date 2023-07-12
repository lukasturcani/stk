from collections.abc import Iterable

import numpy as np

from stk._internal.atom_info import AtomInfo
from stk._internal.bond_info import BondInfo
from stk._internal.constructed_molecule import ConstructedMolecule
from stk._internal.key_makers.inchi_key import InchiKey
from stk._internal.key_makers.molecule import MoleculeKeyMaker
from stk._internal.molecule import Molecule

from .molecule import MoleculeDejsonizer, MoleculeJsonizer
from .utilities import to_atom, to_atom_info, to_bond, to_bond_info


class ConstructedMoleculeJsonizer:
    """
    Abstract base class for creating JSONs of constructed molecules.

    See Also:
        :class:`.MoleculeJsonizer`

    Notes:
        You might notice that the public methods of this abstract base
        class are implemented. These are just default implementations,
        which can be safely ignored or overridden, when implementing
        subclasses. However, the default implementation can be used
        directly, if it suits your needs.

    Examples:
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
        key_makers: Iterable[MoleculeKeyMaker] = (InchiKey(),),
    ) -> None:
        """
        Parameters:
            key_makers (list[MoleculeKeyMaker]):
                Used to make the keys of molecules, which should be
                included in their JSON representations. Keys allow
                molecular data to reference itself when split across
                multiple JSONs.
        """
        self._jsonizer = MoleculeJsonizer(key_makers=())
        self._key_makers = tuple(key_makers)

    def to_json(self, molecule: ConstructedMolecule) -> dict:
        """
        Serialize `molecule`.

        Parameters:
            molecule:
                The constructed molecule to serialize.

        Returns:
            A JSON representation of `molecule`.
        """

        def get_keys(building_block: Molecule) -> dict:
            return {
                key_maker.get_key_name(): key_maker.get_key(building_block)
                for key_maker in self._key_makers
                if isinstance(key_maker, MoleculeKeyMaker)
            }

        building_block_indices: dict[Molecule | None, int | None] = {
            building_block: index
            for index, building_block in enumerate(
                molecule.get_building_blocks()
            )
        }
        building_block_indices[None] = None

        def atom_info_to_json(
            atom_info: AtomInfo,
        ) -> tuple[int | None, int | None, int | None]:
            building_block = atom_info.get_building_block()
            building_block_atom = atom_info.get_building_block_atom()
            if building_block is not None and building_block_atom is not None:
                return (
                    building_block_indices[building_block],
                    atom_info.get_building_block_id(),
                    building_block_atom.get_id(),
                )
            return (
                None,
                None,
                None,
            )

        def bond_info_to_json(
            bond_info: BondInfo,
        ) -> tuple[int | None, int | None]:
            return (
                building_block_indices[bond_info.get_building_block()],
                bond_info.get_building_block_id(),
            )

        molecule_json = self._jsonizer.to_json(molecule)
        constructed_molecule_json = {
            "BB": tuple(
                map(
                    get_keys,
                    molecule.get_building_blocks(),
                )
            ),
            "aI": tuple(
                map(
                    atom_info_to_json,
                    molecule.get_atom_infos(),
                )
            ),
            "bI": tuple(
                map(
                    bond_info_to_json,
                    molecule.get_bond_infos(),
                )
            ),
            "nBB": tuple(
                map(
                    molecule.get_num_building_block,
                    molecule.get_building_blocks(),
                )
            ),
        }
        for key_maker in self._key_makers:
            key_name = key_maker.get_key_name()
            key = key_maker.get_key(molecule)
            molecule_json["molecule"][key_name] = key  # type: ignore
            molecule_json["matrix"][key_name] = key  # type: ignore
            constructed_molecule_json[key_name] = key  # type: ignore

        building_block_jsons = tuple(
            map(
                self._jsonizer.to_json,
                molecule.get_building_blocks(),
            )
        )

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
                json["molecule"][key_name] = key  # type: ignore
                json["matrix"][key_name] = key  # type: ignore

        return {
            "molecule": molecule_json["molecule"],
            "constructedMolecule": constructed_molecule_json,
            "matrix": molecule_json["matrix"],
            "buildingBlocks": building_block_jsons,
        }

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self._key_makers!r})"


class ConstructedMoleculeDejsonizer:
    """
    Abstract base class for creating constructed molecules from JSONs.

    See Also:
        :class:`.MoleculeDejsonizer`

    Notes:
        You might notice that the public methods of this abstract base
        class are implemented. These are just default implementations,
        which can be safely ignored or overridden, when implementing
        subclasses. However, the default implementation can be used
        directly, if it suits your needs.
    """

    def __init__(self) -> None:
        self._dejsonizer = MoleculeDejsonizer()

    def from_json(self, json: dict) -> ConstructedMolecule:
        """
        Get a :class:`.ConstructedMolecule` from a JSON.

        Parameters:
            json:
                A JSON of the constructed molecule.

        Returns:
            The constructed molecule.
        """
        building_blocks = tuple(
            map(
                self._dejsonizer.from_json,
                json["buildingBlocks"],
            )
        )

        num_building_blocks = {
            building_block: num
            for building_block, num in zip(
                building_blocks,
                json["constructedMolecule"]["nBB"],
            )
        }
        atoms = tuple(
            to_atom(atom_id, atom_json)
            for atom_id, atom_json in enumerate(json["molecule"]["a"])
        )
        bonds = tuple(
            to_bond(atoms, bond_json) for bond_json in json["molecule"]["b"]
        )
        return ConstructedMolecule.init(
            atoms=atoms,
            bonds=bonds,
            position_matrix=np.array(json["matrix"]["m"]),
            atom_infos=tuple(
                to_atom_info(
                    building_blocks=building_blocks,
                    atom=atoms[atom_id],
                    json=atom_info_json,
                )
                for atom_id, atom_info_json in enumerate(
                    json["constructedMolecule"]["aI"]
                )
            ),
            bond_infos=tuple(
                to_bond_info(
                    building_blocks=building_blocks,
                    bond=bonds[bond_id],
                    json=bond_info_json,
                )
                for bond_id, bond_info_json in enumerate(
                    json["constructedMolecule"]["bI"]
                )
            ),
            num_building_blocks=num_building_blocks,
        )

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"
