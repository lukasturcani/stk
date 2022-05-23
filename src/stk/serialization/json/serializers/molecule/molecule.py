"""
Molecule JSONizer
=================

"""

from __future__ import annotations

import typing

from stk.molecular import InchiKey, Molecule, MoleculeKeyMaker
from stk.utilities import OneOrMany

from .utilities import AtomJson, BondJson, atom_to_json, bond_to_json


class _MolecularGraphJson(typing.TypedDict):
    a: tuple[AtomJson, ...]
    b: tuple[BondJson, ...]


class _PositionMatrixJson(typing.TypedDict):
    m: list[list[float]]


class _MoleculeJson(typing.TypedDict):
    molecule: _MolecularGraphJson
    matrix: _PositionMatrixJson


class MoleculeJsonizer:
    """
    Abstract base class for creating JSONs of molecules.

    See Also:

        :class:`.ConstructedMoleculeJsonizer`

    Notes:

        You might notice that the public methods of this abstract base
        class are implemented. These are just default implementations,
        which can be safely ignored or overridden, when implementing
        subclasses. However, the default implementation can be used
        directly, if it suits your needs.

    Examples:

        *Converting a Molecule to JSON*

        You want to create a JSON representation of a molecule

        .. testcode:: converting-a-molecule-to-json

            import stk

            jsonizer = stk.MoleculeJsonizer()
            json = jsonizer.to_json(stk.BuildingBlock('NCCN'))

        *Adding Additional Molecular Keys*

        Apart from atoms, bonds and the position matrix, the JSON
        representation holds additional fields, one for each
        :class:`.MoleculeKeyMaker` provided to the initializer

        .. testcode:: adding-additional-molecular-keys

            import stk

            jsonizer = stk.MoleculeJsonizer(
                key_makers=(
                    stk.Inchi(),
                    stk.InchiKey(),
                ),
            )
            json = jsonizer.to_json(stk.BuildingBlock('NCCN'))

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
        a new key will be added to the JSON representation, with its
        name given by :meth:`.MoleculeKeyMaker.get_key_name` and the
        value given by :meth:`.MoleculeKeyMaker.get_key`.

    """

    def __init__(
        self,
        key_makers: OneOrMany[MoleculeKeyMaker] = (InchiKey(), ),
    ) -> None:
        """
        Initialize a :class:`.MoleculeJsonizer` instance.

        Parameters:

            key_makers:
                Used to make the keys of molecules, which should be
                included in their JSON representations. Keys allow
                molecular data to reference itself when split across
                multiple JSONs.

        """

        if isinstance(key_makers, MoleculeKeyMaker):
            key_makers = (key_makers, )

        self._key_makers = tuple(key_makers)

    def to_json(
        self,
        molecule: Molecule,
    ) -> _MoleculeJson:
        """
        Return a JSON representation of `molecule`.

        Parameters:

            molecule:
                The molecule to serialize.

        Returns:

            A JSON representation of `molecule`.

        """

        json: _MolecularGraphJson = {
            'a': tuple(map(atom_to_json, molecule.get_atoms())),
            'b': tuple(map(bond_to_json, molecule.get_bonds())),
        }
        position_matrix: _PositionMatrixJson = {
            'm': molecule.get_position_matrix().tolist(),
        }
        for key_maker in self._key_makers:
            key_name = key_maker.get_key_name()
            key = key_maker.get_key(molecule)
            # TypedDict does not allow keys which are not listed in the
            # class definition. However, because we wish to create
            # indices on keys created by KeyMakers we need to add
            # these additional key-value pairs.
            json[key_name] = key  # type: ignore
            position_matrix[key_name] = key  # type: ignore
        return {
            'molecule': json,
            'matrix': position_matrix,
        }

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self._key_makers!r})'
