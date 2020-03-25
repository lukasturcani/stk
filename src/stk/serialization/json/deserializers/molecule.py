"""
Molecule DeJSONizer
===================

"""

from stk.molecular import Molecule

from .utilities import to_atom, to_bond


class MoleculeDejsonizer:
    """
    Creates :class:`.Molecule` instances from JSON representations.

    """

    # Keep this empty __init__() method for the docstring.
    def __init__(self):
        """
        Initialize a :class`.MoleculeDejsonizer`.

        """

        return

    def from_json(self, json, position_matrix):
        """
        Get a :class:`.Molecule` from a JSON.

        Parameters
        ----------
        json : :class:`dict`
            A JSON representation of a molecule.

        position_matrix : :class:`numpy.ndarray`
            The position matrix of the created molecule.

        Returns
        -------
        :class:`.Molecule`
            The molecule held in `json`.

        """

        atoms = tuple(
            to_atom(atom_id, atom_json)
            for atom_id, atom_json in enumerate(json['a'])
        )
        return Molecule(
            atoms=atoms,
            bonds=tuple(
                to_bond(atoms, bond_json) for bond_json in json['b']
            ),
            position_matrix=position_matrix,
        )

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'{self.__class__.__name__}()'
