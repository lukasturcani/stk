"""
Molecule DeJSONizer
===================

"""

import numpy as np

from stk.molecular import Molecule

from .utilities import to_atom, to_bond


class MoleculeDejsonizer:
    """
    Abstract base class for creating molecules from JSONs.

    See Also
    --------
    :class:`.ConstructedMoleculeDejsonizer`

    Notes
    -----
    You might notice that the public methods of this abstract base
    class are implemented. These are just default implementations,
    which can be safely ignored or overridden, when implementing
    subclasses. However, the default implementation can be used
    directly, if it suits your needs.

    """

    def from_json(self, json):
        """
        Get a :class:`.Molecule` from a JSON.

        Parameters
        ----------
        json : :class:`dict`
            A JSON representation of a molecule.

        Returns
        -------
        :class:`.Molecule`
            The molecule held in `json`.

        """

        atoms = tuple(
            to_atom(atom_id, atom_json)
            for atom_id, atom_json in enumerate(json['molecule']['a'])
        )
        return Molecule(
            atoms=atoms,
            bonds=tuple(
                to_bond(atoms, bond_json)
                for bond_json in json['molecule']['b']
            ),
            position_matrix=np.array(json['matrix']['m']),
        )

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'{self.__class__.__name__}()'
