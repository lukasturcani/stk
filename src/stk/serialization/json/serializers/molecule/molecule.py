from .utilities import atom_to_json, bond_to_json
from ..keys import InchiKey


class MoleculeJsonizer:
    """

    """

    def __init__(
        self,
        molecule_keys=(InchiKey(), ),
    ):
        """

        """

        self._molecule_keys = molecule_keys

    def to_json(self, molecule):
        """
        Return a JSON representation of `molecule`.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            The molecule to serialize.

        Returns
        -------
        :class:`dict`
            A JSON representation of `molecule`.

        """

        json = {
            'atoms': tuple(map(atom_to_json, molecule.get_atoms())),
            'bonds': tuple(map(bond_to_json, molecule.get_bonds())),
        }
        for molecule_key in self._molecule_keys:
            json[molecule_key.get_name()] = (
                molecule_key.get_key(molecule)
            )
        return json
