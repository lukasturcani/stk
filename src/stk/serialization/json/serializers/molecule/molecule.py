from .utilities import atom_to_json, bond_to_json


class MoleculeJsonizer:
    """

    """

    def __init__(self, molecule_key):
        self._molecule_key = molecule_key

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

        return {
            self._molecule_key.get_name():
                self._molecule_key.get_key(molecule),
            'atoms': tuple(map(atom_to_json, molecule.get_atoms())),
            'bonds': tuple(map(bond_to_json, molecule.get_bonds())),
        }
