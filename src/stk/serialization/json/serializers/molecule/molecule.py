"""
Molecule JSONizer
=================

"""

from stk.molecular import InchiKey
from .utilities import atom_to_json, bond_to_json


class MoleculeJsonizer:
    """
    Creates JSON representations of :class:`.Molecule` instances.

    Examples
    --------

    """

    def __init__(
        self,
        key_makers=(InchiKey(), ),
    ):
        """
        Initialize a :class:`.MoleculeJsonizer` instance.

        Parameters
        ----------
        key_makers : :class:`tuple` of \
                :class:`.MoleculeKeyMaker`
            Used to make the keys of molecules, which should be
            included in their JSON representations.

        """

        self._key_makers = key_makers

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
        for key_maker in self._key_makers:
            json[key_maker.get_key_name()] = (
                key_maker.get_key(molecule)
            )
        return json
