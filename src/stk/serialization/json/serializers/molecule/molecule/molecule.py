from .utilities import atom_to_json, bond_to_json


class _MoleculeSerializer:
    """
    Serializes :class:`.Molecule` instance to JSON.

    Notes
    -----
    This class is an implementation detail of
    :class:`.MoleculeSerializer`, use that class directly for your
    serialization needs.

    """

    # Even though this __init__() method does not do anything, keep it
    # for the docstring.
    def __init__(self):
        """
        Initialize a :class:`._MoleculeSerializer` instance.

        """

        return super().__init__()

    def serialize(self, molecule):
        """
        Serialize `molecule`.

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
            'atoms': tuple(map(atom_to_json, molecule.get_atoms())),
            'bonds': tuple(map(bond_to_json, molecule.get_bonds())),
        }
