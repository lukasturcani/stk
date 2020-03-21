"""
Molecule Key
============

#. :class:`.Inchi`
#. :class:`.InchiKey`

"""


class MoleculeKey:
    """
    An abstract base class for molecular keys.

    Notes
    -----
    You might notice that the public methods of this abstract base
    class are implemented. This is purely for convenience when
    implementing subclasses. The implemented public methods are
    simply default implementations, which can be safely ignored or
    overridden, when implementing subclasses.

    Examples
    --------
    *Subclass Implementation*

    The source code of any of the subclasses, listed in
    :mod:`.molecule_key`, can serve as good examples.

    *Usage*

    Apart from using the subclasses provided, a :class:`.MoleculeKey`
    can be used directly , if you don't feel like writing a
    subclass

    .. code-block:: python

        import stk

        # Create a MoleculeKey instance with a custom key method.
        get_num_atoms = stk.MoleculeKey(
            name='num_atoms',
            key=lambda molecule: molecule.get_num_atoms(),
        )

        # Use the MoleculeKey instance.
        jsonizer = stk.MoleculeJsonizer(
            molecule_keys=(get_num_atoms, ),
        )
        # Get the JSON representation of a molecule.
        json = jsonizer.to_json(stk.BuildingBlock('NCCN'))

    """

    def __init__(self, name, key):
        """
        Initialize a :class:`.MoleculeKey` instance.

        Parameters
        ----------
        name : :class:`str`
            The name of the key.

        key : :class:`callable`
            Takes a single parameter, `molecule`, and returns the
            key to use for that molecule. The value passed to the
            parameter must be a :class:`.Molecule` instance.

        """

        self._name = name
        self._key = key

    def get_name(self):
        """
        Get the name of the key.

        Returns
        -------
        :class:`str`
            The name of the key.

        """

        return self._name

    def get_key(self, molecule):
        """
        Get the key of `molecule`.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            The molecule for which a key is needed.

        Returns
        -------
        :class:`object`
            The key of `molecule`.

        """

        return self._key(molecule)
