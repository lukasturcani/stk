"""
Molecule Key Maker
==================

.. toctree::
    :maxdepth: 2

    InChI <stk.molecular.key_makers.inchi>
    InChIKey <stk.molecular.key_makers.inchi_key>
    SMILES <stk.molecular.key_makers.smiles>

"""


class MoleculeKeyMaker:
    """
    An abstract base class for making :class:`.Molecule` keys.

    Keys are used in :mod:`stk` to determine if two molecules are
    duplicates of each other.

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
    :mod:`molecule_key_maker <.key_makers.molecule>`, can
    serve as good examples.

    *Creating a New Key Maker Directly*

    Apart from using the subclasses provided, a
    :class:`.MoleculeKeyMaker` can be used directly , if you don't
    feel like writing a subclass

    .. testcode:: creating-a-new-key-maker-directly

        import stk

        # Create a MoleculeKeyMaker instance with a custom get_key
        # method.
        get_num_atoms = stk.MoleculeKeyMaker(
            key_name='num_atoms',
            get_key=lambda molecule: molecule.get_num_atoms(),
        )

        # Use the MoleculeKeyMaker instance.
        jsonizer = stk.MoleculeJsonizer(
            key_makers=(get_num_atoms, ),
        )
        # Get the JSON representation of a molecule.
        json = jsonizer.to_json(stk.BuildingBlock('NCCN'))

    .. testcode:: creating-a-new-key-maker-directly
        :hide:

        assert json['molecule']['num_atoms'] == 12
        assert json['matrix']['num_atoms'] == 12

    """

    def __init__(self, key_name, get_key):
        """
        Initialize a :class:`.MoleculeKeyMaker` instance.

        Parameters
        ----------
        key_name : :class:`str`
            The name of the key.

        get_key : :class:`callable`
            Takes a single parameter, `molecule`, and returns the
            key to use for that molecule. The value passed to the
            parameter must be a :class:`.Molecule` instance.

        """

        self._key_name = key_name
        self._get_key = get_key

    def get_key_name(self):
        """
        Get the name of the key.

        Returns
        -------
        :class:`str`
            The name of the key.

        """

        return self._key_name

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

        return self._get_key(molecule)

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._key_name!r}, {self._get_key!r}'
            ')'
        )
