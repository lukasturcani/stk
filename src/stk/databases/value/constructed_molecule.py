"""
Constructed Molecule Value Database
===================================

#. :class:`.MoleculeValueMongoDb`
#. :class:`.ConstructedMoleculeValueMongoDb`
#. :class:`.MoleculeValueDatabase`

"""


class ConstructedMoleculeValueDatabase:
    """
    Abstract base class for storing of molecular property values.

    Note that a :class:`.MoleculeValueDatabase` can be used anywhere a
    :class:`.ConstructedMoleculeValueDatabase` is required. However,
    the opposite is not true. If something requires a
    :class:`.MoleculeValueDatabase` you cannot use a
    :class:`.ConstructedMoleculeValueDatabase` in its place.

    Examples
    --------
    *Subclass Implementation*

    The source of any of the subclasses, listed in
    :mod:`constructed_molecule_value_database \
    <.databases.value.constructed_molecule>`,
    can serve as good examples.

    """

    def put(self, molecule, value):
        """
        Put a value into the database.

        Parameters
        ----------
        molecule : :class:`.ConstructedMolecule`
            The molecule which is associated with the `value`.

        value : :class:`object`
            Some value associated with `molecule`.

        Returns
        -------
        None : :class:`NoneType`

        """

        raise NotImplementedError()

    def get(self, molecule):
        """
        Get the stored value for `molecule`.

        Parameters
        ----------
        molecule : :class:`.ConstructedMolecule`
            The molecule, whose value is to be retrieved from the
            database.

        Returns
        -------
        :class:`object`
            The value associated with `molecule`.

        Raises
        ------
        :class:`KeyError`
            If `molecule` is not found in the database.

        """

        raise NotImplementedError()
