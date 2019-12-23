class FunctionalGroup:
    """
    An abstract base class for functional groups.

    It is used to give access to atoms of :class:`.BuildingBlock`
    molecules which are modified during
    :class:`.ConstructedMolecule` construction.

    """

    def get_atoms(self):
        """
        Yield all the atoms in the functional group.

        Yields
        ------
        :class:`.Atom`
            An atom in the functional group.

        """

        raise NotImplementedError()

    def get_atom_ids(self):
        """
        Yield the ids of all atoms in the functional group.

        Yields
        ------
        :class:`int`
            The id of an :class:`.Atom`.

        """

        raise NotImplementedError()

    def get_bonders(self):
        """
        Yield bonder atoms in the functional group.

        These are atoms which have bonds added during
        :class:`.ConstructedMolecule` construction.

        Yields
        ------
        :class:`.Atom`
            A bonder atom.

        """

        raise NotImplementedError()

    def get_bonder_ids(self):
        """
        Yield the ids of bonder atoms.

        Yields
        ------
        :class:`int`
            The id of a bonder :class:`.Atom`.

        """

        raise NotImplementedError()

    def get_deleters(self):
        """
        Yield the deleter atoms in the functional group.

        These are atoms which are removed during
        :class:`.ConstructedMolecule` construction.

        Yields
        ------
        :class:`.Atom`
            A deleter atom.

        """

        raise NotImplementedError()

    def get_deleter_ids(self):
        """
        Yield the ids of deleter atoms.

        Yields
        -------
        :class:`int`
            The id of a deleter :class:`.Atom`.

        """

        raise NotImplementedError()

    def clone(self, atom_map=None):
        """
        Return a clone.

        Public attributes are inherited by the clone but private
        ones are not.

        Parameters
        ----------
        atom_map : :class:`dict`, optional
            If the clone should hold specific :class:`.Atom`
            instances, then a :class:`dict` should be provided, which
            maps atom ids in the current :class:`.FunctionalGroup` to
            the atoms which should be used in the clone.

        Returns
        -------
        :class:`FunctionalGroup`
            A clone.

        Examples
        --------
        .. code-block:: python

            import stk

            n, h1, h2 = stk.N(0), stk.H(1), stk.H(2)
            amine = stk.Amine(n, h1, h2)

            n20 = stk.N(20)
            h100 = stk.H(100)

            # fg_clone is a clone of fg, except that instead of holding
            # a clone of n, fg_clone holds n20, and instead of holding
            # a clone of h1 fg_clone holds h100. fg_clone does hold a
            # clone of h2.
            fg_clone = fg.clone({
                n.id: n20,
                h1.id: h100,
            })

        """

        raise NotImplementedError()
