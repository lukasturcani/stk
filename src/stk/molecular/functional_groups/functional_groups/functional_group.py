"""
Functional Group
================

#. :class:`.Alcohol`
#. :class:`.Aldehyde`
#. :class:`.Alkene`
#. :class:`.Alkyne`
#. :class:`.Amide`
#. :class:`.BoronicAcid`
#. :class:`.Bromo`
#. :class:`.CarboxylicAcid`
#. :class:`.Dibromo`
#. :class:`.Difluoro`
#. :class:`.Diol`
#. :class:`.Fluoro`
#. :class:`.GenericFunctionalGroup`
#. :class:`.Iodo`
#. :class:`.PrimaryAmino`
#. :class:`.RingAmine`
#. :class:`.SecondaryAmino`
#. :class:`.Thioacid`
#. :class:`.Thiol`

Functional groups define which atoms of a :class:`.BuildingBlock` are
modified by during :class:`.ConstructedMolecule` construction.
The class of a :class:`.FunctionalGroup`
affects which :class:`.Reaction` can be used with it.
See the abstract base class :class:`.FunctionalGroup` for more
information.

"""


class FunctionalGroup:
    """
    An abstract base class for functional groups.

    It is used to give access to atoms of :class:`.BuildingBlock`
    molecules which are modified during
    :class:`.ConstructedMolecule` construction.

    """

    def __init__(self, atoms, placers, core_atoms):
        """
        Initialize a :class:`.FunctionalGroup`.

        Parameters
        ----------
        atoms : :class:`tuple` of :class:`.Atom`
            The atoms in the functional group.

        placers : :class:`tuple` of :class:`.Atom`
            The atoms used to calculate the position of the functional
            group.

        core_atoms : :class:`tuple` of :class:`.Atom`
            The atoms of the functional group which also form the core
            of the :class:`.BuildingBlock`.

        """

        self._atoms = atoms
        self._placers = placers
        self._core_atoms = core_atoms

    def get_atoms(self):
        """
        Yield all the atoms in the functional group.

        Yields
        ------
        :class:`.Atom`
            An atom in the functional group.

        """

        yield from self._atoms

    def get_atom_ids(self):
        """
        Yield the ids of all atoms in the functional group.

        Yields
        ------
        :class:`int`
            The id of an :class:`.Atom`.

        """

        yield from (a.get_id() for a in self._atoms)

    def get_placer_ids(self):
        """
        Yield the ids of atoms used for position calculation.

        The ids of atoms yielded by this method should be used to
        calculate the position of the functional group.

        Yields
        ------
        :class:`int`
            The id of an :class:`.Atom`.

        """

        yield from (a.get_id() for a in self._placers)

    def get_core_ids(self):
        """
        Yield the ids of core atoms.

        Yields
        ------
        :class:`int`
            The id of an :class:`.Atom`.

        """

        yield from (a.get_id() for a in self._core_atoms)

    def _with_atoms(self, atom_map):
        """
        Modify the functional group.

        """

        self._atoms = tuple(
            atom_map.get(a.get_id(), a) for a in self._atoms
        )
        self._placers = tuple(
            atom_map.get(a.get_id(), a) for a in self._placers
        )
        return self

    def with_atoms(self, atom_map):
        """
        Return a clone holding different atoms.

        Parameters
        ----------
        atom_map : :class:`dict`
            Maps the id of an atom in the functional group to the new
            atom the clone should hold. If the id of an atom in the
            functional group is not found in `atom_map`, the atom
            will not be replaced in the clone.

        Returns
        -------
        :class:`.FunctionalGroup`
            The clone.

        Examples
        --------
        .. code-block:: python

            import stk

            n, h1, h2 = stk.N(0), stk.H(1), stk.H(2)
            amine = stk.PrimaryAmino(
                nitrogen=n,
                hydrogen1=h1,
                hydrogen2=h2,
                bonders=(n, ),
                deleters=(h1, h2),
            )

            n20 = stk.N(20)
            h100 = stk.H(100)

            # fg_clone is a clone of fg, except that instead of holding
            # n, fg_clone holds n20, and instead of holding h1
            # fg_clone holds h100. fg_clone continues to hold h2.
            fg_clone = fg.with_atoms({
                n.get_id(): n20,
                h1.get_id(): h100,
            })


        """

        return self.clone()._with_atoms(atom_map)

    def clone(self):
        """
        Return a clone.

        Returns
        -------
        :class:`.FunctionalGroup`
            A clone.

        """

        clone = self.__class__.__new__(self.__class__)
        FunctionalGroup.__init__(clone, self._atoms, self._placers)
        return clone

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'atoms={self._atoms}, '
            f'placers={self._placers}'
            ')'
        )
