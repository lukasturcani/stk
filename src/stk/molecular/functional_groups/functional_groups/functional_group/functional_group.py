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
modified by during :class:`.ConstructedMolecule` construction. They
also define which atoms are used to position the
:class:`.BuildingBlock` molecules
during construction. The class of a :class:`.FunctionalGroup`
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

        raise NotImplementedError()

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
            amine = stk.Amine(
                atoms=(n, h1, h2),
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

        raise NotImplementedError()

    def clone(self):
        """
        Return a clone.

        Returns
        -------
        :class:`.FunctionalGroup`
            A clone.

        """

        raise NotImplementedError()
