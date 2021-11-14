"""
Functional Group
================

.. toctree::
    :maxdepth: 2

    Alcohol <stk.molecular.functional_groups.functional_groups.alcohol>
    Aldehyde <\
stk.molecular.functional_groups.functional_groups.aldehyde\
>
    Alkene <stk.molecular.functional_groups.functional_groups.alkene>
    Alkyne <stk.molecular.functional_groups.functional_groups.alkyne>
    Amide <stk.molecular.functional_groups.functional_groups.amide>
    Boronic Acid <\
stk.molecular.functional_groups.functional_groups.boronic_acid\
>
    Bromo <stk.molecular.functional_groups.functional_groups.bromo>
    Carboxylic Acid <\
stk.molecular.functional_groups.functional_groups.carboxylic_acid\
>
    Dibromo <stk.molecular.functional_groups.functional_groups.dibromo>
    Difluoro <\
stk.molecular.functional_groups.functional_groups.difluoro\
>
    Diol <stk.molecular.functional_groups.functional_groups.diol>
    Fluoro <stk.molecular.functional_groups.functional_groups.fluoro>
    Generic Functional Group <\
stk.molecular.functional_groups.functional_groups.\
generic_functional_group\
>
    Iodo <stk.molecular.functional_groups.functional_groups.iodo>
    Primary Amino <\
stk.molecular.functional_groups.functional_groups.primary_amino\
>
    Ring Amine <\
stk.molecular.functional_groups.functional_groups.ring_amine\
>
    Secondary Amino <\
stk.molecular.functional_groups.functional_groups.secondary_amino\
>
    Single Atom <\
stk.molecular.functional_groups.functional_groups.single_atom\
>
    Thioacid <\
stk.molecular.functional_groups.functional_groups.thioacid\
>
    Thiol <stk.molecular.functional_groups.functional_groups.thiol>


Functional groups define which atoms of a :class:`.BuildingBlock` are
modified during :class:`.ConstructedMolecule` construction, and which
are used to position it.
The class of a :class:`.FunctionalGroup`
affects which :class:`.Reaction` can be used with it.
See the abstract base class :class:`.FunctionalGroup` for more
information.

"""


from __future__ import annotations


class FunctionalGroup:
    """
    An abstract base class for functional groups.

    It is used to give access to atoms of :class:`.BuildingBlock`
    molecules which are modified during
    :class:`.ConstructedMolecule` construction, as well as specify
    which atoms of the building block should be used for positioning.

    *Should I use* :meth:`.with_ids` *or* :meth:`.with_atoms` *?*

    That depends on your use case, however, it is generally better to
    default to :meth:`.with_ids` unless you need to actually change
    the atoms held by the functional group. This is because
    :meth:`.with_ids` preserves the most-derived type of the functional
    group, while :meth:`.with_atoms` does not. To give an example

    .. testcode:: with-atoms-vs-with-ids

        import stk

        bromo = stk.Bromo(
            bromine=stk.Br(0),
            atom=stk.C(1),
            bonders=(stk.C(1), ),
            deleters=(stk.Br(0), ),
        )

        bromo2 = bromo.with_ids({
            0: 10,
            1: 100,
        })
        # bromo2 is still a Bromo functional group.
        assert isinstance(bromo2, stk.Bromo)

        not_bromo = bromo.with_atoms({
            0: stk.Br(10),
            1: stk.C(100),
        })
        # not_bromo is not a Bromo functional gorup.
        assert not isinstance(not_bromo, stk.Bromo)
        # However, it is still an instance of FunctionalGroup,
        assert isinstance(not_bromo, stk.FunctionalGroup)
        # and of GenericFunctionalGroup
        assert isinstance(not_bromo, stk.GenericFunctionalGroup)

    The reason that :meth:`.with_atoms` does not produce a
    :class:`.Bromo` instance is to avoid the following pitfall

    .. code-block:: python

        pitfall = bromo.with_atoms({
            0: stk.F(10),
            1: stk.C(100),
        })
        # If with_atoms() returned a Bromo instance then you could
        # call get_bromine() on it, but it would hold a F atom!
        this_is_a_fluorine = pitfall.get_bromine()

    *Why would I want to implement a new subclass?*

    The most common reason you would want to implement a new
    :class:`.FunctionalGroup` subclass, is because you want to
    customize the construction of a :class:`.ConstructedMolecule`.
    Specifically, you want to modify a specific set of atoms in a
    :class:`.BuildingBlock` when doing construction, and you want
    to modify them in a specific way. You will usually accompany the
    creation of the new :class:`.FunctionalGroup` subclass with the
    creation of a new :class:`.Reaction` subclass, which will perform
    the custom modification on the atoms held by your new
    :class:`.FunctionalGroup`
    subclass. Finally, a new :class:`.ReactionFactory` will also be
    created, so that your :class:`.Reaction` subclass instances
    actually get made during construction. Finally, you will pass an
    instance of your :class:`.ReactionFactory` subclass to the
    chosen :class:`.TopologyGraph` you want to make, for example
    :class:`~.polymer.linear.linear.Linear`, and your custom
    modification will take place.

    See Also:

        :mod:`.functional_group_factory`
            Used for automated creation of :class:`.FunctionalGroup`
            instances.

    Notes:

        You might notice that some of the methods of this abstract base
        class are implemented. This is purely for convenience when
        implementing subclasses. The implemented public methods are
        simply default implementations, which can be safely ignored
        or overridden, when implementing subclasses. Any private
        methods are implementation details of these default
        implementations.

    Examples:

        *Subclass Implementation*

        The source code of the subclasses, listed in
        :mod:`.functional_group`, can serve as good examples.

        *Changing the Atoms of a Functional Group*

        You want to substitute the atoms in the functional group for
        other atoms. You can do this by using :meth:`.with_atoms` to
        create a clone of the functional group, which holds the
        replacement atoms

        .. testcode:: changing-the-atoms-of-a-functional-group

            import stk

            c, n, h1, h2 = stk.C(0), stk.N(1), stk.H(2), stk.H(3)
            amine = stk.PrimaryAmino(
                nitrogen=n,
                hydrogen1=h1,
                hydrogen2=h2,
                atom=c,
                bonders=(n, ),
                deleters=(h1, h2),
            )

            n20 = stk.N(20)
            h100 = stk.H(100)

            # amine_clone is a clone of amine, except that instead of
            # holding n, amine_clone holds n20, and instead of holding
            # h1  amine_clone holds h100. amine_clone continues to hold
            # h2.
            amine_clone = amine.with_atoms({
                n.get_id(): n20,
                h1.get_id(): h100,
            })

        .. testcode:: changing-the-atoms-of-a-functional-group
            :hide:

            _atoms = set(amine.get_atom_ids())
            assert n.get_id() in _atoms
            assert n20.get_id() not in _atoms
            assert h1.get_id() in _atoms
            assert h100.get_id()  not in _atoms

            _clone_atoms = set(amine_clone.get_atom_ids())
            assert n.get_id() not in _clone_atoms
            assert n20.get_id() in _clone_atoms
            assert h1.get_id() not in _clone_atoms
            assert h100.get_id() in _clone_atoms

            _bonders = set(amine.get_bonder_ids())
            assert n.get_id() in _bonders
            assert n20.get_id() not in _bonders

            _clone_bonders = set(amine_clone.get_bonder_ids())
            assert n.get_id() not in _clone_bonders
            assert n20.get_id() in _clone_bonders

            _deleters = set(amine.get_deleter_ids())
            assert h1.get_id() in _deleters
            assert h2.get_id() in _deleters
            assert h100.get_id() not in _deleters

            _clone_deleters = set(amine_clone.get_deleter_ids())
            assert h1.get_id() not in _clone_deleters
            assert h2.get_id() in _clone_deleters
            assert h100.get_id() in _clone_deleters

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
            of the :class:`.BuildingBlock`. See
            :meth:`.BuildingBlock.get_core_atom_ids`.

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
        Yield the ids of *placer* atoms.

        *Placer* atoms are those, which should be used to calculate
        the position of the functional group.

        Yields
        ------
        :class:`int`
            The id of an :class:`.Atom`.

        """

        yield from (a.get_id() for a in self._placers)

    def get_core_atom_ids(self):
        """
        Yield the ids of core atoms held by the functional group.

        Yields
        ------
        :class:`int`
            The id of an :class:`.Atom`.

        See Also
        --------
        :meth:`.BuildingBlock.get_core_atom_ids`

        """

        yield from (a.get_id() for a in self._core_atoms)

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
            The clone. Has the same type as the original functional
            group.

        """

        # The clone needs to be downcasted.
        return FunctionalGroup(
            atoms=tuple(
                atom_map.get(a.get_id(), a) for a in self._atoms
            ),
            placers=tuple(
                atom_map.get(a.get_id(), a) for a in self._placers
            ),
            core_atoms=tuple(
                atom_map.get(a.get_id(), a) for a in self._core_atoms
            ),
        )

    def _with_ids(self, id_map: dict[int, int]) -> FunctionalGroup:
        self._atoms = tuple(
            atom.with_id(
                id=id_map.get(
                    atom.get_id(),
                    atom.get_id(),
                ),
            ) for atom in self._atoms
        )
        self._placers = tuple(
            atom.with_id(
                id=id_map.get(
                    atom.get_id(),
                    atom.get_id(),
                ),
            ) for atom in self._placers
        )
        self._core_atoms = tuple(
            atom.with_id(
                id=id_map.get(
                    atom.get_id(),
                    atom.get_id(),
                ),
            ) for atom in self._core_atoms
        )
        return self

    def with_ids(self, id_map: dict[int, int]) -> FunctionalGroup:
        """
        Return a clone holding different atom ids.

        Parameters:

            id_map:
                Maps the id of an atom in the functional group to the
                new id the clone should hold. If the id of an atom in
                the functional group is not found in `id_map`, the
                atom will not be replaced in the clone.

        Returns:

            The clone.

        """

        return self.clone()._with_ids(id_map)

    def clone(self):
        """
        Return a clone.

        Returns
        -------
        :class:`.FunctionalGroup`
            A clone. Has the same type as the original functional
            group.

        """

        clone = self.__class__.__new__(self.__class__)
        FunctionalGroup.__init__(
            self=clone,
            atoms=self._atoms,
            placers=self._placers,
            core_atoms=self._core_atoms,
        )
        return clone

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'atoms={self._atoms}, '
            f'placers={self._placers}, '
            f'core_atoms={self._core_atoms}'
            ')'
        )
