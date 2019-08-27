"""
Defines :class:`.Reactor`.

.. _`adding complex reactions`:

Extending stk: Adding complex reactions.
----------------------------------------

See :class:`.Reactor`.

"""

from collections import Counter
import numpy as np
from scipy.spatial.distance import euclidean

from . import elements
from .bonds import Bond


class _ReactionKey:
    """
    Represents a hash value of a collection of functional group names.

    The functional group names indicate which functional groups are
    involved in a reaction. The order in which the names are provided
    does not change the hash value. The quantity of each functional
    group name does change the hash value. Each name can be provided
    multiple times. Adding the same name again will change the hash
    value.

    Examples
    --------
    The key is constant with respect to permutations

    .. code-block:: python

        key1 = _ReactionKey('amine', 'aldehyde', 'amine')
        key2 = _ReactionKey('amine', 'amine', 'aldehyde')
        key3 = _ReactionKey('aldehyde', 'amine', 'amine')

        key1 == key2  # This is True.
        key1 == key3  # This is True.
        key2 == key3  # This is True.

    The key changes if the number of functional groups changes

    .. code-block:: python

        key4 = _ReactionKey('amine', 'aldehyde')
        key4 != key1  # This is True.


    """

    __slots__ = ['_key']

    def __init__(self, *fg_names):
        """
        Initialize a :class:`_ReactionKey`.

        Parameters
        ----------
        *fg_names : :class:`str`
            The names of functional groups involved in a reaction. If
            multiples of a functional group are present in a reaction,
            the name of the functional group must be present multiple
            times.

        """

        c = Counter(fg_names)
        # _key is the hash value for the reaction. It should change
        # if the different numbers of functional groups are present
        # in the reaction, but it should not change if the order in
        # which the functional groups are given changes.
        self._key = tuple(
            sorted((key, value) for key, value in c.items())
        )

    def __eq__(self, other):
        return self._key == other._key

    def __hash__(self):
        return hash(self._key)

    def __repr__(self):
        fg_names = ', '.join(
            repr(name)
            for name, count in self._key
            for i in range(count)
        )
        return f'_ReactionKey({fg_names})'

    def __str__(self):
        return repr(self)


class Reactor:
    """
    Performs reactions between functional groups of a molecule.

    The :class:`Reactor` is initialized with the
    :class:`.ConstructedMolecule` instance on which it will perform
    reactions.

    .. code-block:: python

        reactor = Reactor(mol)

    The :class:`Reactor` registers which functional groups need to have
    reactions performed with :meth:`add_reaction`

    .. code-block:: python

        fg1, fg2, fg3, fg4, fg5, fg6, fg7 = mol.func_groups

        # Register a reaction between the atoms in fg1 and fg2.
        reactor.add_reaction(fg1, fg2)

        # Register a reaction between atoms in fg3 and fg4.
        reactor.add_reaction(fg3, fg4)

        # Some reactions can take multiple functional groups.
        # You can put in as many functional groups as you like, given
        # an appropriate reaction is defined.
        reactor.add_reaction(fg5, fg6, fg7)

    Once all the reactions have been registered, they are perfomed
    in a single step

    .. code-block:: python

        reactor.finalize()

    An obvious question given this tutorial, is what reaction does
    :meth:`react` carry out? This is documented by :meth:`react`.
    However, :meth:`react` in most cases, will carry out a default
    reaction, which adds a bond between the bonder atoms of two
    functional groups. The bond order of the added bond is single by
    default but can be modified by editing :attr:`_bond_orders`. Here
    you will specify the :class:`_ReactionKey` for a reaction and what
    bond order you want that reaction to use.

    For some reactions you may wish to forgo the default reaction and
    do something more complex. This is neccessary because some
    reactions cannot be described by the simple combination of adding a
    bond while deleting some existing atoms. For example,
    consider the aldol reaction:

        CH3C(=O)CH3 + CH3C(=O)CH3 --> CH3(=O)CH2C(OH)(CH3)CH3

    Here a ketone is converted into an alcohol. If you wish to
    support a complex reaction, add it as a method within this class.
    The method will need take some :class:`FunctionalGroup` instances
    as arguments. These are the functional groups which react.

    Once the method is defined, :attr:`_custom_reactions` needs to
    be updated.

    """

    # When the default reaction is added by add_reaction(),
    # if the bond added between the two functional groups is not
    # single, the desired bond order should be placed in this
    # dictionary. The dictionary maps the reaction's
    # _ReactionKey to the desired bond order.
    _bond_orders = {
        _ReactionKey('amine', 'aldehyde'): 2,
        _ReactionKey('amide', 'aldehyde'): 2,
        _ReactionKey('nitrile', 'aldehyde'): 2,
        _ReactionKey('amide', 'amine'): 2,
        _ReactionKey('terminal_alkene', 'terminal_alkene'): 2,
        _ReactionKey('alkyne2', 'alkyne2'): 3
    }

    def __init__(self, mol):
        """
        Initialize a :class:`Reactor`.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The molecule on which the reactor adds and removes atoms
            and bonds.

        """

        # The molecule from which atoms and bonds are added and
        # removed.
        self._mol = mol
        # The ids of atoms which are to be removed.
        self._deleter_ids = set()
        # The ids of bonds which are to be removed.
        self._deleter_bonds = set()

        # Maps a _ReactionKey for a given reaction to a custom
        # method which carries out the reaction. This means that
        # add_reaction will use that method for carrying out that
        # reaction instead.
        self._custom_reactions = {

            _ReactionKey('boronic_acid', 'diol'):
                self._react_boronic_acid_with_diol,

            _ReactionKey('diol', 'difluorene'):
                self._react_diol_with_dihalogen,

            _ReactionKey('diol', 'dibromine'):
                self._react_diol_with_dihalogen,

            _ReactionKey('ring_amine', 'ring_amine'):
                self._react_ring_amine_with_ring_amine

        }

    def add_reaction(self, func_groups, periodicity):
        """
        Create bonds between functional groups.

        This function first looks at the functional groups provided via
        the `func_groups` argument and checks which functional groups
        are involved in the reaction. If the functional groups are
        handled by one of the custom reactions specified in
        :attr:`_custom_reactions` then that function is executed.

        In all other cases the function is assumed to have received two
        functional groups to react. In these functional groups, the
        bonder atoms have a bond added. The bond is single,
        unless otherwise specified in :attr:`_bond_orders`.

        Parameters
        ----------
        func_groups : :class:`tuple` of :class:`.FunctionalGroup`
            The functional groups to react.

        periodicity : :class:`tuple` of :class:`int`
            Specifies the periodicity of the bonds added by the
            reaction, which bridge the `func_groups`. See
            :attr:`.Bond.periodicity`.

        Returns
        -------
        None : :class:`NoneType`

        """

        names = (fg.fg_type.name for fg in func_groups)
        reaction_key = _ReactionKey(*names)
        reaction = self._custom_reactions.get(
            reaction_key,
            self._react_any
        )
        return reaction(reaction_key, func_groups, periodicity)

    def finalize(self):
        """
        Finish performing reactions.

        Returns
        -------
        None : :class:`NoneType`

        """

        self._mol.atoms = [
            atom for atom in self._mol.atoms
            if atom.id not in self._deleter_ids
        ]
        self._mol.bonds = [
            bond for bond in self._mol.bonds
            if bond.atom1.id not in self._deleter_ids
            and bond.atom2.id not in self._deleter_ids
        ]
        self._mol._position_matrix = [
            row for i, row in enumerate(self._mol._position_matrix)
            if i not in self._deleter_ids
        ]

    def _remove_deleters(self, func_groups):
        """
        Remove deleter atoms from `fgs`.

        Parameters
        ----------
        func_groups : :class:`list` of :class:`.FunctionalGroup`
            The functional groups from which deleter atoms should be
            removed.

        Returns
        -------
        None : :class:`NoneType`

        """

        for fg in func_groups:
            self._deleter_ids.update(fg.get_deleter_ids())
            fg.atoms = tuple(
                a for a in fg.atoms if a.id not in self._deleter_ids
            )
            fg.deleters = ()

    def _react_any(self, reaction_key, func_groups, periodicity):
        """
        Create bonds between functional groups.

        Parameters
        ----------
        reaction_key : :class:`._ReactionKey`
            The key for the reaction.

        func_groups : :class:`list` of :class:`.FunctionalGroup`
            The functional groups from which deleter atoms should be
            removed.

        periodicity : :class:`tuple` of :class:`int`
            Specifies the periodicity of the bonds added by the
            reaction, which bridge the `func_groups`. See
            :attr:`.Bond.periodicity`.

        Returns
        -------
        None : :class:`NoneType`

        """

        self._remove_deleters(func_groups)

        fg1, fg2 = func_groups
        bond_order = self._bond_orders.get(reaction_key, 1)
        bond = Bond(
            atom1=fg1.bonders[0],
            atom2=fg2.bonders[0],
            order=bond_order,
            periodicity=periodicity
        )
        self._mol.bonds.append(bond)
        self._mol.construction_bonds.append(bond)

    def _react_diol_with_dihalogen(
        self,
        reaction_key,
        func_groups,
        periodicity
    ):
        """
        Create bonds between functional groups.

        Parameters
        ----------
        reaction_key : :class:`._ReactionKey`
            The key for the reaction.

        func_groups : :class:`list` of :class:`.FunctionalGroup`
            The functional groups from which deleter atoms should be
            removed.

        periodicity : :class:`tuple` of :class:`int`
            Specifies the periodicity of the bonds added by the
            reaction, which bridge the `func_groups`. See
            :attr:`.Bond.periodicity`.

        Returns
        -------
        None : :class:`NoneType`

        """

        self._remove_deleters(func_groups)

        fg1, fg2 = func_groups
        diol = fg1 if fg1.fg_type.name == 'diol' else fg2
        dihalogen = fg2 if diol is fg1 else fg1

        distances = []
        for carbon in dihalogen.get_bonder_ids():
            c_coord = self._mol._position_matrix[carbon]
            for oxygen in diol.get_bonder_ids():
                o_coord = self._mol._position_matrix[oxygen]
                d = euclidean(c_coord, o_coord)
                distances.append((d, carbon, oxygen))
        distances.sort()

        deduped_pairs = []
        seen_o, seen_c = set(), set()
        for d, c, o in distances:
            if c not in seen_c and o not in seen_o:
                deduped_pairs.append((c, o))
                seen_c.add(c)
                seen_o.add(o)

        (c1, o1), (c2, o2), *_ = deduped_pairs
        assert c1 != c2 and o1 != o2

        bond1 = Bond(
            atom1=self._mol.atoms[c1],
            atom2=self._mol.atoms[o1],
            order=1,
            periodicity=periodicity
        )
        bond2 = Bond(
            atom1=self._mol.atoms[c2],
            atom2=self._mol.atoms[o2],
            order=1,
            periodicity=periodicity
        )
        self._mol.bonds.append(bond1)
        self._mol.construction_bonds.append(bond1)
        self._mol.bonds.append(bond2)
        self._mol.construction_bonds.append(bond2)

    def _react_boronic_acid_with_diol(
        self,
        reaction_key,
        func_groups,
        periodicity
    ):
        """
        Create bonds between functional groups.

        Parameters
        ----------
        reaction_key : :class:`._ReactionKey`
            The key for the reaction.

        func_groups : :class:`list` of :class:`.FunctionalGroup`
            The functional groups from which deleter atoms should be
            removed.

        periodicity : :class:`tuple` of :class:`int`
            Specifies the periodicity of the bonds added by the
            reaction, which bridge the `func_groups`. See
            :attr:`.Bond.periodicity`.

        Returns
        -------
        None : :class:`NoneType`

        """

        self._remove_deleters(func_groups)
        fg1, fg2 = func_groups
        boron = fg1 if fg1.fg_type.name == 'boronic_acid' else fg2
        diol = fg2 if boron is fg1 else fg1

        boron_atom = boron.bonders[0]
        bond1 = Bond(boron_atom, diol.bonders[0], 1, periodicity)
        self._mol.bonds.append(bond1)
        self._mol.construction_bonds.append(bond1)

        bond2 = Bond(boron_atom, diol.bonders[1], 1, periodicity)
        self._mol.bonds.append(bond2)
        self._mol.construction_bonds.append(bond2)

    def _react_ring_amine_with_ring_amine(
        self,
        reaction_key,
        func_groups,
        periodicity
    ):
        """
        Creates bonds between functional groups.

        Parameters
        ----------
        reaction_key : :class:`._ReactionKey`
            The key for the reaction.

        func_groups : :class:`list` of :class:`.FunctionalGroup`
            The functional groups from which deleter atoms should be
            removed.

        periodicity : :class:`tuple` of :class:`int`
            Specifies the periodicity of the bonds added by the
            reaction, which bridge the `func_groups`. See
            :attr:`.Bond.periodicity`.

        Returns
        -------
        None : :class:`NoneType`

        """

        self._remove_deleters(func_groups)
        fg1, fg2 = func_groups
        c1 = next(a for a in fg1.bonders if a.atomic_number == 6)
        n1 = next(a for a in fg1.bonders if a.atomic_number == 7)

        c2 = next(a for a in fg2.bonders if a.atomic_number == 6)
        n2 = next(a for a in fg2.bonders if a.atomic_number == 7)

        n1_coord = np.array(self._mol._position_matrix[n1.id])
        n2_coord = np.array(self._mol._position_matrix[n2.id])
        c1_coord = np.array(self._mol._position_matrix[c1.id])
        c2_coord = np.array(self._mol._position_matrix[c2.id])

        n_joiner = elements.C(len(self._mol.atoms))
        self._mol.atoms.append(n_joiner)
        n_joiner_coord = (n1_coord + n2_coord) / 2
        self._mol._position_matrix.append(n_joiner_coord)

        nh1 = elements.H(len(self._mol.atoms))
        self._mol.atoms.append(nh1)
        nh1_coord = n_joiner_coord + np.array([0, 0, 1])
        self._mol._position_matrix.append(nh1_coord)

        nh2 = elements.H(len(self._mol.atoms))
        self._mol.atoms.append(nh2)
        nh2_coord = n_joiner_coord + np.array([0, 0, -1])
        self._mol._position_matrix.append(nh2_coord)

        nc_joiner1 = elements.C(len(self._mol.atoms))
        self._mol.atoms.append(nc_joiner1)
        nc_joiner1_coord = (c1_coord + n2_coord) / 2
        self._mol._position_matrix.append(nc_joiner1_coord)

        nc1h1 = elements.H(len(self._mol.atoms))
        self._mol.atoms.append(nc1h1)
        nc1h1_coord = nc_joiner1_coord + np.array([0, 0, 1])
        self._mol._position_matrix.append(nc1h1_coord)

        nc1h2 = elements.H(len(self._mol.atoms))
        self._mol.atoms.append(nc1h2)
        nc1h2_coord = nc_joiner1_coord + np.array([0, 0, -1])
        self._mol._position_matrix.append(nc1h2_coord)

        nc_joiner2 = elements.C(len(self._mol.atoms))
        self._mol.atoms.append(nc_joiner2)
        nc_joiner2_coord = (c2_coord + n1_coord) / 2
        self._mol._position_matrix.append(nc_joiner2_coord)

        nc2h1 = elements.H(len(self._mol.atoms))
        self._mol.atoms.append(nc2h1)
        nc2h1_coord = nc_joiner2_coord + np.array([0, 0, 1])
        self._mol._position_matrix.append(nc2h1_coord)

        nc2h2 = elements.H(len(self._mol.atoms))
        self._mol.atoms.append(nc2h2)
        nc2h2_coord = nc_joiner2_coord + np.array([0, 0, -1])
        self._mol._position_matrix.append(nc2h2_coord)

        bonds = (
            Bond(n1, n_joiner, 1),
            Bond(n_joiner, n2, 1, periodicity),
            Bond(n_joiner, nh1, 1),
            Bond(n_joiner, nh2, 1),
            Bond(c1, nc_joiner1, 1),
            Bond(nc_joiner1, n2, 1, periodicity),
            Bond(nc_joiner1, nc1h1, 1),
            Bond(nc_joiner1, nc1h2, 1),
            Bond(nc_joiner2, c2, 1, periodicity),
            Bond(n1, nc_joiner2, 1),
            Bond(nc_joiner2, nc2h1, 1),
            Bond(nc_joiner2, nc2h2, 1)
        )
        for bond in bonds:
            self._mol.bonds.append(bond)
            self._mol.construction_bonds.append(bond)
