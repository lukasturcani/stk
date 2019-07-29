"""
Defines tools for dealing with functional groups and their reactions.

See the documentation of :class:`.Reactor` to see how reactions between
functional groups are performed.

.. _`adding functional groups`:

Extending stk: Adding  more functional groups.
----------------------------------------------

If ``stk`` is to incorporate a new functional group, a new
:class:`FGType` instance should be added to
:data:`functional_groups`, which is defined in this module.

Adding a new :class:`FGType` instance to :data:`functional_groups` will
allow :meth:`.TopologyGraph.construct` to connect the functional group
to all others during construction. In most cases, nothing except adding
this instance should be necessary in order to incorporate new
functional groups.

Note that when adding SMARTS, if you want to make a SMARTS that targets
an atom in an environment, for example, a bromine connected to a
carbon::

    [$([Br][C])]

The atom you are targeting needs to be written first. The above SMARTS
works but::

    [$([C][Br])]

does not.

If a new functional group is to connect to another functional group
with a bond other than a single, the names of the functional groups
should be added to :attr:`Reactor.bond_orders`, along with the desired
bond order.

.. _`adding complex reactions`:

Extending stk: Adding complex reactions.
----------------------------------------

See :class:`Reactor`.

"""

from functools import partial
import numpy as np
from collections import Counter
import rdkit.Chem.AllChem as rdkit

from . import elements
from .bonds import Bond, PeriodicBond
from ..utilities import flatten


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
        Initialize a :class:`ReactionKey`.

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


class FGType:
    """
    Represents a :class:`.FunctionalGroup` type.

    Instances of this class are used to group the
    :attr:`~.Molecule.atoms` of a :class:`.Molecule` into
    :class:`.FunctionalGroup` instances.

    Attributes
    ----------
    name : :class:`str`
        A name for the :class:`.FGType`. Helps with identification.

    """

    __slots__ = ['name', '_func_group', '_bonders', '_deleters']

    def __init__(
        self,
        name,
        func_group_smarts,
        bonder_smarts,
        deleter_smarts
    ):
        """
        Initialize a :class:`FGType` instance.

        Parameters
        ----------
        name : :class:`str`
            A name for the :class:`.FGType`. Helps with identification.

        func_group_smarts : :class:`str`
            A SMARTS string which matches all atoms in a functional
            group.

        bonder_smarts : :class:`list` of :class:`str`
            A :class:`list` of SMARTS strings, each of which matches a
            single atom in a functional group. The matched atom is
            added to :class:`.FunctionalGroup.bonders`. A SMARTS
            string needs to be repeated if a single functional group
            has multiple equivalent atoms, each of which has bonds
            created during construction. The number of times the string
            needs to be repeated, is equal to number of atoms which
            need to be tagged.

        deleter_smarts : :class:`list` of :class:`str`
            A :class:`list` of SMARTS strings, each of which matches a
            single atom in a functional group. The matched atom is
            added to :class:`.FunctionalGroup.deleters`. A SMARTS
            string needs to be repeated if a single functional group
            has multiple equivalent atoms, each of which is deleted
            during construction. The number of times the string needs
            to be repeated, is equal to number of atoms which need to
            be tagged.

        """

        self.name = name
        self._func_group = rdkit.MolFromSmarts(func_group_smarts)
        self._bonders = [
            rdkit.MolFromSmarts(smarts) for smarts in bonder_smarts
        ]
        self._deleters = [
            rdkit.MolFromSmarts(smarts) for smarts in deleter_smarts
        ]

    def get_functional_groups(self, mol):
        """
        Yield the functional groups in `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule which is to have functional groups
            identified.

        Yields
        ------
        :class:`.FunctionalGroup`
            A :class:`.FunctionalGroup` instance for every matched
            functional group in the `mol`.

        Examples
        --------
        .. code-block:: python

            import stk

            mol = stk.BuildingBlock('NCCN')
            amine = stk.FGType(
                func_group_smarts='[N]([H])[H]',
                bonder_smarts=['[$([N]([H])[H])]'],
                deleter_smarts=['[$([H][N][H])]']*2
            )

            # Make mol use amine functional groups during construction
            # of ConstructedMolecule.
            mol.func_groups = tuple(amine.get_functional_groups(mol))

        """

        rdkit_mol = mol.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_mol)

        func_groups = mol.GetSubstructMatches(self._func_group)

        # All the bonder atoms, grouped by fg.
        bonders = [[] for i in range(len(func_groups))]

        for bonder in self._bonders:
            matches = set(flatten(
                rdkit_mol.GetSubstructMatches(bonder)
            ))

            for fg_id, fg in enumerate(func_groups):
                for atom_id in fg:
                    if atom_id in matches:
                        bonders[fg_id].append(atom_id)
                        break

        # All the deleter atoms, grouped by fg.
        deleters = [[] for i in range(len(func_groups))]
        for deleter in self._deleters:
            matches = set(flatten(
                rdkit_mol.GetSubstructMatches(deleter)
            ))

            for fg_id, fg in enumerate(func_groups):
                for atom_id in fg:
                    if atom_id in matches:
                        deleters[fg_id].append(atom_id)
                        break

        for atom_ids in zip(func_groups, bonders, deleters):
            fg, fg_bonders, fg_deleters = atom_ids
            yield FunctionalGroup(
                atoms=tuple(self.atoms[id_] for id_ in fg),
                bonders=tuple(self.atoms[id_] for id_ in fg_bonders),
                deleters=tuple(self.atoms[id_] for id_ in fg_deleters),
                fg_type=self
            )

    def __repr__(self):
        func_group_smarts = rdkit.MolToSmarts(self._func_group)
        bonder_smarts = [
            rdkit.MolToSmarts(mol) for mol in self._bonders
        ]
        deleter_smarts = [
            rdkit.MolToSmarts(mol) for mol in self._deleters
        ]
        return (
            f'FGType(\n'
            f'    name={self.name!r}\n'
            f'    func_group_smarts={func_group_smarts!r},\n'
            f'    bonder_smarts={bonder_smarts!r},\n'
            f'    deleter_smarts={deleter_smarts!r}\n'
            ')'
        )

    def __str__(self):
        return f'FGType({self.name})'


class FunctionalGroup:
    """
    Represents a functional group in a molecule.

    Instances of this class should only by made by using
    :class:`.FGType.get_functional_groups`.

    Attributes
    ----------
    atoms : :class:`tuple` of :class:`.Atom`
        The atoms in the functional group.

    bonders : :class:`tuple` of :class:`.Atom`
        The bonder atoms in the functional group. These are atoms which
        have bonds added during the construction of a
        :class:`.ConstructedMolecule`.

    deleters : :class:`tuple` of :class:`.Atom`
        The deleter atoms in the functional group. These are atoms
        which are deleted during construction of a
        :class:`.ConstructedMolecule`.

    fg_type : :class:`.FGType`
        The :class:`.FGType` instance which created the
        :class:`.FunctionalGroup`.

    """

    def __init__(self, atoms, bonders, deleters, fg_type=None):
        """
        Initialize a :class:`.FunctionalGroup`.

        Parameters
        ----------
        atoms : :class:`tuple` of :class:`.Atom`
            The atoms in the functional group.

        bonders : :class:`tuple` of :class:`.Atom`
            The bonder atoms in the functional group.

        deleters : :class:`tuple` of :class:`.Atom`
            The deleter atoms in the functional group.

        fg_type : :class:`.FGType`
            The :class:`.FGType` instance which created the
            :class:`.FunctionalGroup`.

        """

        self.atoms = atoms
        self.bonders = bonders
        self.deleters = deleters
        self.fg_type = fg_type

    def clone(self, atom_map=None):
        """
        Return a clone.

        Parameters
        ----------
        atom_map : :class:`dict`, optional
            If the clone should hold different :class:`.Atom`
            instances, then a :class:`dict` should be provided, which
            maps atoms in the current :class:`.FunctionalGroup` to the
            atoms which should be used in the clone. Only atoms which
            need to be remapped need to be present in the `atom_map`.

        Returns
        -------
        :class:`FunctionalGroup`
            A clone.

        """

        if atom_map is None:
            atom_map = {}

        return self.__class__(
            atoms=tuple(atom_map.get(a, a) for a in self.atoms),
            bonders=tuple(atom_map.get(a, a) for a in self.bonders),
            deleters=tuple(atom_map.get(a, a)for a in self.deleters),
            fg_type=self.fg_type
        )

    def get_atom_ids(self):
        """
        Get the ids of :attr:`atoms`.

        Returns
        -------
        :class:`generator` of :class:`int`
            The ids of :attr:`atoms`.

        """

        return (a.id for a in self.atoms)

    def get_bonder_ids(self):
        """
        Get the ids of :attr:`bonders`.

        Returns
        -------
        :class:`generator` of :class:`int`
            The ids of :attr:`bonders`.

        """
        return (a.id for a in self.bonders)

    def get_deleter_ids(self):
        """
        Get the ids of :attr:`deleters`.

        Returns
        -------
        :class:`generator` of :class:`int`
            The ids of :attr:`deleters`.

        """

        return (a.id for a in self.deleters)

    def __repr__(self):
        atoms = list(self.atoms)
        if len(atoms) == 1:
            atoms.append('')
        atoms = ', '.join(repr(atom) if atom else '' for atom in atoms)

        bonders = list(self.bonders)
        if len(bonders) == 1:
            bonders.append('')
        bonders = ', '.join(
            repr(atom) if atom else '' for atom in bonders
        )

        deleters = list(self.deleters)
        if len(deleters) == 1:
            deleters.append('')
        deleters = ', '.join(
            repr(atom) if atom else '' for atom in deleters
        )

        return (
            f'FunctionalGroup(\n'
            f'    atoms=( {atoms} ), \n'
            f'    bonders=( {bonders} ), \n'
            f'    deleters=( {deleters} ), \n'
            f'    fg_type={self.fg_type.name}\n'
            ')'
        )

    def __str__(self):
        atoms = list(self.atoms)
        if len(atoms) == 1:
            atoms.append('')
        atoms = ', '.join(str(atom) for atom in atoms)

        bonders = list(self.bonders)
        if len(bonders) == 1:
            bonders.append('')
        bonders = ', '.join(str(atom) for atom in bonders)

        deleters = list(self.deleters)
        if len(deleters) == 1:
            deleters.append('')
        deleters = ', '.join(str(atom) for atom in deleters)

        return (
            f'FunctionalGroup(\n'
            f'    atoms=( {atoms} ), \n'
            f'    bonders=( {bonders} ), \n'
            f'    deleters=( {deleters} ), \n'
            f'    fg_type={self.fg_type.name}\n'
            ')'
        )


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

        bonds_made = reactor.finalize()

    An obvious question given this tutorial, is what reaction does
    :meth:`react` carry out? This is documented by :meth:`react`.
    However, :meth:`react` in most cases, will carry out a default
    reaction, which adds a bond between the bonder atoms of two
    functional groups. The bond order of the added bond is single by
    default but can be modified by editing :attr:`_bond_orders`. Here
    you will specify the :class:`ReactionKey` for a reaction and what
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
    as arguments. These are the functional groups which react. Within
    the method itself, the attribute :attr:`bonds_made` must be
    incremented by the number of bonds added by the reaction.

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
        mol : :class:`.ConstructedMolecule
            The molecule on which the reactor adds and removes atoms
            and bonds.

        """

        # The net number of bonds added.
        self._bonds_made = 0
        # The molecule from which atoms and bonds are added and
        # removed.
        self._mol = mol
        # The ids of atoms which are to be removed.
        self._deleter_ids = set()
        # The ids of bonds which are to be removed.
        self._deleter_bonds = set()

    def add_reaction(self, *fgs):
        """
        Creates bonds between functional groups.

        This function first looks at the functional groups provided via
        the `*fgs` argument and checks which functional groups are
        involved in the reaction. If the functional groups are handled
        by one of the custom reactions specified in
        :attr:`_custom_reactions` then that function is executed.

        In all other cases the function is assumed to have received two
        functional groups to react. In these functional groups, the
        bonder atoms have a bond added. The bond is single,
        unless otherwise specified in :attr:`_bond_orders`.

        Parameters
        ----------
        *fgs : :class:`FunctionalGroup`
            The functional groups to react.

        Returns
        -------
        None : :class:`NoneType`

        """

        names = (fg.fg_type.name for fg in fgs)
        reaction_key = _ReactionKey(*names)
        if reaction_key in self._custom_reactions:
            return self._custom_reactions[reaction_key](self, *fgs)

        for fg in fgs:
            self._deleter_ids.update(fg.get_deleter_ids())
            fg.atoms = tuple(
                a for a in fg.atoms if a.id not in self._deleter_ids
            )
            fg.deleters = ()

        bond = self._bond_orders.get(reaction_key, 1)
        fg1, fg2 = fgs
        self._mol.bonds.append(
            Bond(fg1.bonders[0], fg2.bonders[0], bond)
        )
        self._bonds_made += 1

    def add_periodic_reaction(self, direction, *fgs):
        """
        Like :func:`add_reaction` but adds :class:`.PeriodicBond`.

        Parameters
        ----------
        direction : :class:`list` of :class:`int`
            A 3 member :class:`list` describing the axes along which
            the created bonds are periodic. For example, ``[1, 0, 0]``
            means that he bonds are periodic along the x axis in the
            positive direction.

        *fgs : :class:`FunctionalGroup`
            The functional groups to react.

        Returns
        -------
        None : :class:`NoneType`

        """

        for fg in fgs:
            self._deleter_ids.update(fg.get_deleter_ids())
            fg.deleters = ()

        names = (fg.fg_type.name for fg in fgs)
        reaction_key = _ReactionKey(*names)
        if reaction_key in self._periodic_custom_reactions:
            rxn_fn = self._periodic_custom_reactions[reaction_key]
            return rxn_fn(self, direction, *fgs)

        bond = self._bond_orders.get(reaction_key, 1)

        # Make sure the direction of the periodic bond is maintained.
        fg1, fg2 = fgs
        self._mol.bonds.append(
            PeriodicBond(
                atom1=fg1.bonders[0],
                atom2=fg2.bonders[0],
                order=bond,
                direction=direction
            )
        )
        self._bonds_made += 1

    def finalize(self):
        """
        Finish performing reactions.

        Returns
        -------
        :class:`int`
            The number of bonds added.

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

        return self._bonds_made

    def _diol_with_dihalogen(self, fg1, fg2, dihalogen):
        """
        Creates bonds between functional groups.

        Parameters
        ----------
        fg1 : :class:`.FunctionalGroup`
            A functional group which undergoes the reaction.

        fg2 : :class:`.FunctionalGroup`
            A functional group which undergoes the reaction.

        halogen : :class:`str`
            The name of the dihalogen functional group to use. For
            example, ``'dibromine'`` or ``'difluorene'``.

        Returns
        -------
        None : :class:`NoneType`

        """

        diol = fg1 if fg1.fg_type.name == 'diol' else fg2
        dihalogen = fg2 if diol is fg1 else fg1

        distances = []
        for carbon in dihalogen.get_bonder_ids():
            for oxygen in diol.get_bonder_ids():
                d = self._mol.get_atom_distance(carbon, oxygen)
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
        self._mol.bonds.append(Bond(c1, o1, 1))
        self._mol.bonds.append(Bond(c2, o2, 1))
        self._bonds_made += 2

    def _boronic_acid_with_diol(self, fg1, fg2):
        """
        Creates bonds between functional groups.

        Parameters
        ----------
        fg1 : :class:`.FunctionalGroup`
            A functional group which undergoes the reaction.

        fg2 : :class:`.FunctionalGroup`
            A functional group which undergoes the reaction.

        Returns
        -------
        None : :class:`NoneType`

        """

        for fg in (fg1, fg2):
            self._deleter_ids.update(fg.get_deleter_ids())
            fg.atoms = tuple(
                a for a in fg.atoms if a.id not in self._deleter_ids
            )
            fg.deleters = ()

        boron = fg1 if fg1.fg_type.name == 'boronic_acid' else fg2
        diol = fg2 if boron is fg1 else fg1

        boron_atom = boron.bonders[0]
        self._mol.bonds.append(Bond(boron_atom, diol.bonders[0], 1))
        self._mol.bonds.append(Bond(boron_atom, diol.bonders[1], 1))
        self._bonds_made += 2

    def _ring_amine_with_ring_amine(self, fg1, fg2):
        """
        Creates bonds between functional groups.

        Parameters
        ----------
        fg1 : :class:`.FunctionalGroup`
            A functional group which undergoes the reaction.

        fg2 : :class:`.FunctionalGroup`
            A functional group which undergoes the reaction.

        Returns
        -------
        None : :class:`NoneType`

        """

        c1 = next(a for a in fg1.bonders if a.atomic_number == 6)
        n1 = next(a for a in fg1.bonders if a.atomic_number == 7)

        c2 = next(a for a in fg2.bonders if a.atomic_number == 6)
        n2 = next(a for a in fg2.bonders if a.atomic_number == 7)

        coords = self._mol.get_atom_coords(atom_ids=(c1, n1, c2, n2))
        n1_coord, n2_coord, c1_coord, c2_coord = coords

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

        self._mol.bonds.append(Bond(n1, n_joiner, 1))
        self._mol.bonds.append(Bond(n2, n_joiner, 1))
        self._mol.bonds.append(Bond(n_joiner, nh1, 1))
        self._mol.bonds.append(Bond(n_joiner, nh2, 1))

        self._mol.bonds.append(Bond(c1, nc_joiner1, 1))
        self._mol.bonds.append(Bond(n2, nc_joiner1, 1))
        self._mol.bonds.append(Bond(nc_joiner1, nc1h1, 1))
        self._mol.bonds.append(Bond(nc_joiner1, nc1h2, 1))

        self._mol.bonds.append(Bond(c2, nc_joiner2, 1))
        self._mol.bonds.append(Bond(n1, nc_joiner2, 1))
        self._mol.bonds.append(Bond(nc_joiner2, nc2h1, 1))
        self._mol.bonds.append(Bond(nc_joiner2, nc2h2, 1))

        self._bonds_made += 12

    # Maps a _ReactionKey for a given reaction to a custom
    # method which carries out the reaction. This means that
    # add_reaction will use that method for carrying out that
    # reaction instead.
    _custom_reactions = {

        _ReactionKey('boronic_acid', 'diol'):
            _boronic_acid_with_diol,

        _ReactionKey('diol', 'difluorene'):
            partial(_diol_with_dihalogen, dihalogen='difluorene'),

        _ReactionKey('diol', 'dibromine'):
            partial(_diol_with_dihalogen, dihalogen='dibromine'),

        _ReactionKey('ring_amine', 'ring_amine'):
            _ring_amine_with_ring_amine

    }

    # Maps a _ReactionKey for a given reaction to a custom
    # method which carries out the reaction. This means that
    # add_periodic_reaction() will use that method for carrying out
    # that reaction instead.
    _periodic_custom_reactions = {}


fg_types = {

    'amine': FGType(
        fg_smarts='[N]([H])[H]',
        bonder_smarts=['[$([N]([H])[H])]'],
        del_smarts=['[$([H][N][H])]']*2
    ),

    'aldehyde': FGType(
        fg_smarts='[C](=[O])[H]',
        bonder_smarts=['[$([C](=[O])[H])]'],
        del_smarts=['[$([O]=[C][H])]']
    ),

    'carboxylic_acid': FGType(
        fg_smarts='[C](=[O])[O][H]',
        bonder_smarts=['[$([C](=[O])[O][H])]'],
        del_smarts=[
            '[$([H][O][C](=[O]))]',
            '[$([O]([H])[C](=[O]))]'
        ]
    ),

    'amide': FGType(
        fg_smarts='[C](=[O])[N]([H])[H]',
        bonder_smarts=['[$([C](=[O])[N]([H])[H])]'],
        del_smarts=(
            ['[$([N]([H])([H])[C](=[O]))]'] +
            ['[$([H][N]([H])[C](=[O]))]']*2
        )
    ),

    'thioacid': FGType(
        fg_smarts='[C](=[O])[S][H]',
        bonder_smarts=['[$([C](=[O])[S][H])]'],
        del_smarts=['[$([H][S][C](=[O]))]', '[$([S]([H])[C](=[O]))]']
    ),

    'alcohol': FGType(
        fg_smarts='[O][H]',
        bonder_smarts=['[$([O][H])]'],
        del_smarts=['[$([H][O])]']
    ),

    'thiol': FGType(
        fg_smarts="[S][H]",
        bonder_smarts=['[$([S][H])]'],
        del_smarts=['[$([H][S])]']
    ),

    'bromine': FGType(
        fg_smarts='*[Br]',
        bonder_smarts=['[$(*[Br])]'],
        del_smarts=['[$([Br]*)]']
    ),

    'iodine': FGType(
        fg_smarts='*[I]',
        bonder_smarts=['[$(*[I])]'],
        del_smarts=['[$([I]*)]']
    ),

    'alkyne': FGType(
        fg_smarts='[C]#[C][H]',
        bonder_smarts=['[$([C]([H])#[C])]'],
        del_smarts=['[$([H][C]#[C])]']
    ),

    'terminal_alkene': FGType(
        fg_smarts='[C]=[C]([H])[H]',
        bonder_smarts=['[$([C]=[C]([H])[H])]'],
        del_smarts=(
            ['[$([H][C]([H])=[C])]']*2 +
            ['[$([C](=[C])([H])[H])]']
        )
    ),

    'boronic_acid': FGType(
        fg_smarts='[B]([O][H])[O][H]',
        bonder_smarts=['[$([B]([O][H])[O][H])]'],
        del_smarts=(
            ['[$([O]([H])[B][O][H])]']*2 +
            ['[$([H][O][B][O][H])]']*2
        )
    ),

    # This amine functional group only deletes one of the
    # hydrogen atoms when a bond is formed.
    'amine2': FGType(
        fg_smarts='[N]([H])[H]',
        bonder_smarts=['[$([N]([H])[H])]'],
        del_smarts=['[$([H][N][H])]']
    ),

    'secondary_amine': FGType(
        fg_smarts='[H][N]([#6])[#6]',
        bonder_smarts=[
            '[$([N]([H])([#6])[#6])]'
        ],
        del_smarts=['[$([H][N]([#6])[#6])]']
    ),

    'diol': FGType(
        fg_smarts='[H][O][#6]~[#6][O][H]',
        bonder_smarts=['[$([O]([H])[#6]~[#6][O][H])]']*2,
        del_smarts=['[$([H][O][#6]~[#6][O][H])]']*2
    ),

    'difluorene': FGType(
        fg_smarts='[F][#6]~[#6][F]',
        bonder_smarts=['[$([#6]([F])~[#6][F])]']*2,
        del_smarts=['[$([F][#6]~[#6][F])]']*2
    ),

    'dibromine': FGType(
        fg_smarts='[Br][#6]~[#6][Br]',
        bonder_smarts=['[$([#6]([Br])~[#6][Br])]']*2,
        del_smarts=['[$([Br][#6]~[#6][Br])]']*2
    ),

    'alkyne2': FGType(
        fg_smarts='[C]#[C][H]',
        bonder_smarts=['[$([C]#[C][H])]'],
        del_smarts=['[$([H][C]#[C])]', '[$([C](#[C])[H])]']
    ),

    'ring_amine': FGType(
        fg_smarts='[N]([H])([H])[#6]~[#6]([H])~[#6R1]',
        bonder_smarts=[
            '[$([N]([H])([H])[#6]~[#6]([H])~[#6R1])]',
            '[$([#6]([H])(~[#6R1])~[#6][N]([H])[H])]',
        ],
        del_smarts=(
                ['[$([H][N]([H])[#6]~[#6]([H])~[#6R1])]']*2 +
                ['[$([H][#6](~[#6R1])~[#6][N]([H])[H])]']
        )
    ),

}
