"""
Defines tools for dealing with functional groups and their reactions.

See the documentation of :class:`Reactor` to see how reactions between
functional groups are performed.

.. _`adding functional groups`:

Extending stk: Adding  more functional groups.
----------------------------------------------

If ``stk`` is to incorporate a new functional group, a new
:class:`FGInfo` instance should be added to
:data:`functional_groups`, which is defined in this module.

Adding a new :class:`FGInfo` instance to :data:`functional_groups` will
allow :meth:`.Topology.build` to connect the functional group to
all others during assembly. In most cases, nothing except adding this
instance should be necessary in order to incorporate new functional
groups.

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
from scipy.spatial.distance import euclidean
import rdkit.Chem.AllChem as rdkit
import rdkit.Geometry.rdGeometry as rdkit_geo
from collections import Counter
from ..utilities import AtomicPeriodicBond, flatten


class ReactionKey:
    """
    Used to create a key from a :class:`list` of fg names.

    In effect, this creates a unique id for each reaction, given a
    :class:`list` of functional groups involved in that reaction.

    Attributes
    ----------
    key : :class:`tuple`
        A unique key representing a reaction.

    """

    def __init__(self, *fg_names):
        """
        Intializer.

        Paramters
        ---------
        *fg_names : :class:`str`
            The names of functional groups involved in a reaction.

        """

        c = Counter(fg_names)
        self.key = tuple(
            sorted((key, value) for key, value in c.items())
        )

    def __eq__(self, other):
        return self.key == other.key

    def __hash__(self):
        return hash(self.key)

    def __repr__(self):
        fg_names = ', '.join(repr(name) for name, count in self.key
                             for i in range(count))

        return f'FGInfo({fg_names})'

    def __str__(self):
        return repr(self)


class Match:
    """
    A container for SMARTS queries.

    Attributes
    ----------
    smarts : :class:`str`
        A SMARTS string which matches some atoms.

    n : :class:`int`
        The maximum number of atoms to be matched by :attr:`smarts`,
        per functional group.

    """

    __slots__ = ['smarts', 'n']

    def __init__(self, smarts, n):
        self.smarts = smarts
        self.n = n

    def __eq__(self, other):
        return self.smarts == other.smarts and self.n == other.n

    def __repr__(self):
        return f'Match({self.smarts!r}, {self.n!r})'

    def __str__(self):
        return repr(self)


class FGInfo:
    """
    Contains key information about functional groups.

    The point of this class is to register which atoms of a functional
    group form bonds, and which are deleted during assembly of
    macromolecules.

    Attributes
    ----------
    name : :class:`str`
        The name of the functional group.

    fg_smarts : :class:`str`
        A SMARTS string which matches the functional group.

    bonder_smarts : :class:`list`
        A :class:`list` of the form

        .. code-block:: python

            bonder_smarts = [Match(smarts='[$([N]([H])[H])]', n=1),
                             Match(smarts='[$([H][N][H])]', n=1)]

        Each string is SMARTS string which matches a bonder atom in the
        functional group. The number represents how many matched atoms
        should be used as bonders, per functional group.

        In the example, ``Match(smarts='[$([N]([H])[H])]', n=1)``
        matches the nitrogen atom in the amine functional group. The
        ``n=1`` means that 1 nitrogen atom per functional group will be
        used as a bonder. The second
        ``Match(smarts='[$([H][N][H])]', n=1)``, matches the hydrogen
        atom in the amine functional group. Because ``n=1``, only 1 of
        hydrogen atoms per amine functional group will be used as a
        bonder. If instead ``Match(smarts='[$([H][N][H])]', n=2)`` was
        used, then both of the hydrogen atoms in the functional group
        would be used as bonders.

    del_smarts : :class:`list`
        Same as :attr:`bonder_smarts` but matched atoms are deleters

    """

    __slots__ = ['name', 'fg_smarts', 'bonder_smarts', 'del_smarts']

    def __init__(self, name, fg_smarts, bonder_smarts, del_smarts):
        """
        Initializes a :class:`FGInfo` instnace.

        Parameters
        ---------
        name : :class:`str`
            The name of the functional group.

        fg_smarts : :class:`str`
            A SMARTS string which matches the functional group.

        bonder_smarts : :class:`list`
            See :attr:`bonder_smarts`.

        del_smarts : :class:`list`
            See :attr:`del_smarts`.

        """

        self.name = name
        self.fg_smarts = fg_smarts
        self.bonder_smarts = bonder_smarts
        self.del_smarts = del_smarts

    def __eq__(self, other):
        return (self.name == other.name and
                self.fg_smarts == other.fg_smarts and
                self.bonder_smarts == other.bonder_smarts and
                self.del_smarts == other.del_smarts)

    def __repr__(self):
        return (f'FGInfo(name={self.name!r}, '
                f'fg_smarts={self.fg_smarts!r}, '
                f'bonder_smarts={self.bonder_smarts!r}, '
                f'del_smarts={self.del_smarts!r})')

    def __str__(self):
        return f'FGInfo({self.name!r})'


class FunctionalGroup:
    """
    Represents the functional group of a molecule.

    Attributes
    ----------
    id : :class:`int`
        The id of the functional group.

    atom_ids : :class:`tuple` of :class:`int`
        The ids of atoms in the functional group.

    bonder_ids : :class:`tuple` of :class:`int`
        The ids of bonder atoms in the functional group.

    deleter_ids : :class:`tuple` of :class:`int`
        The ids of deleter atoms in the functional group.

    info : :class:`FGInfo`
        The :class:`FGInfo` of the functional group type.

    """

    def __init__(self, id_, atom_ids, bonder_ids, deleter_ids, info):
        """
        Initialize a functional group.

        Parameters
        ----------
        id_ : :class:`int`
            The id of the functional group.

        atom_ids : :class:`tuple` of :class:`int`
            The ids of atoms in the functional group.

        bonder_ids : :class:`tuple` of :class:`int`
            The ids of bonder atoms in the functional group.

        deleter_ids : :class:`tuple` of :class:`int`
            The ids of deleter atoms in the functional group.

        info : :class:`FGInfo` or :class:`str`
            The :class:`FGInfo` of the functional group to which the
            functional group belongs. Can also be the name of the
            :class:`FGInfo`.

        """

        self.id = id_
        self.atom_ids = atom_ids
        self.bonder_ids = bonder_ids
        self.deleter_ids = deleter_ids

        if isinstance(info, str):
            self.info = next(fg_info for fg_info in functional_groups
                             if fg_info.name == info)
        else:
            self.info = info

    def shifted_fg(self, id_, shift):
        """
        Create a new :class:`FunctionalGroup` with shifted ids.

        Parameters
        ----------
        id_ : :class:`int`
            The id of the new functional group.

        shift : :class:`int`
            The number to shift the atom ids by.

        Returns
        -------
        :class:`FunctionalGroup`
            A :class:`FunctionalGroup` with all the atom ids atoms
            shifted upward by `shift`.

        """

        atom_ids = tuple(id + shift for id in self.atom_ids)
        bonder_ids = tuple(id + shift for id in self.bonder_ids)
        deleter_ids = tuple(id + shift for id in self.deleter_ids)

        return self.__class__(id_=id_,
                              atom_ids=atom_ids,
                              bonder_ids=bonder_ids,
                              deleter_ids=deleter_ids,
                              info=self.info)

    def remove_deleters(self, deleters):
        """
        Update and remove atom ids based on `deleters`.

        For each deleter atom that is smaller than an atom id, the
        atom id is decreased by 1. If the atom id is equal to a
        `deleters` atom id, it is removed.

        Parameters
        ----------
        deleters : :class:`list` of :class:`int`
            Ids of atoms which are being removed from the molecule.
            Must be sorted in ascending order.

        Returns
        -------
        None : :class:`NoneType`

        """

        self.atom_ids = self._remove_deleters(self.atom_ids,
                                              deleters)

        self.bonder_ids = self._remove_deleters(self.bonder_ids,
                                                deleters)

        self.deleter_ids = self._remove_deleters(self.deleter_ids,
                                                 deleters)

    def _remove_deleters(self, atom_ids, deleters):
        """
        Update and remove atom ids based on `deleters`.

        For each deleter atom that is smaller than an atom id, the
        atom id is decreased by 1. If the atom id is equal to a
        `deleters` atom id, it is removed.

        Parameters
        ----------
        atom_ids : :class:`tuple` of :class:`int`
            The atom ids which need to be updated.

        deleters : :class:`list` of :class:`int`
            The ids of atoms which are being removed. Must be sorted in
            ascending order.

        Returns
        -------
        :class:`tuple` of :class:`int`
            The atom ids in `atom_ids` after the update.

        """

        # Map each atom id to the number of smaller deleter ids.
        # If the atom id is to be deleted, map to None.
        id_changes = {atom_id: 0 for atom_id in atom_ids}

        for atom_id in atom_ids:
            for deleter in deleters:
                if atom_id > deleter:
                    id_changes[atom_id] += 1
                elif atom_id == deleter:
                    id_changes[atom_id] = None
                # Because the deleters are sorted, there will be no
                # more deleters which are smaller than atom_id.
                if atom_id < deleter:
                    break

        # Apply the id changes to the original atom ids.
        new_atom_ids = []
        for atom_id in atom_ids:
            change = id_changes[atom_id]
            if change is None:
                continue
            else:
                new_atom_ids.append(atom_id - change)

        return tuple(new_atom_ids)

    def __eq__(self, other):
        return (self.id == other.id and
                self.atom_ids == other.atom_ids and
                self.bonder_ids == other.bonder_ids and
                self.deleter_ids == other.deleter_ids and
                self.info == other.info)

    def __hash__(self):
        return id(self)

    def __repr__(self):
        return (f"FunctionalGroup(id_={self.id!r}, "
                f"atom_ids={self.atom_ids!r}, "
                f"bonder_ids={self.bonder_ids!r}, "
                f"deleter_ids={self.deleter_ids!r}, "
                f"info={self.info.name!r})")

    def __str__(self):
        return repr(self)


class Reactor:
    """
    Performs reactions between functional groups of a molecule.

    This class is responsible for reacting the functional groups in a
    molecule during assembly. First, an instance of this class is
    initialized with a :class:`rdkit.Mol`. This is the molecule which
    is going to have atom and bonds added and removed.

    .. code-block:: python

        mol = rdkit.MolFromMolFile(...)
        reactor = Reactor(mol)

    We force the molecule to have atoms and bonds add or removed
    between certain functional groups by using the :meth:`react`
    method.

    .. code-block:: python

        # Represents a functional group found in mol.
        fg1 = FunctionalGroup(id_=0,
                              atom_ids=[1, 34, 3],
                              bonder_ids=[1],
                              deleter_ids=[34, 3],
                              info=FGInfo('amine')
        )

        # Represents another functional group found in mol.
        fg2 = FunctionalGroup(id_=1,
                              atom_ids=[10, 2, 12],
                              bonder_ids=[12],
                              deleter_ids=[2, 10],
                              info=FGInfo('aldehyde')
        )

        # Carry out a reaction between the atoms in fg1 and fg2.
        reactor.react(fg1, fg2)

        # Lets assume there a further functional groups in mol that
        # we wish to react.
        fg3 = FunctionalGroup(...)
        fg4 = FunctionalGroup(...)
        fg5 = FunctionalGroup(...)
        fg6 = FunctionalGroup(...)
        fg7 = FunctionalGroup(...)

        reactor.react(fg3, fg4)

        # Some reactions can take multiple functional groups.
        # You can put in as many functional groups as you like, given
        # an appropriate reaction is defined.
        reactor.react(fg5, fg6, fg7)

    Once we are done carrying out reactions on the molecule we can
    get the resulting molecule.

    .. code-block:: python

        # product is an rdkit molecule, with the earlier reactions
        # carried out and all deleter atoms removed.
        product = reactor.result(del_atoms=True)

    Note that after :meth:`result` is run, all :class:`FunctionalGroup`
    instances which were passed to :meth:`react` have their atom ids
    updated to account for any atoms which have been deleted. If
    there are functional groups which were not reacted, but have
    atom ids which should be updated, they should be added to
    :attr:`func_groups`.

    An obvious question given this tutorial, is what reaction does
    :meth:`react` carry out? This is documented by :meth:`react`.
    However, react in most cases, will carry out a default reaction,
    which adds a bond between the bonder atoms of two functional
    groups. The bond order of the added bond is single by default but
    can be modified by editing :attr:`bond_orders`. Here you will
    specify the :class:`ReactionKey` for a reaction and what bond order
    you want that reaction to use.

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
    incremented by the number of bonds added by the reaction. Beyond
    that, the method should operate :attr:`emol` and add and remove
    bonds from it as necessary. If the method adds atoms to
    :attr:`emol` it should update :attr:`new_atom_coords`.

    Once the method is defined, :attr:`custom_reactions` needs to
    be updated.

    Attributes
    ----------
    bond_orders : :class:`dict`
        When the default reaction is performed by :meth:`react`,
        if the bond added between the two functional groups is not
        single, the desired bond order should be placed in this
        dictionary. The dictionary maps the reaction's
        :class:`ReactionKey` to the desired bond order.

    custom_reactions : :class:`dict`
        Maps a :class:`ReactionKey` for a given reaction to a custom
        method which carries out the reaction. This means that
        :meth:`react` will use that method for carrying out that
        reaction instead.

    periodic_custom_reactions : :class:`dict`
        Maps a :class:`ReactionKey` for a given reaction to a custom
        method which carries out the reaction. This means that
        :meth:`periodic_react` will use that method for carrying out
        that reaction instead.

    mol : :class:`rdkit.Mol`
        The molecule on which the reactor adds and removes atoms and
        bonds.

    emol : :class:`rdkit.EditableMol`
        An editable version of :attr:`mol`. Used for adding and
        removing atoms and bonds.

    periodic_bonds : :class:`list` of :class:`.AtomicPeriodicBond`
        The periodic bonds added by the reactor.

    bonds_made : :class:`int`
        The number of bonds added.

    new_atom_coords : :class:`list` of :class:`tuple`
        When a new atom is added by the reactor, its desired
        desired coordinates are placed here. The :class:`list`
        has the form

        .. code-block:: python

            new_atom_coords = [
                (32, np.array([12.1, 3.5, 0.1])),
                (41, np.array([13.1, -42.1, 2.]))
            ]

        where the first element of each tuple is the atom id of a
        newly added atom and the second element is a
        :class:`numpy.ndarray` holding its desired coordinate.

    deleters : :class:`list` of :class:`int`
        The ids of atoms which are to be removed.

    func_groups : :class:`list` of :class:`FunctionalGroup`
        The functional groups which the reactor keeps up to date.

    """

    double = rdkit.rdchem.BondType.DOUBLE
    triple = rdkit.rdchem.BondType.TRIPLE
    bond_orders = {
        ReactionKey('amine', 'aldehyde'): double,
        ReactionKey('amide', 'aldehyde'): double,
        ReactionKey('nitrile', 'aldehyde'): double,
        ReactionKey('amide', 'amine'): double,
        ReactionKey('terminal_alkene', 'terminal_alkene'): double,
        ReactionKey('alkyne2', 'alkyne2'): triple
    }

    def __init__(self, mol=None):
        """
        Initialize a :class:`Reactor`.

        Parameters
        ----------
        mol : :class:`rdkit.Mol`, optional
            The molecule on which the reactor adds and removes atoms
            and bonds.

        """

        self.custom_reactions = {

            ReactionKey('boronic_acid', 'diol'):
                self.boronic_acid_with_diol,

            ReactionKey('diol', 'difluorene'):
                partial(self.diol_with_dihalogen,
                        dihalogen='difluorene'),

            ReactionKey('diol', 'dibromine'):
                partial(self.diol_with_dihalogen,
                        dihalogen='dibromine'),

            ReactionKey('ring_amine', 'ring_amine'):
                self.ring_amine_with_ring_amine

        }

        self.periodic_custom_reactions = {}

        if mol is not None:
            self.set_molecule(mol)

        self.periodic_bonds = []
        self.bonds_made = 0
        self.new_atom_coords = []
        self.deleters = []
        self.func_groups = []

    def set_molecule(self, mol):
        """
        Update :attr:`mol` and :attr:`emol`.

        Parameters
        ----------
        mol : :class:`rdkit.Mol`, optional
            The molecule on which the reactor adds and removes atoms
            and bonds.

        Returns
        -------
        None : :class:`NoneType`

        """

        self.mol = mol
        self.emol = rdkit.EditableMol(mol)

    def react(self, *fgs, track_fgs=True):
        """
        Creates bonds between functional groups.

        This function first looks at the functional groups provided via
        the `*fgs` argument and checks which functional groups are
        involved in the reaction. If the functional groups are handled
        by one of the custom reactions specified in
        :attr:`custom_reactions` then that function is executed.

        In all other cases the function is assumed to have received two
        functional groups to react. In these functional groups, the
        bonder atoms have a bond added. The bond is single,
        unless otherwise specified in :attr:`bond_orders`.

        Parameters
        ----------
        *fgs : :class:`FunctionalGroup`
            The functional groups to react.

        track_fgs : :class:`bool`, optional
            Toggles if functional groups are added to
            :attr:`func_groups`.

        Returns
        -------
        None : :class:`NoneType`

        """

        self.deleters.extend(flatten(fg.deleter_ids for fg in fgs))
        if track_fgs:
            self.func_groups.extend(fgs)

        names = (fg.info.name for fg in fgs)
        reaction_key = ReactionKey(*names)
        if reaction_key in self.custom_reactions:
            return self.custom_reactions[reaction_key](*fgs)

        bond = self.bond_orders.get(reaction_key,
                                    rdkit.rdchem.BondType.SINGLE)
        fg1, fg2 = fgs
        self.emol.AddBond(fg1.bonder_ids[0], fg2.bonder_ids[0], bond)
        self.bonds_made += 1

    def periodic_react(self, direction, *fgs, track_fgs=True):
        """
        Like :func:`react` but returns periodic bonds.

        As periodic bonds are returned, no bonds are added to
        :attr:`emol`.

        Parameters
        ----------
        direction : :class:`list` of :class:`int`
            A 3 member list describing the axes along which the created
            bonds are periodic. For example, ``[1, 0, 0]`` means that
            the bonds are periodic along the x axis in the positive
            direction.

        *fgs : :class:`FunctionalGroup`
            The functional groups to react.

        track_fgs : :class:`bool`, optional
            Toggles if functional groups are added to
            :attr:`func_groups`.

        Returns
        -------
        None : :class:`NoneType`

        """

        self.deleters.extend(flatten(fg.deleter_ids for fg in fgs))
        if track_fgs:
            self.func_groups.extend(fgs)

        names = (fg.info.name for fg in fgs)
        reaction_key = ReactionKey(*names)
        if reaction_key in self.periodic_custom_reactions:
            rxn_fn = self.periodic_custom_reactions[reaction_key]
            return rxn_fn(direction, *fgs)

        bond = self.bond_orders.get(reaction_key,
                                    rdkit.rdchem.BondType.SINGLE)

        # Make sure the direction of the periodic bond is maintained.
        fg1, fg2 = fgs
        self.periodic_bonds.append(
            AtomicPeriodicBond(fg1.bonder_ids[0],
                               fg2.bonder_ids[0],
                               bond,
                               direction)
        )
        self.bonds_made += 1

    def result(self, del_atoms):
        """
        Creates the molecule after all reactions have been done.

        This method will also update all the functional groups
        passed with :meth:`react` calls. It will update the atom
        ids to account for the fact that atoms have been deleted.

        Parameters
        ----------
        del_atoms : :class:`bool`
            Toggles if deleter atoms should be removed from the
            product molecule.

        Returns
        -------
        :class:`rdkit.Mol`
            The product molecule.

        """

        # If new atoms were added, update the positions in the
        # conformer.
        if self.new_atom_coords:
            self.mol = self.emol.GetMol()
            conf = self.mol.GetConformer()
            for atom_id, coord in self.new_atom_coords:
                point3d = rdkit_geo.Point3D(*coord)
                conf.SetAtomPosition(atom_id, point3d)
            self.emol = rdkit.EditableMol(self.mol)

        if del_atoms:
            # Needs to be sorted for fg.remove_deleters.
            deleters = sorted(self.deleters)

            # Go in reverse order else atom ids change during loop.
            for atom_id in reversed(deleters):
                self.emol.RemoveAtom(atom_id)

            # Update all the functional groups to account for the fact
            # that atoms have been removed.
            for fg in self.func_groups:
                fg.remove_deleters(deleters)

        return self.emol.GetMol()

    def diol_with_dihalogen(self, fg1, fg2, dihalogen):
        """
        Creates bonds between functional groups.

        Parameters
        ----------
        fg1 : :class:`.FunctionalGroup`
            A functional group which undergoes the reaction.

        fg2 : :class:`.FunctionalGroup`
            A functional group which undergoes the reaction.

        halogen : :class:`str`
            The name of the dihalogen functional group to use. For example,
            ``'dibromine'`` or ``'difluorene'``.

        Returns
        -------
        None : :class:`NoneType`

        """

        conf = self.mol.GetConformer()
        bond = rdkit.rdchem.BondType.SINGLE

        diol = fg1 if fg1.info.name == 'diol' else fg2
        dihalogen = fg2 if diol is fg1 else fg1

        distances = []
        for carbon in dihalogen.bonder_ids:
            c_coord = np.array([*conf.GetAtomPosition(carbon)])
            for oxygen in diol.bonder_ids:
                o_coord = np.array([*conf.GetAtomPosition(oxygen)])
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
        self.emol.AddBond(c1, o1, bond)
        self.emol.AddBond(c2, o2, bond)
        self.bonds_made += 2

    def boronic_acid_with_diol(self, fg1, fg2):
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

        bond = rdkit.rdchem.BondType.SINGLE

        boron = fg1 if fg1.info.name == 'boronic_acid' else fg2
        diol = fg2 if boron is fg1 else fg1

        boron_atom = boron.bonder_ids[0]
        self.emol.AddBond(boron_atom, diol.bonder_ids[0], bond)
        self.emol.AddBond(boron_atom, diol.bonder_ids[1], bond)
        self.bonds_made += 2

    def ring_amine_with_ring_amine(self, fg1, fg2):
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

        c1 = next(a for a in fg1.bonder_ids
                  if self.mol.GetAtomWithIdx(a).GetAtomicNum() == 6)
        n1 = next(a for a in fg1.bonder_ids
                  if self.mol.GetAtomWithIdx(a).GetAtomicNum() == 7)

        c2 = next(a for a in fg2.bonder_ids
                  if self.mol.GetAtomWithIdx(a).GetAtomicNum() == 6)
        n2 = next(a for a in fg2.bonder_ids
                  if self.mol.GetAtomWithIdx(a).GetAtomicNum() == 7)

        conf = self.mol.GetConformer()
        n1_pos = np.array([*conf.GetAtomPosition(n1)])
        n2_pos = np.array([*conf.GetAtomPosition(n2)])

        c1_pos = np.array([*conf.GetAtomPosition(c1)])
        c2_pos = np.array([*conf.GetAtomPosition(c2)])

        n_joiner = self.emol.AddAtom(rdkit.Atom(6))
        n_joiner_pos = (n1_pos + n2_pos) / 2
        self.new_atom_coords.append((n_joiner, n_joiner_pos))

        nh1 = self.emol.AddAtom(rdkit.Atom(1))
        nh1_pos = n_joiner_pos + np.array([0, 0, 1])
        self.new_atom_coords.append((nh1, nh1_pos))

        nh2 = self.emol.AddAtom(rdkit.Atom(1))
        nh2_pos = n_joiner_pos + np.array([0, 0, -1])
        self.new_atom_coords.append((nh2, nh2_pos))

        nc_joiner1 = self.emol.AddAtom(rdkit.Atom(6))
        nc_joiner1_pos = (c1_pos + n2_pos) / 2
        self.new_atom_coords.append((nc_joiner1, nc_joiner1_pos))

        nc1h1 = self.emol.AddAtom(rdkit.Atom(1))
        nc1h1_pos = nc_joiner1_pos + np.array([0, 0, 1])
        self.new_atom_coords.append((nc1h1, nc1h1_pos))

        nc1h2 = self.emol.AddAtom(rdkit.Atom(1))
        nc1h2_pos = nc_joiner1_pos + np.array([0, 0, -1])
        self.new_atom_coords.append((nc1h2, nc1h2_pos))

        nc_joiner2 = self.emol.AddAtom(rdkit.Atom(6))
        nc_joiner2_pos = (c2_pos + n1_pos) / 2
        self.new_atom_coords.append((nc_joiner2, nc_joiner2_pos))

        nc2h1 = self.emol.AddAtom(rdkit.Atom(1))
        nc2h1_pos = nc_joiner2_pos + np.array([0, 0, 1])
        self.new_atom_coords.append((nc2h1, nc2h1_pos))

        nc2h2 = self.emol.AddAtom(rdkit.Atom(1))
        nc2h2_pos = nc_joiner2_pos + np.array([0, 0, -1])
        self.new_atom_coords.append((nc2h2, nc2h2_pos))

        single = rdkit.rdchem.BondType.SINGLE
        self.emol.AddBond(n1, n_joiner, single)
        self.emol.AddBond(n2, n_joiner, single)
        self.emol.AddBond(n_joiner, nh1, single)
        self.emol.AddBond(n_joiner, nh2, single)

        self.emol.AddBond(c1, nc_joiner1, single)
        self.emol.AddBond(n2, nc_joiner1, single)
        self.emol.AddBond(nc_joiner1, nc1h1, single)
        self.emol.AddBond(nc_joiner1, nc1h2, single)

        self.emol.AddBond(c2, nc_joiner2, single)
        self.emol.AddBond(n1, nc_joiner2, single)
        self.emol.AddBond(nc_joiner2, nc2h1, single)
        self.emol.AddBond(nc_joiner2, nc2h2, single)

        self.bonds_made += 6


functional_groups = (

    FGInfo(name="amine",
           fg_smarts="[N]([H])[H]",
           bonder_smarts=[Match(smarts="[$([N]([H])[H])]", n=1)],
           del_smarts=[Match(smarts="[$([H][N][H])]", n=2)]),

    FGInfo(name="aldehyde",
           fg_smarts="[C](=[O])[H]",
           bonder_smarts=[Match(smarts="[$([C](=[O])[H])]", n=1)],
           del_smarts=[Match(smarts="[$([O]=[C][H])]", n=1)]),

    FGInfo(name="carboxylic_acid",
           fg_smarts="[C](=[O])[O][H]",
           bonder_smarts=[Match(smarts="[$([C](=[O])[O][H])]", n=1)],
           del_smarts=[Match(smarts="[$([H][O][C](=[O]))]", n=1),
                       Match(smarts="[$([O]([H])[C](=[O]))]", n=1)]),


    FGInfo(name="amide",
           fg_smarts="[C](=[O])[N]([H])[H]",
           bonder_smarts=[
                Match(smarts="[$([C](=[O])[N]([H])[H])]", n=1)
           ],
           del_smarts=[
                Match(smarts="[$([N]([H])([H])[C](=[O]))]", n=1),
                Match(smarts="[$([H][N]([H])[C](=[O]))]", n=2)
           ]),

    FGInfo(name="thioacid",
           fg_smarts="[C](=[O])[S][H]",
           bonder_smarts=[Match(smarts="[$([C](=[O])[S][H])]", n=1)],
           del_smarts=[Match(smarts="[$([H][S][C](=[O]))]", n=1),
                       Match(smarts="[$([S]([H])[C](=[O]))]", n=1)]),

    FGInfo(name="alcohol",
           fg_smarts="[O][H]",
           bonder_smarts=[Match(smarts="[$([O][H])]", n=1)],
           del_smarts=[Match(smarts="[$([H][O])]", n=1)]),

    FGInfo(name="thiol",
           fg_smarts="[S][H]",
           bonder_smarts=[Match(smarts="[$([S][H])]", n=1)],
           del_smarts=[Match(smarts="[$([H][S])]", n=1)]),

    FGInfo(name="bromine",
           fg_smarts="*[Br]",
           bonder_smarts=[Match(smarts="[$(*[Br])]", n=1)],
           del_smarts=[Match(smarts="[$([Br]*)]", n=1)]),

    FGInfo(name="iodine",
           fg_smarts="*[I]",
           bonder_smarts=[Match(smarts="[$(*[I])]", n=1)],
           del_smarts=[Match(smarts="[$([I]*)]", n=1)]),

    FGInfo(name='alkyne',
           fg_smarts='[C]#[C][H]',
           bonder_smarts=[Match(smarts='[$([C]([H])#[C])]', n=1)],
           del_smarts=[Match(smarts='[$([H][C]#[C])]', n=1)]),

    FGInfo(name='terminal_alkene',
           fg_smarts='[C]=[C]([H])[H]',
           bonder_smarts=[Match(smarts='[$([C]=[C]([H])[H])]', n=1)],
           del_smarts=[Match(smarts='[$([H][C]([H])=[C])]', n=2),
                       Match(smarts='[$([C](=[C])([H])[H])]', n=1)]),

    FGInfo(name='boronic_acid',
           fg_smarts='[B]([O][H])[O][H]',
           bonder_smarts=[Match(smarts='[$([B]([O][H])[O][H])]', n=1)],
           del_smarts=[Match(smarts='[$([O]([H])[B][O][H])]', n=2),
                       Match(smarts='[$([H][O][B][O][H])]', n=2)]),

    # This amine functional group only deletes one of the
    # hydrogen atoms when a bond is formed.
    FGInfo(name="amine2",
           fg_smarts="[N]([H])[H]",
           bonder_smarts=[Match(smarts="[$([N]([H])[H])]", n=1)],
           del_smarts=[Match(smarts="[$([H][N][H])]", n=1)]),

    FGInfo(name="secondary_amine",
           fg_smarts="[H][N]([#6])[#6]",
           bonder_smarts=[
                Match(smarts="[$([N]([H])([#6])[#6])]", n=1)
           ],
           del_smarts=[Match(smarts="[$([H][N]([#6])[#6])]", n=1)]),

    FGInfo(name='diol',
           fg_smarts='[H][O][#6]~[#6][O][H]',
           bonder_smarts=[
               Match(smarts='[$([O]([H])[#6]~[#6][O][H])]', n=2)
           ],
           del_smarts=[
               Match(smarts='[$([H][O][#6]~[#6][O][H])]', n=2)
           ]),

    FGInfo(name='difluorene',
           fg_smarts='[F][#6]~[#6][F]',
           bonder_smarts=[Match(smarts='[$([#6]([F])~[#6][F])]', n=2)],
           del_smarts=[Match(smarts='[$([F][#6]~[#6][F])]', n=2)]),

    FGInfo(name='dibromine',
           fg_smarts='[Br][#6]~[#6][Br]',
           bonder_smarts=[
               Match(smarts='[$([#6]([Br])~[#6][Br])]', n=2)
           ],
           del_smarts=[Match(smarts='[$([Br][#6]~[#6][Br])]', n=2)]),

    FGInfo(name='alkyne2',
           fg_smarts='[C]#[C][H]',
           bonder_smarts=[Match(smarts='[$([C]#[C][H])]', n=1)],
           del_smarts=[Match(smarts='[$([H][C]#[C])]', n=1),
                       Match(smarts='[$([C](#[C])[H])]', n=1)]),

    FGInfo(name='ring_amine',
           fg_smarts='[N]([H])([H])[#6]~[#6]([H])~[#6R1]',
           bonder_smarts=[
               Match(
                   smarts='[$([N]([H])([H])[#6]~[#6]([H])~[#6R1])]',
                   n=1
               ),
               Match(
                   smarts='[$([#6]([H])(~[#6R1])~[#6][N]([H])[H])]',
                   n=1
               ),
            ],
           del_smarts=[
               Match(
                   smarts='[$([H][N]([H])[#6]~[#6]([H])~[#6R1])]',
                   n=2
               ),
               Match(
                   smarts='[$([H][#6](~[#6R1])~[#6][N]([H])[H])]',
                   n=1
               )
           ])

    )

functional_group_infos = {fg.name: fg for fg in functional_groups}
