"""
Defines tools for dealing with functional groups and their reactions.

.. _`adding functional groups`:

Extending stk: Adding  more functional groups.
----------------------------------------------

If ``stk`` is to incorporate a new functional group, a new
:class:`FGInfo` instance should be added to
:data:`functional_groups`.

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
should be added to :data:`bond_orders`, along with the desired bond
order.

Supporting complex reactions.
.............................

During assembly, two functional groups are provided to
:func:`react`. By default, placing an :class:`FGInfo` instance into
:data:`functional_groups` will result in the creation of a single bond
between the atoms tagged as ``'bonder'`` in the two functional groups.
In addtion, any atoms tagged as ``'del'`` will be removed. The bond
order of the created bond can be modified by editing
:data:`bond_orders`.

However, some reactions cannot be described by a simple combination of
adding a bond while deleting some existing atoms. For example, consider
the aldol reaction:

    CH3C(=O)CH3 + CH3C(=O)CH3 --> CH3(=O)CH2C(OH)(CH3)CH3

Here a ketone is converted into an alcohol. In order to support more
complex conversions, a specific function needs to be defined which
modifies the molecule as desired. The function then needs
to be added to :data:`custom_reactions`. See
:func:`boronic_acid_with_diol`
as an example.

"""

from functools import partial
import numpy as np
from scipy.spatial.distance import euclidean
import rdkit.Chem.AllChem as rdkit
import rdkit.Geometry.rdGeometry as rdkit_geo
from collections import Counter
from ..utilities import AtomicPeriodicBond, flatten


class FGKey:
    """
    Used to create a key from a :class:`list` of fg names.

    Used by :data:`bond_orders`, :data:`custom_reactions` and
    :data:`periodic_custom_reactions`.

    Attributes
    ----------
    key : :class:`tuple`
        A unique key based on the functional groups provided to the
        intializer.

    """

    def __init__(self, fgs):
        """
        Intializer.

        Paramters
        ---------
        fgs : :class:`list` of :class:`str`
            A :class:`list` holding the names of functional groups.

        """
        c = Counter(fgs)
        self.key = tuple(sorted((key, value) for key, value in c.items()))

    def __eq__(self, other):
        return self.key == other.key

    def __hash__(self):
        return hash(self.key)

    def __repr__(self):
        fg_names = [name for name, count in self.key
                    for i in range(count)]
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

        Each string is SMARTS string which matches an atom in the
        functional group which is to be tagged as ``'bonder'``. The
        number represents how many matched atoms should be tagged, per
        functional group.

        In the example, ``Match(smarts='[$([N]([H])[H])]', n=1)``
        matches the nitrogen atom in the amine functional group. The
        ``n=1`` means that 1 nitrogen atom per functional group will be
        tagged as ``'bonder'``. The second
        ``Match(smarts='[$([H][N][H])]', n=1)``, matches the hydrogen
        atom in the amine functional group. Because ``n=1``, only 1 of
        hydrogen atom per amine functional group will be tagged
        ``'bonder'``. If instead
        ``Match(smarts='[$([H][N][H])]', n=2)`` was used, then both of
        the hydrogen atoms in the functional group would be tagged.

    del_smarts : :class:`list`
        Same as :attr:`bonder_smarts` but matched atoms are tagged
        as ``'del'``.

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

    """

    def __init__(self, id_, atom_ids, bonder_ids, deleter_ids, info):
        self.id = id_
        self.atom_ids = atom_ids
        self.bonder_ids = bonder_ids
        self.deleter_ids = deleter_ids
        self.info = info

    def __eq__(self, other):
        return (self.id == other.id and
                self.atom_ids == other.atom_ids and
                self.bonder_ids == other.bonder_ids and
                self.deleter_ids == other.deleter_ids and
                self.info == other.info)

    def __hash__(self):
        return id(self)

    def __repr__(self):
        return (f"FunctionalGroup(id={self.id!r}, "
                f"atom_ids={self.atom_ids!r}, "
                f"bonder_ids={self.bonder_ids!r}, "
                f"deleter_ids={self.deleter_ids!r}, "
                f"info={self.info!s})")

    def __str__(self):
        return repr(self)


class Reactor:
    """

    If some functional groups react via a special mechanism not covered
    in by the base "react()" function the function should be placed
    in this dict. The key should be a sorted tuple which holds the name
    of every functional group involved in the reaction along with how
    many such functional groups are invovled.

    """

    double = rdkit.rdchem.BondType.DOUBLE
    triple = rdkit.rdchem.BondType.TRIPLE
    bond_orders = {
        FGKey(['amine', 'aldehyde']): double,
        FGKey(['amide', 'aldehyde']): double,
        FGKey(['nitrile', 'aldehyde']): double,
        FGKey(['amide', 'amine']): double,
        FGKey(['terminal_alkene', 'terminal_alkene']): double,
        FGKey(['alkyne2', 'alkyne2']): triple
    }

    def __init__(self, mol):

        self.custom_reactions = {

            FGKey(['boronic_acid', 'diol']):
                self.boronic_acid_with_diol,

            FGKey(['diol', 'difluorene']):
                partial(self.diol_with_dihalogen,
                        dihalogen='difluorene'),

            FGKey(['diol', 'dibromine']):
                partial(self.diol_with_dihalogen,
                        dihalogen='dibromine'),

            FGKey(['phenyl_amine', 'phenyl_amine']):
                self.phenyl_amine_with_phenyl_amine

        }

        self.periodic_custom_reactions = {}

        self.mol = mol
        self.emol = rdkit.EditableMol(mol)
        self.periodic_bonds = []
        self.bonds_made = 0
        self.new_atom_coords = []
        self.deleters = []

    def react(self, *fgs):
        """
        Creates bonds between functional groups.

        This function first looks at the functional groups provided via
        the `*fgs` argument and checks which functional groups are
        involved in the reaction. If the functional groups are handled
        by one of the custom reactions specified in
        :data:`custom_reactions` then that function is executed.

        In all other cases the function is assumed to have received two
        functional groups to react via `*fgs`. In these functional
        groups the bonder atoms have a bond added. The bond is single,
        unless otherwise specified in :attr:`bond_orders`.

        Parameters
        ----------
        *fgs : :class:`FunctionalGroup`
            The functional groups to react.

        Returns
        -------
        None : :class:`NoneType`

        """

        self.deleters.extend(flatten(fg.deleter_ids for fg in fgs))

        names = [fg.info.name for fg in fgs]
        reaction_key = FGKey(names)
        if reaction_key in self.custom_reactions:
            return self.custom_reactions[reaction_key](*fgs)

        bond = self.bond_orders.get(reaction_key,
                                    rdkit.rdchem.BondType.SINGLE)
        fg1, fg2 = fgs
        self.emol.AddBond(fg1.bonder_ids[0], fg2.bonder_ids[0], bond)
        self.bonds_made += 1

    def periodic_react(self, direction, *fgs):
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

        Returns
        -------
        None : :class:`NoneType`

        """

        self.deleters.extend(flatten(fg.deleter_ids for fg in fgs))

        names = [fg.info.name for fg in fgs]
        reaction_key = FGKey(names)
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

        """

        if self.new_atom_coords:
            self.mol = self.emol.GetMol()
            conf = self.mol.GetConformer()
            for atom_id, coord in self.new_atom_coords:
                point3d = rdkit_geo.Point3D(*coord)
                conf.SetAtomPosition(atom_id, point3d)
            self.emol = rdkit.EditableMol(self.mol)

        if del_atoms:
            for atom_id in sorted(self.deleters, reverse=True):
                self.emol.RemoveAtom(atom_id)

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

    def phenyl_amine_with_phenyl_amine(self, fg1, fg2):
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

    FGInfo(name='phenyl_amine',
           fg_smarts='[N]([H])([H])[#6]~[#6]([H])~[#6]([H])',
           bonder_smarts=[
               Match(
                   smarts='[$([N]([H])([H])[#6]~[#6]([H])~[#6]([H]))]',
                   n=1
               ),
               Match(
                   smarts='[$([#6]([H])(~[#6]([H]))~[#6][N]([H])[H])]',
                   n=1
               ),
            ],
           del_smarts=[
               Match(
                   smarts='[$([H][N]([H])[#6]~[#6]([H])~[#6]([H]))]',
                   n=2
               ),
               Match(
                   smarts='[$([H][#6](~[#6]([H]))~[#6][N]([H])[H])]',
                   n=1
               )
           ])

    )

functional_group_infos = {fg.name: fg for fg in functional_groups}
