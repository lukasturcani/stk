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

import numpy as np
from scipy.spatial.distance import euclidean
import rdkit.Chem.AllChem as rdkit
import rdkit.Geometry.rdGeometry as rdkit_geo
from collections import Counter
from ..utilities import AtomicPeriodicBond


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


def fg_name(mol, fg):
    """
    Retruns the name of the functional group with id `fg`.

    Parameters
    ----------
    mol : :class:`rdkit.Chem.rdchem.Mol`
        An ``rdkit`` molecule with its functional groups tagged.

    fg : :class:`int`
        The id of a functional group as given by the 'fg_id' property.

    Returns
    -------
    :class:`str`
        The name of a functional group.

    """

    for atom in mol.GetAtoms():
        if atom.HasProp('fg_id') and atom.GetIntProp('fg_id') == fg:
            return atom.GetProp('fg')
    raise RuntimeError(f'No functional group with id {fg} found.')


def react(mol, del_atoms, *fgs):
    """
    Crates bonds between functional groups.

    This function first looks at the functional group ids provided via
    the `*fgs` argument and checks which functional groups are
    involved in the reaction. If the functional groups are handled
    by one of the custom reactions specified in
    :data:`custom_reactions` then that function is executed.

    In all other cases the function is assumed to have received two
    functional groups to react via `*fgs`. In these functional groups
    the atoms tagged ``'del'`` are deleted and the atoms tagged
    ``'bonder'`` have a bond added. The bond is a single, unless
    specified otherwise in :data:`bond_orders`.

    Parameters
    ----------
    mol : :class:`rdkit.Chem.rdchem.Mol`
        A molecule being assembled.

    del : :class:`bool`
        Toggles if atoms with the ``'del'`` property are deleted.

    *fgs : :class:`int`
        The ids of the functional groups to react. The ids are held
        by atom of `mol` in the ``'fg_id'`` property.

    Returns
    -------
    :class:`tuple`
        The first element is an :class:`rdkit.Chem.rdchem.Mol`. It is
        the molecule with bonds added between the functional groups.

        The second element is a :class:`int`. It is the number
        of bonds added.

    """

    names = [fg_name(mol, fg) for fg in fgs]
    reaction_key = FGKey(names)
    if reaction_key in custom_reactions:
        return custom_reactions[reaction_key](mol, del_atoms, *fgs)

    emol = rdkit.EditableMol(mol)

    bonders = []
    for atom in mol.GetAtoms():
        if not (atom.HasProp('fg_id') and atom.GetIntProp('fg_id') in fgs):
            continue
        if atom.HasProp('bonder'):
            bonders.append(atom.GetIdx())

    bond = bond_orders.get(reaction_key, rdkit.rdchem.BondType.SINGLE)
    bonder1, bonder2 = bonders
    emol.AddBond(bonder1, bonder2, bond)

    for atom in reversed(mol.GetAtoms()):
        if not (atom.HasProp('fg_id') and atom.GetIntProp('fg_id') in fgs):
            continue

        if atom.HasProp('del') and del_atoms:
            emol.RemoveAtom(atom.GetIdx())

    return emol.GetMol(), 1


def periodic_react(mol, del_atoms, direction, *fgs):
    """
    Like :func:`react` but returns periodic bonds.

    As periodic bonds are returned, no bonds are added to `mol`.

    Parameters
    ----------
    mol : :class:`rdkit.Chem.rdchem.Mol`
        A molecule being assembled.

    del : :class:`bool`
        Toggles if atoms with the ``'del'`` property are deleted.

    direction : :class:`list` of :class:`int`
        A 3 member list describing the axes along which the created
        bonds are periodic. For example, ``[1, 0, 0]`` means that the
        bonds are periodic along the x axis in the positive direction.

    *fgs : :class:`int`
        The ids of the functional groups to react. The ids are held
        by atom of `mol` in the ``'fg_id'`` property.

    Returns
    -------
    :class:`tuple`
        The first element is an :class:`rdkit.Chem.rdchem.Mol`. It is
        the molecule after the reaction.

        The second element is a :class:`int`. It is the number
        of bonds added.

        The third element is a :class:`list` holding
        :class:`.AtomicPeriodicBond`.

    """

    names = [fg_name(mol, fg) for fg in fgs]
    reaction_key = FGKey(names)
    if reaction_key in periodic_custom_reactions:
        return periodic_custom_reactions[reaction_key](mol,
                                                       del_atoms,
                                                       direction,
                                                       *fgs)

    emol = rdkit.EditableMol(mol)

    bonders = {}
    for atom in mol.GetAtoms():
        if not (atom.HasProp('fg_id') and atom.GetIntProp('fg_id') in fgs):
            continue
        if atom.HasProp('bonder'):
            bonders[atom.GetIntProp('fg_id')] = atom.GetIntProp('bonder')

    bond = bond_orders.get(FGKey(names), rdkit.rdchem.BondType.SINGLE)

    # Make sure the direction of the periodic bond is maintained.
    fg1, fg2 = fgs
    bonder1, bonder2 = bonders[fg1], bonders[fg2]
    periodic_bonds = [AtomicPeriodicBond(bonder1,
                                         bonder2,
                                         bond,
                                         direction)]

    for atom in reversed(mol.GetAtoms()):
        if not (atom.HasProp('fg_id') and atom.GetIntProp('fg_id') in fgs):
            continue

        if atom.HasProp('del') and del_atoms:
            emol.RemoveAtom(atom.GetIdx())

    return emol.GetMol(), 1, periodic_bonds


def diol_with_difluorne(mol, del_atoms, fg1, fg2):
    """
    Crates bonds between functional groups.

    Parameters
    ----------
    mol : :class:`rdkit.Chem.rdchem.Mol`
        A molecule being assembled.

    del : :class:`bool`
        Toggles if atoms with the ``'del'`` property are deleted.

    fg1 : :class:`int`
        The id of the first functional group which
        is to be joined, as given by the 'fg_id' property.

    fg2 : :class:`int`
        The id of the second functional group which
        is to be joined, as given by the 'fg_id' property.

    Returns
    -------
    :class:`tuple`
        The first element is an :class:`rdkit.Chem.rdchem.Mol`. It is
        the molecule with bonds added between the functional groups.

        The second element is a :class:`int`. It is the number
        of bonds added.

    """

    bond = rdkit.rdchem.BondType.SINGLE
    fgs = {fg1, fg2}
    oxygens = []
    carbons = []
    deleters = []

    for a in reversed(mol.GetAtoms()):
        if not a.HasProp('fg_id') or a.GetIntProp('fg_id') not in fgs:
            continue

        if a.HasProp('del'):
            deleters.append(a)

        if a.GetProp('fg') == 'difluorene' and a.HasProp('bonder'):
            carbons.append(a)

        if a.GetProp('fg') == 'diol' and a.HasProp('bonder'):
            oxygens.append(a)

    conf = mol.GetConformer()
    distances = []
    for c in carbons:
        cpos = np.array([*conf.GetAtomPosition(c.GetIdx())])
        for o in oxygens:
            opos = np.array([*conf.GetAtomPosition(o.GetIdx())])
            d = euclidean(cpos, opos)
            distances.append((d, c.GetIdx(), o.GetIdx()))
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
    emol = rdkit.EditableMol(mol)
    emol.AddBond(c1, o1, bond)
    emol.AddBond(c2, o2, bond)

    if del_atoms:
        for a in deleters:
            emol.RemoveAtom(a.GetIdx())

    return emol.GetMol(), 2


def boronic_acid_with_diol(mol, del_atoms, fg1, fg2):
    """
    Crates bonds between functional groups.

    Parameters
    ----------
    mol : :class:`rdkit.Chem.rdchem.Mol`
        A molecule being assembled.

    del : :class:`bool`
        Toggles if atoms with the ``'del'`` property are deleted.

    fg1 : :class:`int`
        The id of the first functional group which
        is to be joined, as given by the 'fg_id' property.

    fg2 : :class:`int`
        The id of the second functional group which
        is to be joined, as given by the 'fg_id' property.

    Returns
    -------
    :class:`tuple`
        The first element is an :class:`rdkit.Chem.rdchem.Mol`. It is
        the molecule with bonds added between the functional groups.

        The second element is a :class:`int`. It is the number
        of bonds added.

    """

    bond = rdkit.rdchem.BondType.SINGLE
    fgs = {fg1, fg2}
    oxygens = []
    deleters = []

    for a in reversed(mol.GetAtoms()):
        if not a.HasProp('fg_id') or a.GetIntProp('fg_id') not in fgs:
            continue

        if a.HasProp('del'):
            deleters.append(a)

        if a.GetProp('fg') == 'boronic_acid' and a.HasProp('bonder'):
            boron = a

        if a.GetProp('fg') == 'diol' and a.HasProp('bonder'):
            oxygens.append(a)

    emol = rdkit.EditableMol(mol)
    emol.AddBond(boron.GetIdx(), oxygens[0].GetIdx(), bond)
    emol.AddBond(boron.GetIdx(), oxygens[1].GetIdx(), bond)

    if del_atoms:
        for a in deleters:
            emol.RemoveAtom(a.GetIdx())

    return emol.GetMol(), 2


def amine3_with_amine3(mol, del_atoms, fg1, fg2):
    """
    Crates bonds between functional groups.

    Parameters
    ----------
    mol : :class:`rdkit.Chem.rdchem.Mol`
        A molecule being assembled.

    del : :class:`bool`
        Toggles if atoms with the ``'del'`` property are deleted.

    fg1 : :class:`int`
        The id of the first functional group which
        is to be joined, as given by the 'fg_id' property.

    fg2 : :class:`int`
        The id of the second functional group which
        is to be joined, as given by the 'fg_id' property.

    Returns
    -------
    :class:`tuple`
        The first element is an :class:`rdkit.Chem.rdchem.Mol`. It is
        the molecule with bonds added between the functional groups.

        The second element is a :class:`int`. It is the number
        of bonds added.

    """

    fgs = {fg1, fg2}
    atoms1, atoms2 = {}, {}
    deleters = []

    for a in mol.GetAtoms():
        if not a.HasProp('fg_id') or a.GetIntProp('fg_id') not in fgs:
            continue

        if a.HasProp('bonder') and a.GetIntProp('fg_id') == fg1:
            atoms1[a.GetSymbol()] = a.GetIdx()

        if a.HasProp('bonder') and a.GetIntProp('fg_id') == fg2:
            atoms2[a.GetSymbol()] = a.GetIdx()

        if a.HasProp('del'):
            deleters.append(a.GetIdx())

    conf = mol.GetConformer()
    n1_pos = np.array([*conf.GetAtomPosition(atoms1['N'])])
    n2_pos = np.array([*conf.GetAtomPosition(atoms2['N'])])

    c1_pos = np.array([*conf.GetAtomPosition(atoms1['C'])])
    c2_pos = np.array([*conf.GetAtomPosition(atoms2['C'])])

    emol = rdkit.EditableMol(mol)

    n_joiner = emol.AddAtom(rdkit.Atom(6))
    n_joiner_pos = (n1_pos + n2_pos) / 2
    nh1 = emol.AddAtom(rdkit.Atom(1))
    nh1_pos = n_joiner_pos + np.array([0, 0, 1])
    nh2 = emol.AddAtom(rdkit.Atom(1))
    nh2_pos = n_joiner_pos + np.array([0, 0, -1])

    nc_joiner1 = emol.AddAtom(rdkit.Atom(6))
    nc_joiner1_pos = (c1_pos + n2_pos) / 2
    nc1h1 = emol.AddAtom(rdkit.Atom(1))
    nc1h1_pos = nc_joiner1_pos + np.array([0, 0, 1])
    nc1h2 = emol.AddAtom(rdkit.Atom(1))
    nc1h2_pos = nc_joiner1_pos + np.array([0, 0, -1])

    nc_joiner2 = emol.AddAtom(rdkit.Atom(6))
    nc_joiner2_pos = (c2_pos + n1_pos) / 2
    nc2h1 = emol.AddAtom(rdkit.Atom(1))
    nc2h1_pos = nc_joiner2_pos + np.array([0, 0, 1])
    nc2h2 = emol.AddAtom(rdkit.Atom(1))
    nc2h2_pos = nc_joiner2_pos + np.array([0, 0, -1])

    single = rdkit.rdchem.BondType.SINGLE
    emol.AddBond(atoms1['N'], n_joiner, single)
    emol.AddBond(atoms2['N'], n_joiner, single)
    emol.AddBond(n_joiner, nh1, single)
    emol.AddBond(n_joiner, nh2, single)

    emol.AddBond(atoms1['C'], nc_joiner1, single)
    emol.AddBond(atoms2['N'], nc_joiner1, single)
    emol.AddBond(nc_joiner1, nc1h1, single)
    emol.AddBond(nc_joiner1, nc1h2, single)

    emol.AddBond(atoms2['C'], nc_joiner2, single)
    emol.AddBond(atoms1['N'], nc_joiner2, single)
    emol.AddBond(nc_joiner2, nc2h1, single)
    emol.AddBond(nc_joiner2, nc2h2, single)

    mol = emol.GetMol()
    conf = mol.GetConformer()
    conf.SetAtomPosition(n_joiner, rdkit_geo.Point3D(*n_joiner_pos))
    conf.SetAtomPosition(nh1, rdkit_geo.Point3D(*nh1_pos))
    conf.SetAtomPosition(nh2, rdkit_geo.Point3D(*nh2_pos))

    conf.SetAtomPosition(nc_joiner1, rdkit_geo.Point3D(*nc_joiner1_pos))
    conf.SetAtomPosition(nc1h1, rdkit_geo.Point3D(*nc1h1_pos))
    conf.SetAtomPosition(nc1h2, rdkit_geo.Point3D(*nc1h2_pos))

    conf.SetAtomPosition(nc_joiner2, rdkit_geo.Point3D(*nc_joiner2_pos))
    conf.SetAtomPosition(nc2h1, rdkit_geo.Point3D(*nc2h1_pos))
    conf.SetAtomPosition(nc2h2, rdkit_geo.Point3D(*nc2h2_pos))

    if del_atoms:
        emol = rdkit.EditableMol(mol)
        for a in reversed(deleters):
            emol.RemoveAtom(a)
        mol = emol.GetMol()

    return mol, 6


# If some functional groups react via a special mechanism not covered
# in by the base "react()" function the function should be placed
# in this dict. The key should be a sorted tuple which holds the name
# of every functional group involved in the reaction along with how
# many such functional groups are invovled.
custom_reactions = {
    FGKey(['boronic_acid', 'diol']): boronic_acid_with_diol,
    FGKey(['diol', 'difluorene']): diol_with_difluorne,
    FGKey(['amine3', 'amine3']): amine3_with_amine3

}

periodic_custom_reactions = {}


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
           del_smarts=[Match(smarts='[$([H][O][#6]~[#6][O][H])]', n=2)]),

    FGInfo(name='difluorene',
           fg_smarts='[F][#6]~[#6][F]',
           bonder_smarts=[Match(smarts='[$([#6]([F])~[#6][F])]', n=2)],
           del_smarts=[Match(smarts='[$([F][#6]~[#6][F])]', n=2)]),

    FGInfo(name='alkyne2',
           fg_smarts='[C]#[C][H]',
           bonder_smarts=[Match(smarts='[$([C]#[C][H])]', n=1)],
           del_smarts=[Match(smarts='[$([H][C]#[C])]', n=1),
                       Match(smarts='[$([C](#[C])[H])]', n=1)]),

    FGInfo(name='amine3',
           fg_smarts='[N]([H])([H])[#6]~[#6]([H])~[#6]([H])',
           bonder_smarts=[
            Match(smarts='[$([N]([H])([H])[#6]~[#6]([H])~[#6]([H]))]', n=1),
            Match(smarts='[$([#6]([H])(~[#6]([H]))~[#6][N]([H])[H])]', n=1),
            ],
           del_smarts=[
            Match(smarts='[$([H][N]([H])[#6]~[#6]([H])~[#6]([H]))]', n=2),
            Match(smarts='[$([H][#6](~[#6]([H]))~[#6][N]([H])[H])]', n=1)])

    )

double = rdkit.rdchem.BondType.DOUBLE
triple = rdkit.rdchem.BondType.TRIPLE
bond_orders = {
    FGKey(['amine', 'aldehyde']): double,
    FGKey(['amide', 'aldehyde']): double,
    FGKey(['nitrile', 'aldehyde']): double,
    FGKey(['amide', 'amine']): double,
    FGKey(['terminal_alkene', 'terminal_alkene']): double,
    FGKey(['alkyne2', 'alkyne2']): triple}
