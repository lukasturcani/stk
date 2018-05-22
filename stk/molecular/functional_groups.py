"""
Defines tools for dealing with functional groups and their reactions.

.. _`adding functional groups`:

Extending stk: Adding  more functional groups.
----------------------------------------------

If ``stk`` is to incorporate a new functional group, a new
:class:``FGInfo`` instance should be added to
:data:`functional_groups`.

Adding a new ``FGInfo`` instance to :data:`functional_groups` will
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

import rdkit.Chem.AllChem as rdkit
from collections import Counter


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
        A :class:`tuple` of the form

        .. code-block:: python

            bonder_smarts = (('[$([N]([H])[H])]', 1),
                             ('[$([H][N][H])], 1))

        Each string is SMARTS string which matches an atom in the
        functional group which is to be tagged as ``'bonder'``. The
        number represents how many matched atoms should be tagged, per
        functional group.

        In the example, ``('[$([N]([H])[H])]', 1)`` matches the
        nitrogen atom in the amine functional group. The ``1`` means
        that 1 nitrogen atom per functional group will be tagged as
        ``'bonder'``. The second :class:`tuple`,
        ``('[$([H][N][H])], 1))``, matches the hydrogen atom in the
        amine functional group. Because the number in the
        :class:`tuple` is ``1``, only 1 of hydrogen atom per
        amine functional group will be tagged ``'bonder'``. If instead
        the tuple was ``('[$([H][N][H])], 2)``, then both of the
        hydrogen atoms in the functional group would be tagged.

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

        The seoncd element is a :class:`int`. It is the number
        of bonds added.

    """

    names = [fg_name(mol, fg) for fg in fgs]
    reaction_key = tuple(sorted(Counter(names).items()))
    if reaction_key in custom_reactions:
        return custom_reactions[reaction_key](mol, del_atoms, *fgs)

    emol = rdkit.EditableMol(mol)

    bonders = []
    for atom in mol.GetAtoms():
        if not (atom.HasProp('fg_id') and atom.GetIntProp('fg_id') in fgs):
            continue
        if atom.HasProp('bonder'):
            bonders.append(atom.GetIdx())

    bond = bond_orders.get(frozenset(names), rdkit.rdchem.BondType.SINGLE)
    bonder1, bonder2 = bonders
    emol.AddBond(bonder1, bonder2, bond)

    for atom in reversed(mol.GetAtoms()):
        if not (atom.HasProp('fg_id') and atom.GetIntProp('fg_id') in fgs):
            continue

        if atom.HasProp('del') and del_atoms:
            emol.RemoveAtom(atom.GetIdx())

    return emol.GetMol(), 1


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


# If some functional groups react via a special mechanism not covered
# in by the base "react()" function the function should be placed
# in this dict. The key should be a sorted tuple which holds the name
# of every functional group involved in the reaction along with how
# many such functional groups are invovled.
custom_reactions = {
    tuple(sorted((('boronic_acid', 1), ('diol', 1)))): boronic_acid_with_diol}


_amine = FGInfo(name="amine",
                fg_smarts="[N]([H])[H]",
                bonder_smarts=[("[$([N]([H])[H])]", 1)],
                del_smarts=[("[$([H][N][H])]", 1)])


_aldehyde = FGInfo(name="aldehyde",
                   fg_smarts="[C](=[O])[H]",
                   bonder_smarts=[("[$([C](=[O])[H])]", 1)],
                   del_smarts=[("[$([O]=[C][H])]", 1)])

_carboxylic_acid = FGInfo(name="carboxylic_acid",
                          fg_smarts="[C](=[O])[O][H]",
                          bonder_smarts=[("[$([C](=[O])[O][H])]", 1)],
                          del_smarts=[("[$([H][O][C](=[O]))]", 1),
                                      ("[$([O]([H])[C](=[O]))]", 1)])

_amide = FGInfo(name="amide",
                fg_smarts="[C](=[O])[N]([H])[H]",
                bonder_smarts=[("[$([C](=[O])[N]([H])[H])]", 1)],
                del_smarts=[("[$([N]([H])([H])[C](=[O]))]", 1),
                            ("[$([H][N]([H])[C](=[O]))]", 2)])

_thioacid = FGInfo(name="thioacid",
                   fg_smarts="[C](=[O])[S][H]",
                   bonder_smarts=[("[$([C](=[O])[S][H])]", 1)],
                   del_smarts=[("[$([H][S][C](=[O]))]", 1),
                               ("[$([S]([H])[C](=[O]))]", 1)])

_alcohol = FGInfo(name="alcohol",
                  fg_smarts="[O][H]",
                  bonder_smarts=[("[$([O][H])]", 1)],
                  del_smarts=[("[$([H][O])]", 1)])

_thiol = FGInfo(name="thiol",
                fg_smarts="[S][H]",
                bonder_smarts=[("[$([S][H])]", 1)],
                del_smarts=[("[$([H][S])]", 1)])

_bromine = FGInfo(name="bromine",
                  fg_smarts="*[Br]",
                  bonder_smarts=[("[$(*[Br])]", 1)],
                  del_smarts=[("[$([Br]*)]", 1)])

_iodine = FGInfo(name="iodine",
                 fg_smarts="*[I]",
                 bonder_smarts=[("[$(*[I])]", 1)],
                 del_smarts=[("[$([I]*)]", 1)])

_alkyne = FGInfo(name='alkyne',
                 fg_smarts='[C]#[C][H]',
                 bonder_smarts=[('[$([C]([H])#[C])]', 1)],
                 del_smarts=[('[$([H][C]#[C])]', 1)])

_terminal_alkene = FGInfo(name='terminal_alkene',
                          fg_smarts='[C]=[C]([H])[H]',
                          bonder_smarts=[('[$([C]=[C]([H])[H])]', 1)],
                          del_smarts=[('[$([H][C]([H])=[C])]', 2),
                                      ('[$([C](=[C])([H])[H])]', 1)])

functional_groups = (

                ,

,

,

,

,

,

,

,

,


,

,

                FGInfo('boronic_acid',
                       '[B]([O][H])[O][H]',
                       '[$([B]([O][H])[O][H])]',
                       ('[$([O]([H])[B][O][H])].'
                        '[$([H][O][B][O][H])]')),

                # This amine functional group only deletes one of the
                # hydrogen atoms when a bond is formed.
                FGInfo("amine2",
                       "[N]([H])[H]",
                       "[$([N]([H])[H])]",
                       "[$([H][N][H])]"),

                FGInfo("secondary_amine",
                       "[H][N]([#6])[#6]",
                       "[$([N]([H])([#6])[#6])]",
                       "[$([H][N]([#6])[#6])]"),

                FGInfo('diol',
                       '[H][O][#6][#6][O][H]',
                       ('[$([O]([H])[#6][#6][O][H])].'
                        '[$([O]([H])[#6][#6][O][H])]'),
                       ('[$([H][O][#6][#6][O][H])].'
                        '[$([H][O][#6][#6][O][H])]'))


                    )

double = rdkit.rdchem.BondType.DOUBLE
bond_orders = {
    frozenset(('amine', 'aldehyde')): double,
    frozenset(('amide', 'aldehyde')): double,
    frozenset(('nitrile', 'aldehyde')): double,
    frozenset(('amide', 'amine')): double,
    frozenset(('terminal_alkene', )): double}
