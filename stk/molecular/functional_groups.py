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

    bonder_smarts : :class:`str`
        A SMARTS string which matches the atom on the functional group
        which forms bonds during reactions.

    del_smarts : :class:`str`
        A SMARTS string, which matches the atoms removed when the
        functional group reacts.

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

        bonder_smarts : :class:`str`
            A SMARTS string which matches the atom on the functional
            group which forms bonds during reactions.

        del_smarts : :class:`str`
            A SMARTS string, which matches the atoms removed when the
            functional group reacts.

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


functional_groups = (

                FGInfo("amine",
                       "[N]([H])[H]",
                       "[$([N]([H])[H])]",
                       "[$([H][N][H])].[$([H][N][H])]"),

                FGInfo("aldehyde",
                       "[C](=[O])[H]",
                       "[$([C](=[O])[H])]",
                       "[$([O]=[C][H])]"),

                FGInfo("carboxylic_acid",
                       "[C](=[O])[O][H]",
                       "[$([C](=[O])[O][H])]",
                       "[$([H][O][C](=[O]))][O]"),

                FGInfo("amide",
                       "[C](=[O])[N]([H])[H]",
                       "[$([C](=[O])[N]([H])[H])]",
                       "[$([N]([H])([H])[C](=[O]))]([H])[H]"),

                FGInfo("thioacid",
                       "[C](=[O])[S][H]",
                       "[$([C](=[O])[O][H])]",
                       "[$([H][O][C](=[O]))][S]"),

                FGInfo("alcohol",
                       "[O][H]",
                       "[$([O][H])]",
                       "[$([H][O])]"),

                FGInfo("thiol",
                       "[S][H]",
                       "[$([S][H])]",
                       "[$([H][S])]"),

                FGInfo("bromine",
                       "*[Br]",
                       "[$(*[Br])]",
                       "[$([Br]*)]"),

                FGInfo("iodine",
                       "*[I]",
                       "[$(*[I])]",
                       "[$([I]*)]"),

                FGInfo("nitrile",
                       "[C][C]#[N]",
                       "[$([C]([H])([H])[C]#[N])]",
                       "[$([H][C][H])].[$([H][C][H])]"),

                FGInfo('alkyne',
                       '[C]#[C][H]',
                       '[$([C]([H])#[C])]',
                       '[$([H][C]#[C])]'),

                FGInfo('terminal_alkene',
                       '[C]=[C]([H])[H]',
                       '[$([C]=[C]([H])[H])]',
                       ('[$([H][C]([H])=[C])].'
                        '[$([H][C]([H])=[C])].'
                        '[$([C](=[C])([H])[H])]')),

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
