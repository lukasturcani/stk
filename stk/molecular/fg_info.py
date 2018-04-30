"""
Defines tools for dealing with functional groups and their reactions.

.. _`adding functional groups`:

Extending stk: Adding  more functional groups.
----------------------------------------------

If ``stk`` is to incorporate a new functional group, a new
:class:``FGInfo`` instance should be added to
:data:`functional_groups`. This is a :class:`dict` defined in this
module.

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
with a double other than a single, the names of the functional groups
should be added to :data:`bond_orders`, along with the desired bond
order.

Supporting complex reactions.
.............................

During assembly, two functional groups are provided to
:func:`join_fgs`. By default, placing an :class:`FGInfo` instance into
:data:`functional_groups` will result in the creation of a single bond
between the atoms tagged as 'bonder' in the two functional groups.
In addtion, any atoms tagged as 'del' will be removed. The bond order
of the created bond can be modified by editing :data:`bond_orders`.

However, some reactions cannot be described by a simple combination of
adding a bond while deleting some existing atoms. For example, consider
the aldol reaction:

    CH3C(=O)CH3 + CH3C(=O)CH3 --> CH3(=O)CH2C(OH)(CH3)CH3

Here a ketone is converted into an alcohol. In order to support more
complex conversions, a specific function needs to be defined which
modifies the molecule as desired. The function then needs
to be added to :data:`custom_joins`. See
:func:`join_boronic_acid_with_diol`
as an example.

"""

import rdkit.Chem.AllChem as rdkit


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

    fg : :class:`str`
        The id of a functional group as given by the 'fg_id' property.

    Returns
    -------
    :class:`str`
        The name of a functional group.

    """

    for atom in mol.GetAtoms():
        if atom.HasProp('fg_id') and atom.GetProp('fg_id') == fg:
            return atom.GetProp('fg')


def join_fgs(mol, fg1, fg2):
    """
    Crates bonds between functional groups.

    Except for cases in :data:`custom_joins`, this function operates
    on atoms where the 'fg_id' tag is either `fg1` or `fg2`. All such
    atoms with the tag 'del' are deleted and the ones with the tag
    'bonder' are joined with a bond. The bond order is ``1`` unless
    listed otherwise in :data:`bond_orders`.

    Parameters
    ----------
    mol : :class:`rdkit.Chem.rdchem.Mol`
        A molecule being assembled.

    fg1 : :class:`str`
        The id of the first functional group which
        is to be joined, as given by the 'fg_id' property.

    fg2 : :class:str`
        The id of the second functional group which
        is to be joined, as given by the 'fg_id' property.

    Returns
    -------
    :class:`rdkit.Chem.rdchem.Mol`
        The molecule with bonds added between the functional groups.

    """

    ids = {fg1, fg2}
    name1, name2 = fg_name(mol, fg1), fg_name(mol, fg2)
    join_key = frozenset((name1, name2))
    if join_key in custom_joins:
        return custom_joins[join_key](fg1, fg2)

    emol = rdkit.EditableMol(mol)
    bonders = []
    for atom in reversed(mol.GetAtoms()):
        if not (atom.HasProp('fg_id') and atom.GetProp('fg_id') in ids):
            continue

        if atom.HasProp('del'):
            emol.RemoveAtom(atom.GetIdx())

        if atom.HasProp('bonder'):
            bonders.append(atom)

    bond = bond_orders.get(join_key, rdkit.rdchem.BondType.SINGLE)
    bonder1, bonder2 = bonders
    emol.AddBond(bonder1.GetIdx(), bonder2.GetIdx(), bond)
    return emol.GetMol()


def join_boronic_acid_with_diol(mol, fg1, fg2):
    """
    Crates bonds between functional groups.

    Parameters
    ----------
    mol : :class:`rdkit.Chem.rdchem.Mol`

    fg1 : :class:`str`
        The id of the first functional group which
        is to be joined, as given by the 'fg_id' property.

    fg2 : :class:str`
        The id of the second functional group which
        is to be joined, as given by the 'fg_id' property.

    Returns
    -------
    None : :class:`NoneType`

    """

    ...


custom_joins = {
    frozenset(('boronic_acid', 'diol')): join_boronic_acid_with_diol}


functional_groups = {

                'amine': FGInfo("amine",
                                "[N]([H])[H]",
                                "[$([N]([H])[H])]",
                                "[$([H][N][H])].[$([H][N][H])]"),

                'aldehyde': FGInfo("aldehyde",
                                   "[C](=[O])[H]",
                                   "[$([C](=[O])[H])]",
                                   "[$([O]=[C][H])]"),

                'carboxylic_acid': FGInfo("carboxylic_acid",
                                          "[C](=[O])[O][H]",
                                          "[$([C](=[O])[O][H])]",
                                          "[$([H][O][C](=[O]))][O]"),

                'amide': FGInfo("amide",
                                "[C](=[O])[N]([H])[H]",
                                "[$([C](=[O])[N]([H])[H])]",
                                "[$([N]([H])([H])[C](=[O]))]([H])[H]"),

                'thioacid': FGInfo("thioacid",
                                   "[C](=[O])[S][H]",
                                   "[$([C](=[O])[O][H])]",
                                   "[$([H][O][C](=[O]))][S]"),

                'alcohol': FGInfo("alcohol",
                                  "[O][H]",
                                  "[$([O][H])]",
                                  "[$([H][O])]"),

                'thiol': FGInfo("thiol",
                                "[S][H]",
                                "[$([S][H])]",
                                "[$([H][S])]"),

                'bromine': FGInfo("bromine",
                                  "*[Br]",
                                  "[$(*[Br])]",
                                  "[$([Br]*)]"),

                'iodine': FGInfo("iodine",
                                 "*[I]",
                                 "[$(*[I])]",
                                 "[$([I]*)]"),

                'nitrile': FGInfo("nitrile",
                                  "[C][C]#[N]",
                                  "[$([C]([H])([H])[C]#[N])]",
                                  "[$([H][C][H])].[$([H][C][H])]"),

                'alkyne': FGInfo('alkyne',
                                 '[C]#[C][H]',
                                 '[$([C]([H])#[C])]',
                                 '[$([H][C]#[C])]'),

                'terminal_alkene': FGInfo('terminal_alkene',
                                          '[C]=[C]([H])[H]',
                                          '[$([C]=[C]([H])[H])]',
                                          ('[$([H][C]([H])=[C])].'
                                           '[$([H][C]([H])=[C])].'
                                           '[$([C](=[C])([H])[H])]')),

                'boronic_acid': FGInfo('boronic_acid',
                                       '[B]([O][H])[O][H]',
                                       '[$([B]([O][H])[O][H])]',
                                       ('[$([O]([H])[B][O][H])].'
                                        '[$([H][O][B][O][H])]')),

                # This amine functional group only deletes one of the
                # hydrogen atoms when a bond is formed.
                'amine2': FGInfo("amine2",
                                 "[N]([H])[H]",
                                 "[$([N]([H])[H])]",
                                 "[$([H][N][H])]"),

                'secondary_amine': FGInfo("secondary_amine",
                                          "[H][N]([#6])[#6]",
                                          "[$([N]([H])([#6])[#6])]",
                                          "[$([H][N]([#6])[#6])]")


                    }

double = rdkit.rdchem.BondType.DOUBLE
bond_orders = {
    frozenset(('amine', 'aldehyde')): double,
    frozenset(('amide', 'aldehyde')): double,
    frozenset(('nitrile', 'aldehyde')): double,
    frozenset(('amide', 'amine')): double,
    frozenset(('terminal_alkene', 'terminal_alkene')): double}
