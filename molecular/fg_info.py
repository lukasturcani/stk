"""
A module for the FGInfo class.

Extending MMEA: Adding  more functional groups.
-----------------------------------------------

If MMEA is to incorporate a new functional group, a new ``FGInfo``
instance should be added to the `functional_groups` list.

Adding a new ``FGInfo`` instance to `functional_groups` will allow the
``Topology.join_mols()`` method to connect this functional group to
(all) others during assembly. Nothing except adding this instance
should be necessary in order to incorporate new functional groups.

Note that when adding SMARTS if you want to make a SMARTS that targets
an atom in an environment, for example a bromine connected to a carbon

    [$([Br][C])]

The atom you are targeting needs to be written first. The above SMARTS
works but

    [$([C][Br])]

does not.

If this new functional group is to connect to another functional group
with a double bond during assembly, the names of the functional groups
should be added to the `double_bond_combs` list. The order in
which they are placed in the tuple does not matter. Again, this is all
that needs to be done for MMEA to create double bonds between given
functional groups.

"""


class FGInfo:
    """
    Contains key information about functional groups.

    The point of this class is to register which atoms of a functional
    group form bonds, and which are deleted during assembly of
    macromolecules.

    Attributes
    ----------
    name : str
        The name of the functional group.

    fg_smarts : str
        A SMARTS string which matches the functional group.

    bonder_smarts :
        A SMARTS string which matches the atom on the functional group
        which forms bonds during reactions.

    del_smarts : str
        A SMARTS string, which matches the atoms removed when the
        functional group reacts.

    """

    __slots__ = ['name', 'fg_smarts', 'bonder_smarts', 'del_smarts']

    def __init__(self, name, fg_smarts, bonder_smarts, del_smarts):
        self.name = name
        self.fg_smarts = fg_smarts
        self.bonder_smarts = bonder_smarts
        self.del_smarts = del_smarts


functional_groups = [

                FGInfo("amine",
                       "[N]([H])[H]",
                       "[$([N]([H])[H])]",
                       "[$([H][N][H])].[$([H][N][H])]"),

                FGInfo("aldehyde", "[C](=[O])[H]",
                                   "[$([C](=[O])[H])]",
                                   "[$([O]=[C][H])]"),

                FGInfo("carboxylic_acid",
                       "[C](=[O])[O][H]",
                       "[$([C](=[O])[O][H])]",
                       "[$([H][O][C](=[O]))][O]"),

                FGInfo("amide", "[C](=[O])[N]([H])[H]",
                       "[$([C](=[O])[N]([H])[H])]",
                       "[$([N]([H])([H])[C](=[O]))]([H])[H]"),

                FGInfo("thioacid",
                       "[C](=[O])[S][H]",
                       "[$([C](=[O])[O][H])]",
                       "[$([H][O][C](=[O]))][S]"),

                FGInfo("alcohol", "[O][H]",
                                  "[$([O][H])]",
                                  "[$([H][O])]"),

                FGInfo("thiol", "[S][H]",
                                "[$([S][H])]",
                                "[$([H][S])]"),

                FGInfo("bromine", "*[Br]",
                                  "[$(*[Br])]",
                                  "[$([Br]*)]"),

                FGInfo("iodine", "*[I]",
                                 "[$(*[I])]",
                                 "[$([I]*)]"),

                # This amine functional group only deletes one of the
                # hydrogen atoms when a bond is formed.
                FGInfo("amine2",
                       "[N]([H])[H]",
                       "[$([N]([H])[H])]",
                       "[$([H][N][H])]"),

                FGInfo("secondary_amine",
                       "[H][N]([#6])[#6]",
                       "[$([N]([H])([#6])[#6])]",
                       "[$([H][N]([#6])[#6])]")


                    ]

double_bond_combs = [('amine', 'aldehyde'),
                     ('amide', 'aldehyde'),
                     ('amide', 'amine')]
