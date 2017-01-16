"""
A module for the FGInfo class.

Extending MMEA: Adding  more functional groups.
-----------------------------------------------

If MMEA is to incorporate a new functional group, a new ``FGInfo`` 
instance should be added to the `functional_groups` list. 

Adding a new ``FGInfo`` instace to `functional_groups` will allow the
``Topology.join_mols()`` method to connect this functional group to 
(all) others during assembly. Nothing except adding this instance should 
be necessary in order to incorporate new functional groups.

If this new functional group is to connect to another functional group 
with a double bond during assembly, the names of the functional groups
should be added to the `double_bond_combs` list. The order in 
which they are placed in the tuple does not matter. Again, this is all 
that needs to be done for MMEA to create double bonds between given 
functional groups.  

"""

from ..convenience_tools import bond_dict

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

    target_smarts : 
        A SMARTS string which matches the atom on the functional group
        which forms bonds during reactions.        
        
    del_smarts : str
        A SMARTS string, which match the atoms removed the functional 
        group reacts.
    
    """
    
    __slots__ = ['name', 'fg_smarts', 'target_smarts', 'del_smarts'] 
    
    def __init__(self, name, fg_smarts, target_smarts, del_smarts):
         self.name = name
         self.fg_smarts = fg_smarts
         self.target_smarts = target_smarts
         self.del_smarts = del_smarts

functional_groups = [

                FGInfo("amine", 
                       "[N]([H])[H]", 
                       "[$([N]([H])[H]);$([N])]", 
                       "[$([H][N][H]);$([H])].[$([H][N][H]);$([H])]"),                        

                FGInfo("aldehyde", "[C](=[O])[H]", 
                                   "[$([C](=[O])[H]);$([C])]", 
                                   "[$([O]=[C][H]);$([O])]"), 

                FGInfo("carboxylic_acid", 
                       "[C](=[O])[O][H]", 
                       "[$([C](=[O])[O][H]);$([C])]", 
                       "[$([H][O][C](=[O]));$([H])][O]"),

                FGInfo("amide", "[C](=[O])[N]([H])[H]",
                       "[$([C](=[O])[N]([H])[H]);$([C])]", 
                       "[$([N]([H])([H])[C](=[O]));$([N])]([H])[H]"),

                FGInfo("thioacid", 
                       "[C](=[O])[S][H]", 
                       "[$([C](=[O])[O][H]);$([C])]", 
                       "[$([H][O][C](=[O]));$([H])][S]"),

                FGInfo("alcohol", "[O][H]", 
                                  "[$([O][H]);$([O])]",
                                  "[$([H][O]);$([H])]"),

                FGInfo("thiol", "[S][H]", 
                                "[$([S][H]);$([S])]",   
                                "[$([H][S]);$([H])]")
                           
                    ]

double_bond_combs = [('amine','aldehyde'), 
                     ('amide','aldehyde'), 
                     ('amide','amine')]
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     