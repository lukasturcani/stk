from collections import namedtuple

from ..convenience_tools import bond_dict

class FGInfo:
    """
    Contains key information for functional group substitutions.
    
    The point of this class is to register which atoms are substituted
    for which, when an atom in a functional group is substituted with a 
    heavy metal atom. If MMEA is to incorporate a new functional group, 
    a new ``FGInfo`` instance should be added to the 
    `functional_group_list` class attribute of ``FGInfo``. 
    
    Adding a new ``FGInfo`` instace to `functional_group_list` will 
    allow the `Topology.join_mols` method to connect this functional 
    group to (all) others during assembly. Nothing except adding this
    instance should need to be done in order to incorporate new 
    functional groups.
    
    If this new functional group is to connect to another functional 
    group with a double bond during assembly, the symbols of the heavy 
    atoms of both functional groups should be added to the 
    `double_bond_combs` class attribute. The order in which the heavy 
    symbols are placed in the tuple does not matter. Again, this is all
    that needs to be done for MMEA to create double bonds between
    certain functional groups.  
    
    Class attributes
    ----------------
    functional_groups_list : list of FGInfo instances
        This list holds all ``FGInfo`` instances used by MMEA. If a new
        functional group is to be used by MMEA, a new ``FGInfo`` 
        instance must be added to this list.
        
    double_bond_combs : list of tuples of strings
        When assembly is carried out, if the heavy atoms being joined
        form a tuple in this list, they will be joined with a double
        rather than single bond. If a single bond is desired there is no
        need to change this variable.
        
    heavy_symbols : set of str
        A set of all the heavy symbols used by ``FGInfo`` instances in 
        MMEA. This set updates itself automatically. There is no need to
        modify it when changes are made to any part of MMEA.
        
    heavy_atomic_nums : set of ints
        A set of all atomic numbers of heavy atoms used by ``FGInfo``
        instances in MMEA. This set updates itself automatically. There
        is no need to modify it when chagnes are made to any part of
        MMEA.

    Attributes
    ----------
    name : str
        The name of the functional group.
    
    smarts_start : str
        A ``SMARTS`` string describing the functional group before 
        substitution by a heavy atom.
        
    del_tags : list of DelAtom instances
        Every member of this list represents an atom on the functional
        group which should be deleted during assembly. One atom in each
        functional group is removed for each list member.
    
    target_atomic_num : int
        The atomic number of the atom being substituted by a heavy atom.
    
    heavy_atomic_num : int
        The atomic number of the heavy atom which replaces the target 
        atom.
    
    target_symbol : str
        The atomic symbol of the atom, being substituted by a heavy 
        atom.       
    
    heavy_symbol : str
        The atomic symbol of the heavy atom which replaces the target 
        atom.
    
    """
    
    __slots__ = ['name', 'smarts_start', 'del_tags', 
                 'target_atomic_num', 'heavy_atomic_num', 
                 'target_symbol', 'heavy_symbol'] 
    
    def __init__(self, name, smarts_start, del_tags, target_atomic_num, 
                 heavy_atomic_num, target_symbol, heavy_symbol):
         self.name = name
         self.smarts_start = smarts_start
         self.del_tags = del_tags
         self.target_atomic_num = target_atomic_num
         self.heavy_atomic_num = heavy_atomic_num
         self.target_symbol = target_symbol
         self.heavy_symbol = heavy_symbol

# An atom is deleted based on what type of bond connects it to the
# substituted functional group atom. The element of the atom is ofcourse
# a factor as well. When both of these are satisfied the atom is
# removed. The ``DelAtom`` class conveniently stores this information.
# Bond type is an rdkit bond type (see the bond dictionary above for
# the two possible values it may take) and atomic num in an integer.
DelAtom = namedtuple('DelAtom', ['bond_type', 'atomic_num'])

FGInfo.functional_group_list = [
                        
    FGInfo("aldehyde", "C(=O)[H]", [ DelAtom(bond_dict['2'], 8) ], 
                                                       6, 39, "C", "Y"), 
    
    FGInfo("carboxylic_acid", "C(=O)O[H]", 
           [ DelAtom(bond_dict['1'], 8) ], 6, 40, "C", "Zr"),
    
    FGInfo("amide", "C(=O)N([H])[H]", [ DelAtom(bond_dict['1'], 7) ], 
                                                      6, 41, "C", "Nb"),
    
    FGInfo("thioacid", "C(=O)S[H]", [ DelAtom(bond_dict['1'], 16) ], 
                                                      6, 42, "C", "Mo"),
    
    FGInfo("alcohol", "O[H]", [], 8, 43, "O", "Tc"),
    FGInfo("thiol", "[S][H]", [], 16, 44, "S", "Ru"),
    FGInfo("amine", "[N]([H])[H]", [], 7, 45, "N", "Rh"),       
    FGInfo("nitroso", "N=O", [], 7, 46, "N", "Pd"),
    FGInfo("boronic_acid", "[B](O[H])O[H]", [], 5, 47, "B", "Ag")
                             
                             ]

FGInfo.double_bond_combs = [("Rh","Y"), ("Nb","Y"), ("Mb","Rh")]

FGInfo.heavy_symbols = {x.heavy_symbol for x 
                                        in FGInfo.functional_group_list}
                        
FGInfo.heavy_atomic_nums = {x.heavy_atomic_num for x 
                                        in FGInfo.functional_group_list}