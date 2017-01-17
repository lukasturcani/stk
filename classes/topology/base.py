from collections import Counter
import rdkit
import rdkit.Chem as chem
from scipy.spatial.distance import euclidean

from ..fg_info import double_bond_combs


class Topology:
    """
    Represents the topology of an assembled molecule.
    
    The ``Topology`` class is concerned with how individual building 
    blocks are placed and connected in space to form an assembled 
    molecule used by MMEA. It also takes care of assembling these
    molecules. The general process of building macromolecules is 
    discussed in detail in the `build` method documentation.
    
    This class directly defines any operations and attributes that are
    needed by any topology, be it a tetrahedron, octahedron or even a
    polymer. However, this class is not used directly by MMEA. It is
    intended to be inherited from. Any individual within MMEA will have
    a `topology` attribute which refers to an instance of a class 
    derived from this. Derived classes of ``Topology`` define things
    specific to that one topology. For example, each derived class must 
    define which ``pair_up`` function defined in the ``Topology`` class 
    it uses for pairing up its heavy atoms. This is done by placing
    the function in the `pair_up` attribute of the derived class. See
    the included derived classes as examples. In addition, each class 
    derived from ``Topology`` must define methods which place building 
    blocks in the correct positions, such as chosen edges or vertices.

    Instances of this class should not be created directly. Only via a
    derived class. Multiple inheritance can be useful when creating 
    derived classes. For example, all topologies describing cages will
    share some characteristics. This means a class ``CageTopology`` can
    be created which holds all information required by all cage 
    topologies. This class, ``CageTopology``, will inherit ``Topology``. 
    A specific cage topology such as ``FourPlusSix`` or 
    ``EightPlusTwelve`` will then inherit ``CageTopology`` and add any 
    information specific to that one topology.
    
    Extending MMEA: Adding new topologies
    -------------------------------------
    > Cages
    To add a new cage topology a new class should be created, named
    after the topology. This class should inhertic the ``CageTopology``
    class. This will give access to various methods which are necessary
    for dealing with any cage molecule. See the documenation of 
    ``CageTopology`` for more details.
    
    The new class will only need to have five class attributes added:
        1) a list called `vertices` 
        2) a list called `edges`
        3) `n_windows` which holds the number of windows the cage 
           topology has
        4) `n_window_types` which holds the number of different window
           types. For example, if `n_window_types` is 2 then the 
           topology will have two kinds of windows, each with a 
           different expected size even in a perfectly symmetrical case. 
           Windows of the same type are expected to be of the same size.
        
    The `vertices` list holds instances of the class ``Vertex``. Each
    instance represents a vertex of a cage and needs to be initialized
    with the coordinates of that vertex. Vertices of a cage are where
    building-blocks* of cages are placed.
    
    The `edges` list holds instances of the class ``Edge``. Each
    instance represents an edge of a cage and needs to be initialized
    with two instnaces of the ``Vertex`` class. The ``Vertex`` instances
    should be held in the `vertices` list mentioned above. These are 
    the two vertices which the edge connects. Linkers of cages are 
    placed on edges. The edge instances automatically derive their 
    positions from the vertices supplied during initialization.

    The vertices need to be positioned such that the center of the
    topology is at the origin.
    
    Attributes
    ----------
    macro_mol : MacroMolecule
        The ``MacroMolecule`` instance which has this topology. This 
        gives easy access to the macromolecule's attributes to the 
        ``Topology``instance.
        
    bonds_made : int
        The number of bonds created during assembly. This should be
        incremened for each new bond made during assembly. Used in some
        fitness functions.
    
    bb_counter : Counter
        A counter keeping track of how much of each building block was
        used during assembly. The ``StructUnit`` instance acts as the
        counter key.
        
        
    """
    
    def __init__(self, macro_mol):
        self.macro_mol = macro_mol
        self.bonds_made = 0
        self.bb_counter = Counter()         
        
    def build(self):
        """
        Assembles rdkit instances of the macromolecules.
        
        To carry out `build` an instance of a class derived from 
        ``Topology`` must be used. This is because instances of such
        classes define a `pair_up` attribute during initialization.
        (This should be done by default, not passed as an argument to
        the initializer.) The `pair_up` attribute holds the pair up 
        function defined within ``Topology``, which should be used by
        `build`. (`pair_up` is used within the `join_mols` subroutine of
        `build`.)
        
        Modifies
        --------
        self.macro_mol.mol
            Adds an rdkit instance of the assembled molecule into this 
            attribute.
            
        self.bonds_made : int
            This counter is updated with each bond made during assembly.            
            
        self.bb_counter : Counter
            The counter is updated with the number of building blocks of
            each StructUnit that make up the MacroMolecule.               
            
        Returns
        -------
        None : NoneType
        
        """

        self.place_mols()
        self.join_mols()

    def determine_bond_type(self, atom1_id, atom2_id):
        """
        Returns the bond order to be formed between the atoms.
        
        Some atoms will need to have a double bond created between them.
        This is defined in the `double_bond_combs` list. If the
        atom ids provided as paramters belong to atoms of elements
        found in this list, the rdkit double bond type will be returned.
        If not the rdkit single bond type will be returned. These types
        are needed when adding bonds using ``EditableMol`` instances.
        
        Parameters
        ----------
        atom1_id : int
            The id number of the first atom.
        
        atom2_id : int
            The id number of the second atom.
        
        Returns
        -------
        rdkit.Chem.rdchem.BondType.SINGLE
            If the combination of heavy atoms passed as arguments is not 
            in `double_bond_combs`.
        
        rdkit.Chem.rdchem.BondType.DOUBLE
            If the combination of heavy atoms passed as arguments is in
            `double_bond_combs`.
            
        """
        
        # Get the functional groups of the of the atoms whose atom ids 
        # were supplied as arguments. If the groups form a tuple in 
        # `double_bond_combs` return a rdkit double bond type. If they 
        # do not, return a rdkit single bond type.
        
        atom1 = self.macro_mol.mol.GetAtomWithIdx(atom1_id)
        atom1_grp = atom1.GetProp('fg')
        atom2 = self.macro_mol.mol.GetAtomWithIdx(atom2_id)
        atom2_grp= atom2.GetProp('fg')      
        
        double_bond_present = ((atom1_grp, atom2_grp) == tup or 
                               (atom2_grp, atom1_grp) == tup for 
                               tup in double_bond_combs)
        
        if any(double_bond_present):
            return rdkit.Chem.rdchem.BondType.DOUBLE
        else:
            return rdkit.Chem.rdchem.BondType.SINGLE

