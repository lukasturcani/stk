from MMEA.classes.topology import Topology

class MolFileData:
    """
    Used to store data collected from a ``.mol`` file during assembly.
    
    This class is used during assembly to keep track of the content
    which should be written to the final ``.mol`` file of the assembled
    molecule.
    
    Attributes
    ----------
    mol_file_content : str
        This attribute holds a string with all of the ``.mol`` file 
        content. This content is edited to include the additional bonds
        need to assemble a larger molecule.
    
    count_line : str
        Represents the ``count line`` in the ``.mol`` file. This line
        holds the atom and bond number. The bond number needs to be
        edited during course of assembly.
        
    at_num : str
        This attribute holds the number of atoms defined in the ``.mol``
        file. It is held a string for convenience as it is read from a
        string and later written to one.

    bond_number : int
        Holds the number of bonds in the ``.mol`` file. This will be 
        written to the final ``.mol`` file generated during assembly.
        As a result it will be added to for each new bond created during
        assembly.

    atom1_list : list of ``Atom`` instances
        All the heavy atoms representing one of the functional groups in
        a ``.mol`` file are placed into ``Atom`` instances which are 
        stored here. Each functional group has its own such 
        attribute/list.
        
    atom2_list : list of ``Atom`` instances      
        All the heavy atoms representing one of the functional groups in
        a ``.mol`` file are placed into ``Atom`` instances which are 
        stored here. Each functional group has its own such 
        attribute/list.
    
    """
    
    __slots__ = ['mol_file_content', 'count_line', 'at_num', 
                 'bond_number', 'atom1_list', 'atom2_list']
    
    def __init__(self, mol_file_content, count_line, at_num, 
                 bond_number, atom1_list, atom2_list):
         self.mol_file_content = mol_file_content
         self.count_line = count_line
         self.at_num = at_num
         self.bond_number = bond_number
         self.atom1_list = atom1_list
         self.atom2_list = atom2_list
    
class Atom(object):
    """
    Used to store and manipulate data collected from a ``.mol`` file.    

    This class is used during assembly to keep track of which atoms need
    to be linked together by writing bonds to the ``.mol`` file. It is 
    also used to find the distances between heavy atoms in the ``.mol``
    file. This allows the creation of bonds between the correct atoms. A
    full list of operations performed by this class can be found by 
    examining its entire definition.
    
    Class attributes
    ----------------
    bb :
    
    bb_heavy_atoms_per_molecule :     
    
    lk :
    
    lk_heavy_atoms_per_molecule :
    
    linked_mols :    
    
    
    Attributes
    ----------
    element : 
    
    number :
    
    heavy_atom_num :
    
    x :
    
    y : 
    
    z :
    
    min_parter :
    
    min_distance :
    
    distances : set
    
    paired : bool

    """
    bb = None
    bb_heavy_atoms_per_molecule = None
    lk = None    
    lk_heavy_atoms_per_molecule = None
    linked_mols = []    
    
    def __init__(self, element, number, heavy_atom_num, x, y, z):
        self.element = element
        self.number = number
        self.heavy_atom_num = heavy_atom_num
        self.x = x
        self.y = y
        self.z = z
        self.min_partner = 10**9
        self.min_distance = 10**9
        self.distances = {}        
        self.paired = False

    @staticmethod
    def extract_mol_file_data(mol_file):
        """
        Reads a ``.mol`` file and stores its data.
        
        This function reads the ``.mol`` containing heavy building
        blocks before they are connected during assembly. The data
        extracted by this function is used to connect the right
        molecules.
        
        Parameters
        ----------
        mol_file : str
            Full path of the ``.mol`` file which holds the heavy
            assembled molecule. This is the file that needs to be edited
            so that bonds are created to form an assembled molecule.
        
        Returns
        -------
        MolFileData
            An object used for storing the content of a ``.mol`` file
            being edited during assembly.
        
        """
    
        # This function goes through the ``.mol`` file holding 
        # disconnected heavy molecules, line by line. The first time a 
        # line describing a heavy atom is found its element is saved as
        # the ``atom1_element``. An ``Atom`` instance representing this 
        # atom is added to ``atom1_list`` and the count, 
        # ``type1_heavy_atom_number``, is incremented.  The next time an 
        # atom of this element is encountered an ``Atom`` object is
        # again created and added to the list and the count incremented.
    
        # After the first heavy atom is encountered, the next heavy atom
        # of a different element is used the initiate the ``atom2`` 
        # variables. The process is otherwise the same for heavy atoms
        # of this element as was the case for the first.
    
        # In the meantime, every line that is iterated through is saved
        # in the string ``mol_file__content``. This content is used by
        # other functions to write the final cage structure to a 
        # ``.mol`` file.
    
        # Once the file is iterated through, the data collected is 
        # returned in a ``MolFileData`` object.
        mol_file_content = ""
        with open(mol_file, "r") as mol_file:
            
            # Store the atomic symbols of the heavy elements in the 
            # ``.mol`` file as a string.
            atom1_element = None    
            atom2_element = None
            
            # Store ``Atom`` instances desrcibing heavy atoms in the 
            # ``.mol`` file.
            atom1_list = []
            atom2_list = []
        
            # This flag is ``True`` when the lines being iterated 
            # through belong to the block describing atoms in the 
            # ``.mol`` file. It is ``False`` otherwise.
            take_atom = False
            
            # This flag is ``True`` when the lines being iterated 
            # through belong to the block describing bonds in the 
            # ``.mol`` file. It is ``False`` otherwise.            
            take_bond = False
            
            # This flag is ``True`` when the lines being iterated
            # will be added to the string holding the ``.mol`` file 
            # content. Some lines do not need to be copied and are added
            # to the final ``.mol`` later. This is because they need to
            # reflect the changes made during the course of the assembly
            # so copying them makes little sense.
            write_line = True
        
            # Count the total of heavy atoms of each type in the 
            # ``.mol`` file.
            type1_heavy_atom_number = 1
            type2_heavy_atom_number = 1
        
            for raw_line in mol_file:
                
                # Split the line into a list of ``words``.
                line = raw_line.split()
                
                # Indicates that the bond block has ended. After this
                # point the lines do not need to be saved.
                if "M  V30 END BOND" in raw_line:
                    write_line = False            
                
                # This string indicates the ``count line`` which holds
                # the total number of atoms in the ``.mol`` file. This
                # number needs to be saved, as the does line itself.
                if "M  V30 COUNTS" in raw_line:
                    count_line = raw_line
                    at_num = line[3]            
                
                # If the ``write_line`` flag is ``True`` copy the line.
                if write_line == True:            
                    mol_file_content += raw_line
                
                # This string indicates that the atom block begins on
                # the next line. Set the appropriate flag to ``True`` so
                # that atomic data starts being collected.
                if "M  V30 BEGIN ATOM" in raw_line:
                    take_atom = True
                    continue
                
                # This string indicates taht the atom block has ended.
                # Set the appropriate flag to ``False`` so that atomic
                # data stops being collected.
                if "M  V30 END ATOM" in raw_line:
                    take_atom = False
                    continue
                
                # This string indicates that the bond block is about to
                # start. Set the appropriate flag to ``True`` so that
                # bond data is collected.
                if "M  V30 BEGIN BOND" in raw_line:
                    take_bond = True
                    continue
                
                if take_atom == True:                 
                    # If the flag was set to true the line describes an
                    # atom. As a result the second word in the line is
                    # the atoms id number, the next is the elemental 
                    # symbol and the next 3 are the x, y and z 
                    # coordinates, respectively.
                    atom_id = line[2]
                    atomic_symbol = line[3]
                    atom_x = float(line[4])
                    atom_y = float(line[5])
                    atom_z = float(line[6])

                    # This flag is ``True`` if the atom described in the
                    # line has the same element as one of the elements
                    # belonging to the functional groups used by MMEA.
                    atom_is_heavy = (atomic_symbol in 
                                                Topology.heavy_symbols)
                    
                    # These flags check to see if a heavy atom element
                    # was previously found.                    
                    atom1_found = atom1_element != None
                    atom2_found = atom2_element != None
                    
                    # If the atom in the line is heavy and is the first
                    # heavy atom of any type found in the ``.mol`` file.
                    if( atom_is_heavy and not atom1_found 
                                      and not atom2_found):
                            
                            # Registers its atomic symbol as the atomic
                            # symbol of the first heavy element.
                            atom1_element = atomic_symbol
                            
                            # Register a ``Atom`` instance for the atom.
                            atom1_list.append(
                                        Atom(atomic_symbol, atom_id, 
                                             type1_heavy_atom_number, 
                                              atom_x, atom_y, atom_z))
                                              
                            # Increment the relevant count.
                            type1_heavy_atom_number += 1                    
                            continue
                    
                    # If the atom in the line is heavy and is the second
                    # type of heavy element to be found - but not
                    # necessarily the second heavy atom to be found.
                    if (atom_is_heavy and atom1_found 
                                    and not atom2_found 
                                    and atomic_symbol != atom1_element):
                            
                            atom2_element = atomic_symbol
                            atom2_list.append(
                                        Atom(atomic_symbol, atom_id, 
                                             type2_heavy_atom_number, 
                                              atom_x, atom_y, atom_z))
                            type2_heavy_atom_number += 1                        
                            continue
                    
                    # If the atom in the line is heavy and is the second
                    # or later atom of the first heavy element found.
                    if (atom_is_heavy and atom1_found and 
                                    atom1_element == atomic_symbol):
                            
                            atom1_list.append(
                                            Atom(atomic_symbol, atom_id, 
                                             type1_heavy_atom_number, 
                                              atom_x, atom_y, atom_z))
                            type1_heavy_atom_number += 1
                            continue
                    # If the atom in the line is heavy and is the second
                    # or later atom of the second heavy element found.
                    if (atom_is_heavy and atom2_found and 
                                    atom2_element == atomic_symbol):
                             
                            atom2_list.append(
                                            Atom(atomic_symbol, atom_id, 
                                               type2_heavy_atom_number, 
                                               atom_x, atom_y, atom_z))
                            type2_heavy_atom_number += 1
                            continue  
                # If this conditions are true the line holds the number
                # of bonds in the ``.mol`` file. Save this number.
                if take_bond == True and len(line) == 6:                   
                    bond_number = int(line[2])
        
        return MolFileData(mol_file_content, count_line, at_num, 
                           bond_number, atom1_list, atom2_list)
        
    def assign_molecule_number(self):
        """
        Provides the id of the molecule in which the atom is found.
        
        Before the building blocks are connected during assembly, the 
        building block molecules are placed together in a ``.mol`` file.
        Each of these molecules will have it own id, to prevent heavy
        atoms in the same mocule being connected during assembly.        
        
        """
        
        # If a building block (note not building-block*) has ``x`` 
        # functional groups,the first ``x`` heavy atom molecules of one 
        # element belong to the first building block molecule, the next 
        # ``x`` heavy atoms belong to the second building block molecule
        # (of that type) and so on. As a result, the id of the molecule
        # the heavy atom belongs to can be found by taking which heavy
        # atom of that element it was and dividing it by the number 
        # of functional groups per building block of that type. The 
        # whole number part of the resulting float is the id of the 
        # molecule.
        
        if self.element == Atom.bb:
            self.mol_number = math.ceil(self.heavy_atom_num / 
                                    Atom.bb_heavy_atoms_per_molecule)  
        if self.element == Atom.lk:
            self.mol_number = math.ceil(self.heavy_atom_num / 
                                    Atom.lk_heavy_atoms_per_molecule)
        
        return 1
    
    def distance(self, atom2):
        """
        Finds the distance between two atoms.
        
        """
        
        x_diff_sq = (self.x - atom2.x) ** 2
        y_diff_sq = (self.y - atom2.y) ** 2
        z_diff_sq = (self.z - atom2.z) ** 2
        r = np.sqrt(x_diff_sq + y_diff_sq + z_diff_sq)
        self.distances[atom2] = r
        atom2.distances[self] = r
        
        return r
        
    def pair_up_nitro(self):
        """
        Finds atoms pairs to join in Nitroso cages.        
        
        This case is specific for the Nitroso cages, for which a double 
        counting for the atoms needs to be avoided.
        
        """
        if self.paired == True:
            return None
        
        if self.paired == False: 
            while True:
                min_partner = min(self.distances, 
                                  key=self.distances.get)
                
                if (min_partner.paired == False and 
                (self.mol_number, min_partner.mol_number) not in 
                Atom.linked_mols and 
                self.mol_number != min_partner.mol_number and 
                ((self.number, self.element, min_partner.number, 
                                              min_partner.element) 
                or (min_partner.number, min_partner.element, 
                    self.number, self.element)) not in 
                    Atom.linked_atoms):

                    self.min_partner = min_partner
                    self.paired = True
                    min_partner.min_partner = self
                    min_partner.paired = True
                    Atom.linked_mols.append((self.mol_number, 
                                            min_partner.mol_number))
                    Atom.linked_atoms.append((self.number, 
                                              self.element, 
                                              min_partner.number, 
                                              min_partner.element))
                    
                    return None
                else:
                    del self.distances[min_partner]

    def pair_up_v3_v2(self):
        if self.paired == True:
            return None
        
        if self.paired == False:
            
#                    print("  DISTANCES: \n  {} {}".format(self.number, self.distances.values()))
            min_partner = min(self.distances, 
                              key=self.distances.get)
#                    print("MIN PARTNER PAIRED? {}".format(min_partner.paired))
            
            if (min_partner.paired == False and 
            (self.mol_number, min_partner.mol_number) not in 
            Atom.linked_mols):
#                    if min_partner.paired == False:
                self.min_partner = min_partner
                self.paired = True
                min_partner.min_partner = self
                min_partner.paired = True
                Atom.linked_mols.append((self.mol_number, 
                                        min_partner.mol_number))
                Atom.linked_atoms.append((self.number, 
                                          self.element, 
                                          min_partner.number, 
                                          min_partner.element))
                
                return None
            else:
                del self.distances[min_partner]
                    
    def pair_up_v4_v2(self):
        if self.paired == True:
            return None
        
        if self.paired == False:
            while True:

                min_partner = min(self.distances, 
                                  key=self.distances.get)
                
                if (min_partner.paired == False and 
                (self.mol_number, min_partner.mol_number) not in 
                Atom.linked_mols):

                    self.min_partner = min_partner
                    self.paired = True
                    min_partner.min_partner = self
                    min_partner.paired = True
                    Atom.linked_mols.append((self.mol_number, 
                                        min_partner.mol_number))
                    Atom.linked_atoms.append((self.number, 
                                              self.element, 
                                             min_partner.number, 
                                           min_partner.element))
                    
                    return None
                else:
                    del self.distances[min_partner]        

    def pair_up_1plus1(self):
        if self.paired == True:
            return None
            
        if self.paired == False:
            while True:

                min_partner = min(self.distances, 
                                  key=self.distances.get)
                
                if (min_partner.paired == False and 
                ((self.number, self.element, min_partner.number, 
                  min_partner.element) 
                or (min_partner.number, min_partner.element, 
                    self.number, self.element)) not in 
                    Atom.linked_atoms):

                    self.min_partner = min_partner
                    self.paired = True
                    min_partner.min_partner = self
                    min_partner.paired = True
                    Atom.linked_mols.append((self.mol_number, 
                                        min_partner.mol_number))
                    Atom.linked_atoms.append((self.number, 
                                              self.element,
                                             min_partner.number, 
                                           min_partner.element))
                    
                    return None
                else:
                    del self.distances[min_partner]

    def pair_up_2plus2(self):
        """
        For 2+2 we have been experiencing ValueError for specific bb and lks.
        """
        
        if self.paired == True:
            return None

        if self.paired == False:
            min_partner = min(self.distances, 
                              key=self.distances.get)                    
            
            bonds_between_mols = (Atom.linked_mols.count(
            (self.mol_number, min_partner.mol_number)) + 
            Atom.linked_mols.count((min_partner.mol_number, 
                                            self.mol_number)))

            #In this topology each molecule cannot bind more than twice to an other molecule
            if (min_partner.paired == False and
                    bonds_between_mols <= 2):
                self.min_partner = min_partner
                self.paired = True
                min_partner.min_partner = self
                min_partner.paired = True
                Atom.linked_mols.append((self.mol_number, 
                                         min_partner.mol_number))
                Atom.linked_atoms.append((self.number, 
                                          self.element, 
                                          min_partner.number, 
                                          min_partner.element))
                
                return None
            else:
                del self.distances[min_partner]

    def pair_up_4plus4(self):
        if self.paired == True:
            return None
        
        if self.paired == False: 
            while True:
                min_partner = min(self.distances, 
                                  key=self.distances.get)
                
                if (min_partner.paired == False and 
                (self.mol_number, min_partner.mol_number) not in 
                Atom.linked_mols):

                    self.min_partner = min_partner
                    self.paired = True
                    min_partner.min_partner = self
                    min_partner.paired = True
                    Atom.linked_mols.append((self.mol_number, 
                                        min_partner.mol_number))
                    Atom.linked_atoms.append((self.number, 
                                              self.element, 
                                             min_partner.number, 
                                        min_partner.element))
                    
                    return None
                else:
                    del self.distances[min_partner]