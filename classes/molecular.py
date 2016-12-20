import numpy as np
from functools import total_ordering, partial
import itertools
from rdkit import Chem as chem
from rdkit.Chem import AllChem as ac
import rdkit.Geometry.rdGeometry as rdkit_geo
from rdkit import DataStructs
import os
import networkx as nx
from scipy.spatial.distance import euclidean
import pickle
import subprocess as sp

from ..convenience_tools import (flatten, periodic_table, 
                                 normalize_vector, rotation_matrix,
                                 vector_theta, LazyAttr,
                                 rotation_matrix_arbitrary_axis,
                                 mol_from_mol_file, 
                                 mol_from_mae_file, ChargedMolError,
                                 FGInfo)

from .exception import MacroMolError

class CachedMacroMol(type):
    """
    A metaclass for creating classes which create cached instances.
    
    This class is tailored to the needs of createding cached
    ``MacroMolecule`` instances.
    
    Extending MMEA
    --------------
    If a MacroMolecule class is added such that one of its initializer
    args of kwargs should not be used for caching, the `__call__()` 
    method in this class will need to be modified.    
    
    """    
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)        
        self._cache = dict()
    
    def __call__(self, *args, **kwargs):
        # Sort the first argument which corresponds to an iterable of
        # StructUnit instances. You want (bb1, lk1) to register as the
        # same as (lk1, bb1). Because this creates the same cage.         
        _, *other_args = args
        args = [sorted(args[0])]              
        args.extend(other_args)
        # Do not take the last arg because this is the name of the
        # ``.mol`` file. This will likely be different even if the cage
        # is the same. However even in this case you want to return the
        # cached instance.
        key = str(args[:-1]) + str(kwargs)
        if key in self._cache.keys():
            return self._cache[key]
        else:
            obj = super().__call__(*args, **kwargs)
            obj.key = key
            self._cache[key] = obj            
            return obj

    def _update_cache(self, macro_mol):
        """
        Updates the cache of stored molecule.

        Parallel processes, such as optimization, return a copy of the
        optimized molecule. The original molecule, stored in `_cache`,
        does not have its attributes updated. In order to replace the 
        molecule with the optimized copy within the cache this function 
        should be used.
        
        Parameters
        ----------
        population : iterable of MacroMolecule instances
            A population holding molecules which should replace the
            ones held in `_cache` which share the same key.
        
        Returns
        -------
        None : NoneType
        
        """
           
        self._cache[macro_mol.key] = macro_mol


class Cached(type):
    """
    A metaclass for creating classes which create cached instances.
    
    """    
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)        
        self.__cache = dict()
    
    def __call__(self, *args, **kwargs):
        key = str(args) + str(kwargs)
        if key in self.__cache.keys():
            return self.__cache[key]
        else:
            obj = super().__call__(*args, **kwargs)
            self.__cache[key] = obj
            return obj

class Energy:
    """
    Handles all things related to a ``Molecule``'s energy.

    An instance of this class will be placed in the `energy` attribute
    of a ``Molecule`` instance.    

    Attributes
    ----------
    molecule : Molecule    
        The energetic information held by an instance of ``Energy`` 
        concerns the molecule held in this attribute.
    
    values : dict
        The keys in the dict code for the function and parameters which
        were used to calculate the energy. The values are the energies.

    """

    def __init__(self, molecule):
        self.molecule = molecule
        self.values = {}
        
    def formation(self, key, products, 
                  building_blocks=None, force_e_calc=False):
        """
        Calculates the formation energy.
        
        The formation energy is calculated under the assumption that the
        molecule in `self.molecule` is composed of the molecules in 
        `building_blocks` and that during formation molecules in 
        `products` are formed in addition to `self.molecule`.     
        
        Parameters
        ----------
        key : tuple
            The first member of the tuple is a string holding the name 
            of a method used to calculate energies. For exmaple 'rdkit'
            or 'macromodel'. The remaning elements in the tuple are the
            parameters that the user wishes to pass to the function.
            
        products : tuple of form (float, Molecule)
            This tuple holds the molecules produced in addition to 
            `self.molecule`, when  a single `self.molecule` is made. The
            ``int`` represents the number made per `self.molecule`.

        building_blocks : tuple (default = None)
            This argument should be a tuple of the form 
            (float, Molecule). It holds the number of a given Molecule
            required to build a single molecule held in `self.molecule`.
            This argument can be omitted when the formation energy of a 
            MacroMolecule instance is being found, as they keep this 
            data stored elsewhere already.

        force_e_calc : bool (default = False)
            If the this is ``True`` then all building blocks, products
            and `self.molecule` will have their energies recalculated.
            Even if the energy values have already been found if the 
            chosen forcefield/method. If ``False`` the energy is only
            calculated if the value has not already been foud.
        
        Modifies
        --------
        self.values : dict
            Adds an entry to this dictionary. The key is a tuple of the
            form ('formation', key[0], key[1]). The value is the 
            calculated formation energy.            
        
        Returns
        -------
        float
            The formation energy. Note that this value is also stored
            in the dictionary `self.values`.
        
        """

        func_name, *params = key
        
        # Recalculate energies if requested.
        if force_e_calc:
            for _, mol in products:
                getattr(mol.energy, func_name)(*params)
        
        e_products = 0
        for n, mol in products:
            if (func_name, params[0]) not in mol.energy.values.keys():
                getattr(mol.energy, func_name)(*params)
            
            e_products += n * mol.energy.values[(func_name, params[0])]
        
        eng = self.pseudoformation(key, building_blocks, force_e_calc) 
        eng -= e_products
        print(e_products)
        print(products[0][1].energy.values)
        self.values[('formation', func_name, params[0])] = eng         
        return eng

    def pseudoformation(self, key, 
                        building_blocks=None, force_e_calc=False):
        """
        Calculates the formation energy, sans other products.

        This is the formation energy if the energy of the other products
        of the reaction is not taken into account.
        
        Parameters
        ----------
        key : tuple
            The first member of the tuple is a string holding the name 
            of a method used to calculate energies. For exmaple 'rdkit'
            or 'macromodel'. The remaning elements in the tuple are the
            parameters that the user wishes to pass to the function.

        building_blocks : tuple (default = None)
            This argument should be a tuple of the form 
            (float, Molecule). It holds the number of a given Molecule
            required to build a single molecule held in `self.molecule`.
            This argument can be omitted when the formation energy of a 
            MacroMolecule instance is being found, as they keep this 
            data stored elsewhere already.

        force_e_calc : bool (default = False)
            If the this is ``True`` then all building blocks, products
            and `self.molecule` will have their energies recalculated.
            Even if the energy values have already been found if the 
            chosen forcefield/method. If ``False`` the energy is only
            calculated if the value has not already been foud.        
        
        Modifies
        --------
        self.values : dict
            Adds an entry to this dictionary. The key is a tuple of the
            form ('formation', key[0], key[1]). The value is the 
            calculated formation energy sans products (pseudoformation).            
        
        Returns
        -------
        float
            The pseudoformation energy. Note that this value is also 
            stored in the dictionary `self.values`.
        
        """
        
        if building_blocks is None:
            building_blocks = ((n, mol) for mol, n in 
                              self.molecule.topology.bb_counter.items())
        
        func_name, *params = key
        
        # Recalculate energies if requested.
        if force_e_calc:
            for _, mol in building_blocks:
                getattr(mol.energy, func_name)(*params)
        
            getattr(self, func_name)(*params)

        # Calculate the energy of building blocks and products using the
        # chosen force field, if it has not been found already.
        e_reactants = 0
        for n, mol in building_blocks:
            if (func_name, params[0]) not in mol.energy.values.keys():
                getattr(mol.energy, func_name)(*params)
            
            e_reactants += n * mol.energy.values[(func_name, params[0])]
        
        e_products = (self.values[(func_name, params[0])] if 
                    (func_name, params[0]) in self.values.keys() else
                    getattr(self, func_name)(*params))

        eng = e_reactants - e_products
        self.values[('pseudoformation', func_name, params[0])] = eng         
        return eng        

    def rdkit(self, forcefield):
        """
        Uses rdkit to calculate the energy of `self.molecule`.
        
        Parameters
        ----------
        forcefield : str
            The name of the forcefield to be used.
            
        Modifies
        --------
        self.values : dict
            Adds an entry to this dictionary. The key is a tuple of the
            form ('rdkit', `forcefield`). The value is the caculated 
            energy.
            
        Returns 
        -------
        float
            The calculated energy. Note that this value is also stored
            in the dictionary `self.values`.
        
        """
        
        if forcefield == 'uff':
            self.molecule.prist_mol.UpdatePropertyCache()
            ff = ac.UFFGetMoleculeForceField(self.molecule.prist_mol)
        if forcefield == 'mmff':
            chem.GetSSSR(self.molecule.prist_mol)      
            ff = ac.MMFFGetMoleculeForceField(self.molecule.prist_mol,
                  ac.MMFFGetMoleculeProperties(self.molecule.prist_mol))

        eng = ff.CalcEnergy()        
        self.values[('rdkit', forcefield)] = eng
        return eng
        
    def macromodel(self, forcefield, macromodel_path):
        """
        Calculates the energy of `self.molecule` using macromodel.

        Note that this requires macromodel to be installed and have a 
        valid license.

        Parameters
        ----------
        forcefield : int
            The id number of the forcefield to be used by macromodel.  
            
        macromodel_path : str
            The full path of the ``Schrodinger`` suite within the user's 
            machine. For example, in a default Microsoft installation 
            the folder will probably be something like
            ``C:\Program Files\Schrodinger2016-2``.
 
        Modifies
        --------
        self.values : dict
            Adds an entry to this dictionary. The key is a tuple of the
            form ('macromodel', `forcefield`). The value is the 
            caculated energy.
            
        Returns 
        -------
        float
            The calculated energy. Note that this value is also stored
            in the dictionary `self.values`.
        
        """
        
        file_root, ext = os.path.splitext(self.molecule.prist_mol_file)

        # Create a .mae file of the molecule if it does not exist
        # already.
        if ext != '.mae':
            convrt_app = os.path.join(macromodel_path, 'utilities', 
                                                        'structconvert')
            convrt_cmd = [convrt_app, 
                         self.molecule.prist_mol_file, file_root+'.mae']
            sp.call(convrt_cmd)

        # Create an input file and run it.        
        input_script = (
        "{0}.mae\n"
        "{0}-out.maegz\n"
        " MMOD       0      1      0      0     0.0000     0.0000     "
        "0.0000     0.0000\n"
        " FFLD{1:8}      1      0      0     1.0000     0.0000     "
        "0.0000     0.0000\n"
        " BGIN       0      0      0      0     0.0000     0.0000     "
        "0.0000     0.0000\n"
        " READ      -1      0      0      0     0.0000     0.0000     "
        "0.0000     0.0000\n"
        " ELST      -1      0      0      0     0.0000     0.0000     "
        "0.0000     0.0000\n"
        " WRIT       0      0      0      0     0.0000     0.0000     "
        "0.0000     0.0000\n"
        " END       0      0      0      0     0.0000     0.0000     "
        "0.0000     0.0000\n\n"
        ).format(file_root, forcefield)
        
        with open(file_root+'.com', 'w') as f:
            f.write(input_script)
        
        cmd = [os.path.join(macromodel_path,'bmin'), 
               file_root, "-WAIT", "-LOCAL"]
        sp.call(cmd)
        
        # Read the .log file and return the energy.
        with open(file_root+'.log', 'r') as f:
            for line in f:
                if "                   Total Energy =" in line:
                    eng = float(line.split()[-2].replace("=", ""))
    
        self.values[('macromodel', forcefield)] = eng
        return eng
        
class Molecule:
    """
    The most basic class representing molecules.

    This class defines defines the operations which any class describing
    molecules should inherit or may find useful. Examples of such 
    classes are ``StructUnit`` and ``MacroMolecule``. This calls is
    unlikely to be useful as in and of itself. It lacks an `__init__()`
    method because it does not need one. Child classes should define it.
    
    It is assumed that child classes will have some basic attributes.
    
    Minimum required attributes of child classes
    --------------------------------------------
    prist_mol : rdkit.Chem.rdchem.Mol
        A rdkit molecule instance representing the molecule.
        
    heavy_mol : rdkit.Chem.rdchem.Mol
        A rdkit molecule instance representing the molecule where the
        functional groups have been substited with heavy atoms.
        
    prist_mol_file : str
        The full path of the molecular structure file holding the 
        structure of the pristine molecule.
        
    heavy_mol_file : str
        The full path of the molecular structure file holding the
        structure of the substituted molecule.
        
    energy : Energy
        An instance of the ``Energy`` class. It handles all things
        energy.
        
    
    """

    def all_atom_coords(self, mol_type):
        """
        Yields the coordinates of atoms in `heavy_mol` or `prist_mol`.        

        This yields the coordinates of every atom in the pristine or
        heavy molecule, depending on `mol_type`.

        The `heavy_mol` and `prist_mol` attributes hold a 
        ``rdkit.Chem.rdchem.Mol`` instance. This instance holds a 
        ``rdkit.Chem.rdchem.Conformer`` instance. The conformer instance
        holds the positions of the atoms within that conformer. This
        generator yields those coordinates.

        Parameters
        ----------
        mol_type : str (allowed values = 'heavy' or 'prist')
            A string which defines whether the pristine or heavy
            molecule is used.   
       
        Yields
        ------
        tuple of (int, numpy.array)
            The ``int`` represents the atom id of the atom whose
            coordinates are being yielded.
        
            The array represents the complete position in space. Each 
            float within the array represents the value of the x, y or z 
            coordinate of an atom. The x, y and z coordinates are 
            located in the array in that order. 
            
        Raises
        ------
        ValueError
            If the `mol_type` string is not 'heavy' or 'prist'.
        
        """
        
        if mol_type not in {'prist', 'heavy'}:
            raise ValueError(("`mol_type` must be either 'prist'"
                                " or 'heavy'."))
        
        if mol_type == 'prist':
            rdkit_mol = self.prist_mol
            
        elif mol_type == 'heavy':
            rdkit_mol = self.heavy_mol
        
        # Get the conformer from the rdkit instance. 
        conformer = rdkit_mol.GetConformer()
        
        # Go through all the atoms and ask the conformer to return
        # the position of each atom. This is done by supplying the 
        # conformers `GetAtomPosition` method with the atom's id.
        for atom in rdkit_mol.GetAtoms():
            atom_id = atom.GetIdx()
            atom_position = conformer.GetAtomPosition(atom_id)
            yield atom_id, np.array([*atom_position]) 

    def all_heavy_atom_distances(self):
        """
        Yield distances between all pairs of heavy atoms in `heavy_mol`.
        
        All distances are only yielded once. This means that if the 
        distance between atoms with ids ``1`` and ``2``is yielded as
        ``(12.4, 1, 2)``, no tuple of the form ``(12.4, 2, 1)`` will be 
        yielded.
        
        Only distances between heavy atoms used for functional group
        substitutions are considered. Distances between heavy atoms
        and regular atoms or between two regular atoms are not yielded.
        
        Yields
        ------
        tuple of form (scipy.double, int, int)
            This tuple holds the distance between two heavy atoms. The 
            first element is the distance and the next two are the 
            relevant atom ids.

        """
                
        # Iterate through each pair of atoms - do not allow iterating
        # through the same pair twice. In otherwords, recombinations of 
        # the same pair are not allowed.
        for atom1, atom2 in itertools.combinations(
                                        self.heavy_mol.GetAtoms(), 2):
 
            # Only yield if both atoms are heavy. 
            if (atom1.GetAtomicNum() in FGInfo.heavy_atomic_nums and 
                atom2.GetAtomicNum() in FGInfo.heavy_atomic_nums):               
                
                # Get the atom ids, use them to calculate the distance
                # and yield the resulting data.
                atom1_id = atom1.GetIdx()
                atom2_id = atom2.GetIdx()
                yield (self.atom_distance('heavy', atom1_id, atom2_id), 
                      atom1_id, atom2_id)

    def atom_coords(self, mol_type, atom_id):
        """
        Return coordinates of an atom in `heavy_mol` or `prist_mol`.

        Parameters
        ----------
        mol_type : str (allowed values = 'heavy' or 'prist')
            A string which defines whether the pristine or heavy
            molecule is used.             
        
        atom_id : int
            The id of the atom whose coordinates are desired.
            
        Returns
        -------
        numpy.array
            An array hodling the x, y and z coordinates (respectively) 
            of atom with id of `atom_id` in `heavy_mol`.
        
        Raises
        ------
        ValueError
            If the `mol_type` string is not 'heavy' or 'prist'.
        
        """
        
        if mol_type not in {'prist', 'heavy'}:
            raise ValueError(("`mol_type` must be either 'prist'"
                                " or 'heavy'."))
        
        if mol_type == 'prist':
            rdkit_mol = self.prist_mol
            
        elif mol_type == 'heavy':
            rdkit_mol = self.heavy_mol
        
        conf = rdkit_mol.GetConformer()
        atom_position = conf.GetAtomPosition(atom_id)
        return np.array([*atom_position])

    def atom_distance(self, mol_type, atom1_id, atom2_id):
        """
        Return the distance between atoms in `heavy_mol` or `prist_mol`.
        
        Parameters
        ----------
        mol_type : str (allowed values = 'heavy' or 'prist')
            A string which defines whether the pristine or heavy
            molecule is used.
            
        atom1_id : int
            The id of the first atom.
        
        atom2_id : int
            The id of the second atom.
            
        Returns 
        -------
        scipy.double
            The distance between the first and second atoms.

        Raises
        ------
        ValueError
            If the `mol_type` string is not 'heavy' or 'prist'.
        
        """
        
        if mol_type not in {'prist', 'heavy'}:
            raise ValueError(("`mol_type` must be either 'prist'"
                                " or 'heavy'."))

        # Get the atomic positions of each atom and use the scipy 
        # function to calculate their distance in Euclidean space.              
        atom1_coords = self.atom_coords(mol_type, atom1_id)
        atom2_coords = self.atom_coords(mol_type, atom2_id)
        
        return euclidean(atom1_coords, atom2_coords)

    def atom_symbol(self, mol_type, atom_id):
        """
        Returns the symbol of the atom with id `atom_id`.
        
        Parameters
        ----------
        mol_type : str (allowed values = 'heavy' or 'prist')
            A string which defines whether the pristine or heavy
            molecule is used.        

        Returns
        -------
        str
            The atomic symbol of the atom with id number `atom_id`.        
        
        Raises
        ------
        ValueError
            If the `mol_type` string is not 'heavy' or 'prist'.
        
        """
        
        if mol_type not in {'prist', 'heavy'}:
            raise ValueError(("`mol_type` must be either 'prist'"
                                " or 'heavy'."))        
        if mol_type == 'prist':
            rdkit_mol = self.prist_mol            
        elif mol_type == 'heavy':
            rdkit_mol = self.heavy_mol
            
        atom = rdkit_mol.GetAtomWithIdx(atom_id)
        atomic_num = atom.GetAtomicNum()
        return periodic_table[atomic_num] 

    def center_of_mass(self, mol_type):
        """
        Returns the centre of mass of the pristine or heavy molecule.

        Parameters
        ----------
        mol_type : str (allowed values = 'heavy' or 'prist')
            A string which defines whether the pristine or heavy
            molecule is used.
            
        Returns
        -------
        numpy.array
            The array holds the x, y and z coordinates of the center of
            mass, in that order.
                
        Raises
        ------
        ValueError
            If the `mol_type` string is not 'heavy' or 'prist'.

        References
        ----------
        https://en.wikipedia.org/wiki/Center_of_mass
        
        """
        
        if mol_type not in {'prist', 'heavy'}:
            raise ValueError(("`mol_type` must be either 'prist'"
                                " or 'heavy'."))        
        if mol_type == 'prist':
            rdkit_mol = self.prist_mol            
        elif mol_type == 'heavy':
            rdkit_mol = self.heavy_mol
              
        center = np.array([0,0,0])        
        total_mass = 0
        for atom_id, coord in self.all_atom_coords(mol_type):
            mass = rdkit_mol.GetAtomWithIdx(atom_id).GetMass()
            total_mass += mass
            center = np.add(center, np.multiply(mass, coord))

        return np.divide(center, total_mass)

    def centroid(self, mol_type):
        """
        Returns the centroid of the heavy or pristine rdkit molecule.

        Parameters
        ----------
        mol_type : str (allowed values = 'heavy' or 'prist')
            A string which defines whether the pristine or heavy
            molecule is used.

        Returns
        -------
        numpy.array
            A numpy array holding the position of the centroid.
        
        Raises
        ------
        ValueError
            If the `mol_type` string is not 'heavy' or 'prist'.
        
        """
        
        if mol_type not in {'prist', 'heavy'}:
            raise ValueError(("`mol_type` must be either 'prist'"
                                " or 'heavy'."))
        
        if mol_type == 'prist':
            rdkit_mol = self.prist_mol
            
        elif mol_type == 'heavy':
            rdkit_mol = self.heavy_mol
            
        centroid = sum(x for _, x in self.all_atom_coords(mol_type)) 
        return np.divide(centroid, rdkit_mol.GetNumAtoms())

    def dump(self, file_name):
        """
        Write the MacroMolecule object to a file.
        
        Parameters
        ----------
        file_name : str
            The full path of the file to which the macromolecule should
            be written.
            
        Returns
        -------
        None : NoneType

        """
        
        with open(file_name, 'wb') as dump_file:    
            pickle.dump(self, dump_file)     

    def graph(self, mol_type):
        """
        Returns a mathematical graph representing the molecule.        
        
        Parameters
        ----------        
        mol_type : str (allowed values = 'heavy' or 'prist')
            A string which defines whether the pristine or heavy
            molecule is used. 
       
        Returns
        -------
        networkx.Graph
            A graph where the nodes are the ids of the atoms in the
            rdkit molecule and the edges are the bonds.
        
        Raises
        ------
        ValueError
            If the `mol_type` string is not 'heavy' or 'prist'.
        
        """
        
        if mol_type not in {'prist', 'heavy'}:
            raise ValueError(("`mol_type` must be either 'prist'"
                                " or 'heavy'."))        
        if mol_type == 'prist':
            rdkit_mol = self.prist_mol            
        elif mol_type == 'heavy':
            rdkit_mol = self.heavy_mol
        
        # Create a graph instance and add the atom ids as nodes. Use the
        # the atom ids from each end of a bond to define edges. Do this
        # for all bonds to account for all edges.
        
        graph = nx.Graph()        
        
        for atom in rdkit_mol.GetAtoms():
            graph.add_node(atom.GetIdx())
        
        for bond in rdkit_mol.GetBonds():
            graph.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
             
        return graph

    def heavy_atom_position_matrix(self):
        """
        Returns a matrix holding positions of heavy atoms.

        Returns
        -------
        numpy.matrix
            The matrix is 3 x n. Each column holds the x, y and z
            coordinates of a heavy atom. The index of the column 
            corresponds to the index of the atom id in `heavy_ids`.    
        
        """
        
        pos_array = np.array([])

        for atom_id in self.heavy_ids:
            pos_vect = np.array([*self.atom_coords('heavy', atom_id)])
            pos_array = np.append(pos_array, pos_vect)

        return np.matrix(pos_array.reshape(-1,3).T)

    def heavy_direction_vector_atoms(self):
        """
        Yields the atoms which form `heavy_direction_vectors`.

        The method `heavy_direction_vectors` yields direction vectors.
        This method yields the atom ids of atoms which form those 
        direction vectors. The atoms ids are yielded in the same order
        as the vectors.

        Yield
        -------
        tuple of ints
            The tuple holds the start and end atom ids of a
            `heavy_direction_vector`. The first tuple will correspond 
            to the atoms forming the first `heavy_direction_vector` and
            so on. The atom id at index 0 is the start atom and the atom
            id at index 1 is the end atom.
        
        """

        for atom1_id, atom2_id in itertools.combinations(self.heavy_ids,
                                                                     2):
            yield atom2_id, atom1_id  

    def heavy_direction_vectors(self):
        """
        Yields the direction vectors between all pairs of heavy atoms.
                
        The yielded vector is normalized. If a pair (1,2) is yielded, 
        the pair (2,1) will not be.
        
        Yields
        ------
        The next direction vector running between a pair of heavy atoms.
        
        """
        
        for atom1_id, atom2_id in itertools.combinations(self.heavy_ids, 
                                                                     2):
            p1 = self.atom_coords('heavy', atom1_id)
            p2 = self.atom_coords('heavy', atom2_id)
        
            yield normalize_vector(p1-p2)

    @staticmethod
    def load(file_name):
        """
        Initializes a MacroMolecule from one dumped to a file.
        
        Parameters
        ----------
        file_name : str
            The full path of the file holding the dumped molecule.
            
        Returns
        -------
        MacroMolecule
            The macromolecule stored in the dump file.
            
        """
        
        with open(file_name, 'rb') as dump_file:
            return pickle.load(dump_file)

    def position_matrix(self, mol_type):
        """
        Return position of all atoms in the building block as a matrix. 
        
        Positions of the atoms in the heavy or pristine rdkit instance
        will be returned based on `mol_type`.

        Parameters
        ----------
        mol_type : str (allowed values = 'heavy' or 'prist')
            A string which defines whether the pristine or heavy
            molecule is used.           
        
        Returns
        -------
        numpy.matrix
            The matrix is 3 x n. Each column holds the x, y and z
            coordinates of an atom. The index of the column corresponds
            to the id of the atom in the rdkit molecule whose 
            coordinates it holds.
        
        Raises
        ------
        ValueError
            If the `mol_type` string is not 'heavy' or 'prist'.
        
        """
        
        if mol_type not in {'prist', 'heavy'}:
            raise ValueError(("`mol_type` must be either 'prist'"
                                " or 'heavy'."))       
        if mol_type == 'prist':
            rdkit_mol = self.prist_mol            
        elif mol_type == 'heavy':
            rdkit_mol = self.heavy_mol

        pos_array = np.array([])
        for atom in rdkit_mol.GetAtoms():
            atom_id = atom.GetIdx()
            pos_vect = np.array([*self.atom_coords(mol_type, atom_id)])
            pos_array = np.append(pos_array, pos_vect)

        return np.matrix(pos_array.reshape(-1,3).T)

    def rotate(self, mol_type, theta, axis):
        """
        Rotates the heavy rdkit molecule by `theta` about `axis`.
        
        The rotation occurs about the molecular centroid.        
        
        Parameters
        ----------    
        mol_type : str (allowed values = 'heavy' or 'prist')
            A string which defines whether the pristine or heavy
            molecule is used.           
        
        theta : float
            The size of the rotation in radians.
        
        axis : numpy.array
            The axis about which rotation happens.
        
        Modifies
        --------
        prist_mol : rdkit.Chem.rdchem.Mol
            The atoms in this molecule are rotated if `mol_type` is 
            'prist'.
        
        heavy_mol : rdkit.Chem.rdchem.Mol
            The atoms in this molecule are rotated if `mol_type` is 
            'heavy'.
    
        Returns
        -------
        None : NoneType
            
        """

        if mol_type not in {'prist', 'heavy'}:
            raise ValueError(("`mol_type` must be either 'prist'"
                                " or 'heavy'."))       
        
        og_position = self.centroid(mol_type)
        self.set_position(mol_type, [0,0,0])
        rot_mat = rotation_matrix_arbitrary_axis(theta, axis)
        new_pos_mat = np.dot(rot_mat, self.position_matrix(mol_type))
        self.set_position_from_matrix(mol_type, new_pos_mat)
        self.set_position(mol_type, og_position)
        
    def rotate2(self, theta, axis):
        """
        Rotates the heavy rdkit molecule by `theta` about `axis`.
        
        The rotation occurs about the heavy atom centroid.        
        
        Parameters
        ----------        
        theta : float
            The size of the rotation in radians.
        
        axis : numpy.array
            The axis about which rotation happens.
        
        Modifies
        --------
        heavy_mol : rdkit.Chem.rdchem.Mol
            The atoms in this molecule are rotated.
    
        Returns
        -------
        None : NoneType
            
        """
        
        og_position = self.heavy_atom_centroid()
        self.set_heavy_atom_centroid([0,0,0])
        rot_mat = rotation_matrix_arbitrary_axis(theta, axis)
        new_pos_mat = np.dot(rot_mat, self.position_matrix('heavy'))
        self.set_position_from_matrix('heavy', new_pos_mat)
        self.set_heavy_atom_centroid(og_position)

    def set_orientation(self, mol_type, start, end):
        """
        Rotates heavy molecule by rotation of `start` to `end`.
 
        Note: The difference between this method and 
        `_set_heavy_mol_orientation()` is about which point the rotation
        occurs (centroid of entire molecule versus centroid of heavy 
        atoms, respectively). This method works on both heavy and 
        pristine molecules while `_set_heavy_mol_orientation()` doesn't.
        
        Given two direction vectors, `start` and `end`, this method
        applies the rotation required transform `start` to `end` on 
        the molecule. The rotation occurs about the centroid of the
        molecule.
        
        For example, if the `start` and `end` vectors
        are 45 degrees apart, a 45 degree rotation will be applied to
        the molecule. The rotation will be along the appropriate axis.
        
        The great thing about this method is that you as long as you can 
        associate a gemotric feature of the molecule with a vector, then 
        the molecule can be roatated so that this vector is aligned with 
        `end`. The defined vector can be virtually anything. This means 
        that any geomteric feature of the molecule can be easily aligned 
        with any arbitrary axis.
        
        Parameters
        ----------
        mol_type : str
            Must be either 'prist' or 'heavy' to signify which rdkit
            molecule this method should operate on.
        
        start : numpy.array
            A vector which is to be rotated so that it transforms to the
            `end` vector.
        
        end : numpy.array
            This array holds the directional vector along which the 
            heavy atoms in the linker should be placed.
            
        Modifies
        --------
        prist_mol : rdkit.Chem.rdchem.Mol   
            When `mol_type` is 'prist', the conformer in this rdkit 
            instance is changed due to rotation of the molecule about
            its centroid.
        
        heavy_mol : rdkit.Chem.rdchem.Mol   
            When `mol_type` is 'heavy', the conformer in this rdkit 
            instance is changed due to rotation of the molecule about
            its centroid.
        
        Returns
        -------
        rdkit.Chem.rdchem.Mol
            An rdkit molecule instance of the rotated molecule. This is 
            a copy of the rdkit molecule in `heavy_mol` or `prist_mol`.
            
        Raises
        ------
        ValueError
            If the `mol_type` string is not 'heavy' or 'prist'.
            
        """
        
        if mol_type not in {'prist', 'heavy'}:
            raise ValueError(("`mol_type` must be either 'prist'"
                                " or 'heavy'."))
        
        if mol_type == 'prist':
            rdkit_mol = self.prist_mol
            
        elif mol_type == 'heavy':
            rdkit_mol = self.heavy_mol

        # Normalize the input direction vectors.
        start = normalize_vector(start)
        end = normalize_vector(end)
        
        # Record the position of the molecule then translate the heavy
        # atom centroid to the origin. This is so that the rotation
        # occurs about this point.
        og_center = self.centroid(mol_type)
        self.set_position(mol_type, [0,0,0]) 
        
        # Get the rotation matrix.
        rot_mat = rotation_matrix(start, end)
        
        # Apply the rotation matrix to the atomic positions to yield the
        # new atomic positions.
        new_pos_mat = np.dot(rot_mat, self.position_matrix(mol_type))

        # Set the positions in the heavy rdkit molecule.
        self.set_position_from_matrix(mol_type, new_pos_mat)
        self.set_position(mol_type, og_center)

        return chem.Mol(rdkit_mol)        

            
    def set_position(self, mol_type, position):
        """
        Changes the position of the rdkit molecule. 
        
        This method translates the heavy or pristine building block 
        molecule so that its centroid is on the point `position`. Which
        one is set depends on `mol_type`.
        
        Unlike `shift_mol` this method does change the original
        rdkit molecule found in `prist_mol` or `heav_mol`. A copy is NOT 
        created. This method does not change the ``.mol`` files however.
        
        Parameters
        ----------
        mol_type : str (allowed values = 'heavy' or 'prist')
            A string which defines whether the pristine or heavy
            molecule is used.        
        
        position : numpy.array
            This array holds the position on which the centroid of the 
            building block should be placed.
            
        Modifies
        --------
        prist_mol : rdkit.Chem.rdchem.Mol       
            If `mol_type` is 'prist',  the conformer in this rdkit 
            instance is changed so that its centroid falls on 
            `position`.         
        
        heavy_mol : rdkit.Chem.rdchem.Mol   
            If `mol_type` is 'heavy', the conformer in this rdkit 
            instance is changed so that its centroid falls on 
            `position`.
        
        Returns
        -------
        rdkit.Chem.rdchem.Mol
            The heavy rdkit molecule with the centroid placed at 
            `position`. This is the same instance as that in 
            `heavy_mol`.
        
        Raises
        ------
        ValueError
            If the `mol_type` string is not 'heavy' or 'prist'.
        
        """
        
        if mol_type not in {'prist', 'heavy'}:
            raise ValueError(("`mol_type` must be either 'prist'"
                                " or 'heavy'."))
        
        if mol_type == 'prist':
            rdkit_mol = self.prist_mol
            
        elif mol_type == 'heavy':
            rdkit_mol = self.heavy_mol
        
        centroid = self.centroid(mol_type)
        shift = position - centroid
        new_conf = self.shift(mol_type, shift).GetConformer()

        rdkit_mol.RemoveAllConformers()
        rdkit_mol.AddConformer(new_conf)
        
        return rdkit_mol        

    def set_position_from_matrix(self, mol_type, pos_mat):
        """
        Set atomic positions of rdkit molecule to those in `pos_mat`.
        
        The form of the `pos_mat` is 3 x n. The column of `pos_mat` 
        represents the x, y and z coordinates to which an atom
        in the heavy molecule should be set. The 1st column sets the
        coordinate of the atom in `heavy_mol` with id 0. The next column
        sets the coordinate of the atom in `heavy_mol` with id 1, and
        so on.
        
        Parameters
        ----------
        mol_type : str (allowed values = 'heavy' or 'prist')
            A string which defines whether the pristine or heavy
            molecule is used.       
        
        pos_mat : numpy.array
            The matrix holds the coordinates to which atoms of the heavy 
            molecule should be set. The index of the column in
            the matrix corresponds to the id of the atom who's
            coordinates it holds.            
            
        Modifies
        --------
        heavy_mol : rdkit.Chem.rdchem.Mol
            If `mol_type` is 'heavy', the coordinates of atoms in this 
            molecule are set to the coordinates in `pos_mat`.
    
        prist_mol : rdkit.Chem.rdchem.Mol   
            if `mol_type` is 'prist', the coordinates of atoms in this 
            molecule are set to the coordinates in `pos_mat`.
    
        Returns
        -------
        None : NoneType
        
        Raises
        ------
        ValueError
            If the `mol_type` string is not 'heavy' or 'prist'.
        
        """
        
        if mol_type not in {'prist', 'heavy'}:
            raise ValueError(("`mol_type` must be either 'prist'"
                                " or 'heavy'."))
        
        if mol_type == 'prist':
            rdkit_mol = self.prist_mol
            
        elif mol_type == 'heavy':
            rdkit_mol = self.heavy_mol
        
        conf = rdkit_mol.GetConformer()
        for i, coord_mat in enumerate(pos_mat.T):
            coord = rdkit_geo.Point3D(coord_mat.item(0), 
                                      coord_mat.item(1), 
                                      coord_mat.item(2))
            conf.SetAtomPosition(i, coord)
        
    def shift(self, mol_type, shift):
        """
        Shifts the coordinates of all atoms.
        
        This method will shift the rdkit instance of either the pristine 
        or heavy molecule depending on `mol_type`.
        
        The `heavy_mol` and `prist_mol` attributes hold a 
        ``rdkit.Chem.rdchem.Mol`` instance. This instance holds holds a 
        ``rdkit.Chem.rdchem.Conformer`` instance. The conformer instance
        holds the positions of the atoms in the molecule. This function 
        creates a new conformer with all the coordinates shifted by the
        values supplied in `shift`. This function does not change the 
        existing conformer.
        
        To be clear, consider the following code:
        
            >>> b = a.shift('heavy', [10,10,10])
            >>> c = a.shift('heavy', [10,10,10])
        
        In the preceeding code where ``a`` is a ``StructUnit`` instance, 
        ``b`` and ``c`` are two new ``rdkit.Chem.rdchem.Mol`` instances. 
        The ``rdkit.Chem.rdchem.Mol`` instances held by ``a`` in 
        `prist_mol` and `heavy_mol` are completely unchanged. As are any 
        other attributes of ``a``. Both ``b`` and ``c`` are rdkit 
        molecule instances with conformers which are exactly the same.
        Both of these conformers are exactly like the conformer of the 
        heavy rdkit molecule in ``a`` except all the atomic positions
        are increased by 10 in the x, y and z directions. 
        
        Because ``a`` was not modified by runnig the method, running it
        again with the same arguments leads to the same result. This is 
        why the conformers in ``b`` and ``c`` are the same.

        Parameters
        ----------
        mol_type : str (allowed values = 'heavy' or 'prist')
            A string which defines whether the pristine or heavy
            molecule is used.
            
        shift : numpy.array
            A numpy array holding the value of the shift along the
            x, y and z axes in that order.
        
        Returns
        -------
        rdkit.Chem.rdchem.Mol
            An rdkit molecule instance which has a modified version of 
            the conformer found in either `heavy_mol` or `prist_mol`.
        
        Raises
        ------
        ValueError
            If the `mol_type` string is not 'heavy' or 'prist'.
        
        """
        
        if mol_type not in {'prist', 'heavy'}:
            raise ValueError(("`mol_type` must be either 'prist'"
                                " or 'heavy'."))
        
        if mol_type == 'prist':
            rdkit_mol = self.prist_mol
            
        elif mol_type == 'heavy':
            rdkit_mol = self.heavy_mol
        
        # The function does not modify the existing conformer, as a 
        # result a new instance is created and used for modification.
        conformer = chem.Conformer(rdkit_mol.GetConformer())
        
        # For each atom, get the atomic positions from the conformer 
        # and shift them. Create a new geometry instance from these new
        # coordinate values. The geometry instance is used by rdkit to
        # store the coordinates of atoms. Finally, set the conformers
        # atomic position to the values stored in this newly generated
        # geometry instance.
        for atom in rdkit_mol.GetAtoms():
            
            # Remember the id of the atom you are currently using. It 
            # is used to change the position of the correct atom at the
            # end of the loop.
            atom_id = atom.GetIdx()
            
            # `atom_position` in an instance holding in the x, y and z 
            # coordinates of an atom in its 'x', 'y' and 'z' attributes.
            atom_position = np.array(conformer.GetAtomPosition(atom_id))
            
            # Inducing the shift.
            new_atom_position = atom_position + shift
            
            # Creating a new geometry instance.
            new_coords = rdkit_geo.Point3D(*new_atom_position)            
            
            # Changes the position of the atom in the conformer to the
            # values stored in the new geometry instance.
            conformer.SetAtomPosition(atom_id, new_coords)
        
        # Create a new copy of the rdkit molecule instance representing
        # the molecule - the original instance is not to be modified.
        new_mol = chem.Mol(rdkit_mol)
        
        # The new rdkit molecule was copied from the one held in the
        # `heavy_mol` or `prist_mol` attribute, as result it has a copy 
        # of its conformer. To prevent the rdkit molecule from holding 
        # multiple conformers the `RemoveAllConformers` method is run 
        # first. The shifted conformer is then given to the rdkit 
        # molecule, which is returned.
        new_mol.RemoveAllConformers()
        new_mol.AddConformer(conformer)
        return new_mol 

    def update_from_mae(self, mae_path=None):
        """
        Updates data in attributes to match what is held a .mae file.
        
        The molecule can be updated from a random .mae file held
        anywhere. Alterntatively, if no path is specified it is assumed
        the .mae file has the same path and name as 
        `self.prist_mol_file` only with .mol replaced with .mae.
        
        Parameters
        ----------
        mae_path : str (default = None)
            The full path of the .mae file from which the attributes
            should be updated. If ``None`` the .mae file is assumed to
            have the same path and name as `self.prist_mol_file` only
            with .mol replaced by .mae.
        
        Modifies
        --------
        self.prist_mol : rdkit.Chem.rdchem.Mol
            The rdkit molecule held in this attribute is changed so that
            it matches the moleclue held in the .mae file.
            
        self.prist_mol_file's content
            The content in this file is changed to match the content in
            the .mae file.
            
        Returns
        -------
        None : NoneType
        
        """
        
        if mae_path is None:
            mae_path = os.path.splitext(self.prist_mol_file)[0] + '.mae'
        
        self.prist_mol = mol_from_mae_file(mae_path)
        self.write_mol_file('prist')

    def write_mol_file(self, rdkit_mol_type, path=None):
        """
        Writes a V3000 ``.mol`` file of the macromolecule.

        The heavy or pristine molecule can be written via the argument
        `rdkit_mol_type`. The molecule is written to the location in the
        `prist_mol_file`/`heavy_mol_file` attribute. The function uses
        the structure of the rdkit molecules held in the `prist_mol` and
        `heavy_mol` attributes as the basis for what is written to the
        file.

        This bypasses the need to use rdkit's writing functions, which
        have issues with macromolecules due to poor ring finding and
        sanitization issues.

        Parameters
        ----------
        rdkit_mol_type : str
            Allowed values for this parameter 'prist' and 'heavy'.
            
        path : str (default = None)
            If the .mol file is to be written to a directory other than
            the one in `prist_mol_file` or `heavy_mol_file`, it should
            be written here.
        
        Modifies
        --------
        prist_mol_file's content
            If the string 'prist' was passed, the content in this
            ``.mol`` file will be replaced with the structure of the 
            current rdkit molecule in `prist_mol`.
            
        heavy_mol_file's content
            If the string 'heavy' was passed, the content in this
            ``.mol`` file will be replaced with the structure of the 
            current rdkit molecule in `heavy_mol`.
                
        Returns
        -------
        None : NoneType
        
        Raises
        ------
        ValueError
            If the `rdkit_mol_type` value is not either 'prist' or 
            'heavy'.
        
        """
        
        if rdkit_mol_type == 'prist':
            rdkit_mol = self.prist_mol
            file_name = self.prist_mol_file

        elif rdkit_mol_type == 'heavy':
            rdkit_mol= self.heavy_mol
            file_name = self.heavy_mol_file
            
        else:
            raise ValueError(("The argument `rdkit_mol_type` must be "
                              "either 'prist' or 'heavy'."))

        file_name = os.path.splitext(file_name)[0] + '.mol'
            
        main_string = ("\n"
                       "     RDKit          3D\n"
                       "\n"
                       "  0  0  0  0  0  0  0  0  0  0999 V3000\n"
                       "M  V30 BEGIN CTAB\n"
                       "M  V30 COUNTS {0} {1} 0 0 0\n"
                       "M  V30 BEGIN ATOM\n"
                       "!!!ATOM!!!BLOCK!!!HERE!!!\n"
                       "M  V30 END ATOM\n"
                       "M  V30 BEGIN BOND\n"
                       "!!!BOND!!!BLOCK!!!HERE!!!\n"
                       "M  V30 END BOND\n"
                       "M  V30 END CTAB\n"
                       "M  END\n"
                       "\n"
                       "$$$$\n")

        # id atomic_symbol x y z
        atom_line = "M  V30 {0} {1} {2:.4f} {3:.4f} {4:.4f} 0\n"
        atom_block = ""        
        
        # id bond_order atom1 atom2
        bond_line = "M  V30 {0} {1} {2} {3}\n"
        bond_block = ""
        
        
        main_string = main_string.format(rdkit_mol.GetNumAtoms(),
                                         rdkit_mol.GetNumBonds())
                                         
        for atom in rdkit_mol.GetAtoms():
            atom_id = atom.GetIdx()
            atom_sym = periodic_table[atom.GetAtomicNum()]
            x, y, z = self.atom_coords(rdkit_mol_type, atom_id)
            atom_block += atom_line.format(atom_id+1, atom_sym, x, y, z)
            
        for bond in rdkit_mol.GetBonds():
            bond_id = bond.GetIdx()
            atom1_id = bond.GetBeginAtomIdx() + 1
            atom2_id = bond.GetEndAtomIdx() + 1
            bond_order = int(bond.GetBondTypeAsDouble())
            bond_block += bond_line.format(bond_id, bond_order, 
                                           atom1_id, atom2_id)

        main_string = main_string.replace("!!!ATOM!!!BLOCK!!!HERE!!!\n",
                                          atom_block)
        main_string = main_string.replace("!!!BOND!!!BLOCK!!!HERE!!!\n",
                                          bond_block)
        
        if path:
            base_name = os.path.basename(file_name)
            file_name = os.path.join(path, base_name)
        
        with open(file_name, 'w') as f:
            f.write(main_string)

@total_ordering        
class StructUnit(Molecule, metaclass=Cached):
    """
    Represents the building blocks of macromolecules examined by MMEA.
    
    ``Building blocks`` in this case refers to the smallest molecular 
    unit of the assembled molecules (such as cages) examined by MMEA. 
    This is not the be confused with building-blocks* of cages. 
    Building-blocks* of cages are examples of the ``building blocks`` 
    referred to here. To be clear, the ``StructUnit`` class represents 
    all building blocks of macromolecules, such as both the linkers and 
    building-blocks* of cages.
    
    To avoid confusion, in the documentation general building blocks 
    represented by ``StructUnit`` are referred to as `building blocks`` 
    while building-blocks* of cages are always referred to as 
    ``building-blocks*``. 
    
    The goal of this class is to conveniently store information about, 
    and perform operations on, single instances of the building blocks 
    used to form macromolecules. The class stores information regarding
    the rdkit instance of the building block and the location of its 
    molecular structure file. See the attributes section of this 
    docstring for a full list of information stored.
    
    This class also takes care of perfoming substitutions of the 
    functional groups in the building blocks via the 
    `_generate_functional_group_atoms` method. This method is 
    automatically invoked by the initialzer, so each initialized
    instance of ``StructUnit`` should automatically have all of the 
    attributes associated with the substituted version of the molecular 
    building block.
    
    More information regarding what operations the class supports can be
    found by examining the methods below.
    
    The ``StructUnit`` class is intended to be inherited from. As 
    mentioned before, ``StructUnit`` is a general building block. If one 
    wants to represent a specific building block, such as a linker or 
    building-block* (of a cage) a new class should be created. This new
    class will inherit ``StructUnit``. In this way, any operations 
    which apply generally to building blocks can be stored here and any
    which apply specifically to one kind of building block such as a 
    linker or building-block* can be placed within its own class. Due to
    inheritance, instances of a derived class are able to use operations
    of the parent class.
    
    Consider a useful result of this approach. When setting the 
    coordinates of linkers or building-blocks* during assembly, it is 
    necessary to know if the molecule you are placing is a 
    building-block* or linker. This is because a building-block* will go 
    on vertex (in this example, this may or may not be generally true)
    and a linker will go on an edge. 
    
    Assume that there is a ``Linker`` and a ``BuildingBlock`` class 
    which inherit from ``StructUnit``. As luck would have it, these 
    classes are implemented in MMEA. Even if nothing is present in the
    class definition itself, both classes will have all the 
    attributes and methods associated with ``StructUnit``. This means
    the positions of the rdkit molecules held in instances of those 
    classes can be shifted with the `shift_heavy_mol` method.
    
    By running:
    
        >>> isinstance(your_struct_unit_instance, Linker)
        
    you can determine if the molecule you are dealing with is an example
    of a building-block* or linker of a cage. As a result you can easily
    choose to run the correct function which shifts the coordinates 
    either to a vertex or an edge.
    
    A final note on the intended use. Each instance of an assembled 
    molecule class (such as an instance of the ``Cage`` class) will have
    1 instance of a ``StructUnit`` derived class for each type of 
    building block present. For example, if a cage is made by comining 4 
    of one kind of building-block* with 6 of some kind of linker only 
    one instance of ``BuildingBlock`` and one instance of ``Linker`` is 
    to be held. The fact that there are 4 of one kind of building-block*
    arranged in some way to form a  cage is the ``Cage`` instance's 
    `topology` attribute's problem. If however the cage was for some 
    reason built from 2 of one kind of building-block*, 2 of another 
    kind of building-block* and 6 of 1 type of linker, then 2 
    ``BuildingBlock`` instances would need to be held (and 1 
    ``Linker``).
    
    In summary, the intended use of this class is to answer questions
    such as (not exhaustive):
        
        > What basic structural units were used in the assembly of this 
          cage?
        > Which functional group was substituted in building-blocks*
          of this cage? 
        > Which atom was substituted for which in the linker? (Note that
          this question is delegated to the ``FGInfo`` instance held in 
          the `func_grp` attribute of a ``StructUnit`` instance)
        > Where is the ``.mol`` file represnting a single 
          building-block* of the cage located?
        > Where is the ``.mol`` file represnting the a single 
          building-block* of the cage, after it has been substituted 
          with a heavy atom, located?
        > Give me an rdkit instance of the molecule which represents the
          building-block* of a cage. Before and after 
          it has been substituted.
        > Give me an rdkit instance of the molecule which represents a
          a single linker of a cage, at postion ``(x,y,z)``.
          
    Questions which this class should not answer include:
    
        > How many building-blocks* does this cage have? (Ask the 
          ``Cage`` instance.)
        > What is the position of a linker within this cage? (Ask the 
          ``Cage`` instance.)
        > Create a bond between this ``Linker`` and ``BuildingBlock``. 
          (Ask the ``Cage`` instance.)
          
    A good guide is to ask ``Can this question be answered by examining
    a single building block molecule in and of itself?``. 
    
    This should be kept in mind when extending MMEA as well. If a 
    functionality which only requires a building block ``in a vaccuum`` 
    is to be added, it should be placed here. If it requires the 
    building blocks relationship to other objects there should be a 
    better place for it (if not, make one). 

    PS. The ``StructUnit`` class supports ordering so that the parameter
    `building_blocks` in the ``MacroMolecule`` initializer can be
    sorted. This is necessary to correctly implement caching of
    ``MacroMolecule`` instances. See the ``CachedMacroMol`` class for 
    more details. Do not use comparison operations of this class outside 
    of this function.

    PPS. The ``StructUnit`` class is itself cached via the ``Cached``
    metaclass.
    
    Attributes
    ----------
    prist_mol_file : str
        The full path of the molecular structure file holding the 
        unsubstituted molecule. The supported file formats are the keys 
        in `init_funcs` dictionary in the ``__init__()`` method. As long 
        as a file of one of these types is provided, MMEA will 
        automatically use the correct initialization function.
        
    prist_mol : rdkit.Chem.rdchem.Mol
        This is an ``rdkit molecule type``. It is the rdkit instance
        of the molecule held in `prist_mol_file`.
        
    heavy_mol_file : str
        The full path of the ``.mol`` file (V3000) holding the 
        substituted molecule. This attribute is initialized by the 
        initializer indirectly when it calls the `_generate_heavy_attrs` 
        method. 
    
    heavy_mol : rdkit.Chem.rdchem.Mol
        The rdkit instance of the substituted molecule. Generated by 
        the initializer when it calls the `_generate_heavy_attrs` 
        method.
        
    func_grp : FGInfo
        This attribute holds an instance of ``FGInfo``. The ``FGInfo``
        instance holds the information regarding which functional group
        was substituted in the pristine molecule and which atom was 
        substituted for which. Furthermore, it also holds the atomic 
        numbers of the atom which was substitued and the one used in its 
        place. For details on how this information is stored see the 
        ``FGInfo`` class string.
        
    heavy_ids : list of ints
        A list holding the atom ids of the substituted atoms in the 
        heavy molecule. The ids correspond to the ids in the rkdit 
        molecule.

    energy : Energy
        This attribute handles information about the energy of the 
        instance. See the documentation of ``Energy`` to see what is
        available.  
        
    optimized : bool (default = False)
        A flag to monitor whether an optimization has been performed on
        the molecule.
    
    """
    
    def __init__(self, prist_mol_file, minimal=False):
        """
        Initializes a ``StructUnit`` instance.
        
        Parameters
        ----------
        prist_mol_file : str
            The full path of the molecular structure file holding the 
            building block.
            
        minimal : bool (default = False)
            If ``True`` the full initialization is not carried out. It
            stops the initialization before the functional groups in the
            molecule are replaced by heavy atoms.
            
        """
        
        # Heavy versions of ``StructUnit`` instances will be placed
        # in a folder named ``HEAVY``. This folder is placed in the 
        # current working directory. The iterator checks to see if the
        # folder already exists. If not, it returns ``False`` instead
        # of the folder name.
        heavy_dir = next((name for name in os.listdir(os.getcwd()) if 
                        os.path.isdir(name) and "HEAVY" == name), False)
        # If ``HEAVY`` folder not created yet, make it.
        if not heavy_dir:        
            os.mkdir("HEAVY")
        
        # This dictionary maps file extensions to the initialization 
        # functions.
        init_funcs = {'.mol' : partial(chem.MolFromMolFile, 
                                    sanitize=False, removeHs=False), 
                   '.mol2' : partial(chem.MolFromMol2File, 
                                     sanitize=False, removeHs=False),
                    '.mae' : mol_from_mae_file,
                    '.pdb' : partial(chem.MolFromPDBFile,
                                     sanitize=False, removeHs=False),
                    'mol2_block' : partial(chem.MolFromMol2Block,
                                           sanitize=False,
                                           removeHs=False)}
 
        self.prist_mol_file = prist_mol_file
        _, ext = os.path.splitext(prist_mol_file)
#        ext = 'mol2_block'
        if ext not in init_funcs:
            raise TypeError(
            'Unable to initialize from "{}" files.'.format(ext))
                                     
        self.prist_mol = init_funcs[ext](prist_mol_file)
        print(init_funcs[ext], ext)
        self.energy = Energy(self)
        self.optimized = False        

        # Check for minimal initialization.
        if minimal:
            return
        
        # Define a generator which yields an ``FGInfo`` instance from
        # the `FGInfo.functional_group_list`. The yielded ``FGInfo``
        # instance represents the functional group found on the pristine
        # molecule used for initialization. The generator determines 
        # the functional group of the molecule from the path of its 
        # ``.mol`` file. 
        
        # The database of precursors should be organized such that any 
        # given ``.mol`` file has the name of its functional group in
        # its path. Ideally, this will happen because the ``.mol`` file
        # is in a folder named after the functional group the molecule 
        # in the ``.mol`` file contains. This means each ``.mol`` file 
        # should have the name of only one functional group in its path. 
        # If this is not the case, the generator will return the 
        # functional group which appears first in 
        # `FGInfo.functional_group_list`.
        
        # Calling the ``next`` function on this generator causes it to
        # yield the first (and what should be the only) result. The
        # generator will return ``None`` if it does not find the name of
        # a functional group in the path of the ``.mol`` file.
        self.func_grp = next((x for x in 
                                FGInfo.functional_group_list if 
                                x.name in prist_mol_file), None)
        self.heavy_ids = []        
        
        # Calling this function generates all the attributes assciated
        # with the molecule the functional group has been subtituted
        # with heavy atoms.
        self._generate_heavy_attrs()

    def _delete_tag_ids(self, heavy_atom):
        """
        Returns the ids of neighbor atoms tagged for deletion.
        
        These are neighbors to `heavy_atom`.

        Parameters
        ----------
        heavy_atom : rdkit.Chem.rdchem.Atom
            The heavy atom who's neighbors should be deleted.
            
        Returns
        -------
        list of ints
            The ids of neighbor atoms which have been tagged for
            deletion.
        
        """
        
        heavy_id = heavy_atom.GetIdx()
        
        del_ids = []
        
        # For each deletion tag, go through all the neighors that
        # `heavy_atom` has. If the atom's element and bond type matches
        # the deletion tag, add the atom's id to the ``del_ids`` list.        
        for del_atom in self.func_grp.del_tags:
            for n in heavy_atom.GetNeighbors():
                n_id = n.GetIdx()

                bond = self.heavy_mol.GetBondBetweenAtoms(heavy_id, 
                                                          n_id)
                bond_type = bond.GetBondType()

                if (n.GetAtomicNum() == del_atom.atomic_num and
                    bond_type == del_atom.bond_type and
                    n_id not in del_ids):
                    del_ids.append(n_id)
                    break
        
        return del_ids

    def _generate_heavy_attrs(self):
        """
        Adds attributes associated with a substituted functional group.
        
        This function is private because it should not be used outside 
        of the initializer.
        
        The function creates rdkit molecule where all of the functional
        group atoms have been replaced by the desired atom as indicated
        in the ``FGInfo`` instance of `self`. It also creates ``.mol``
        file holding this molecule and a ``SMILES`` string representing
        it.

        Modifies
        --------
        self : StructUnit
            Adds the `heavy_mol`, `heavy_mol_file` and `heavy_smiles`
            attributes to ``self``. It also adds the ids of heavy atoms
            to `heavy_ids` by calling the `_make_atoms_heavy_in_heavy`
            method.
        
        Returns
        -------
        None : NoneType                

        """
        
        # In essence, this function first finds all atoms in the 
        # molecule which form a functional group. It then switches the 
        # target atoms in the functional groups for heavy atoms and
        # deletes any tagged for deletion. This new molecule is then 
        # stored in the ``StructUnit`` instance in the form of an 
        # ``rdkit.Chem.rdchem.Mol``, a SMILES string and a ``.mol`` file 
        # path.
        
        # First create a copy of the ``rdkit.Chem.rdchem.Mol`` instance
        # representing the pristine molecule. This is so that after 
        # any changes are made, the pristine molecule's data is not 
        # corrupted. This second copy which will turn into the 
        # substituted ``rdkit.Chem.rdchem.Mol`` will be operated on.
        self.heavy_mol = chem.Mol(self.prist_mol)      

        # Subtitutes the relevant functional group atoms in `heavy_mol`
        # for heavy atoms and deletes Hydrogen atoms of the functional
        # group as well as any other atoms tagged for deletion.
        self._make_atoms_heavy()

        # Change the pristine ``.mol`` file name to include the word
        # ``HEAVY_`` at the end. This generates the name of the 
        # substituted version of the ``.mol`` file.
        heavy_file_name = list(os.path.splitext(self.prist_mol_file))
        heavy_file_name.insert(1,'HEAVY')
        heavy_file_name.insert(2, self.func_grp.name)
        self.heavy_mol_file = '_'.join(heavy_file_name)
        self.heavy_mol_file = os.path.split(self.heavy_mol_file)[1]
        self.heavy_mol_file = os.path.join(os.getcwd(), "HEAVY", 
                                           self.heavy_mol_file)
        #self.write_mol_file('heavy')      

    def _make_atoms_heavy(self):
        """
        Converts functional group in `heavy_mol` to substituted version.

        The function changes the element of the central atom of the
        functional group and deletes any hydrogen atoms in the
        functional group. It also deletes any atoms in the functional
        group that have been tagged for deletion.        
        
        Modifies
        --------
        heavy_mol : rdkit.Chem.rdchem.Mol
            The rdkit molecule has atoms removed and converted to 
            different elements, as described in the docstring.

        heavy_ids : list
            Adds the ids of heavy atoms to this list.

        Returns
        -------
        None : NoneType

        """        
         
        # Go through all atom ids corresponding to a functional group
        # atoms. If the element of the atom corresponds to the atom
        # which must be substituted, do so. For each substituted atom
        # check if any of its neighbors are tagged for deletion and
        # tag all the neighboring Hydrogen atoms for deletion. An atom
        # tagged for deletion will have its id added to ``del_ids``.
        del_ids = []
        for atom_id in flatten(self.find_functional_group_atoms()):
            atom = self.heavy_mol.GetAtomWithIdx(atom_id)
            if atom.GetAtomicNum() == self.func_grp.target_atomic_num:            
                atom.SetAtomicNum(self.func_grp.heavy_atomic_num)
                
                del_ids.extend(self._delete_tag_ids(atom))
                for n in atom.GetNeighbors():
                    if n.GetAtomicNum() == 1:
                        del_ids.append(n.GetIdx())

        # Make an EditableMol and delete all of its atoms which have
        # been tagged for deletion. Delete ids with larger ids first to
        # prevent rdkit from raising range errors.      
        editable_mol = chem.EditableMol(self.heavy_mol)
        for del_id in sorted(del_ids, reverse=True):
            editable_mol.RemoveAtom(del_id)
               
        self.heavy_mol = editable_mol.GetMol()
        
        # Take note of heavy atom ids.
        for atom in self.heavy_mol.GetAtoms():
            if atom.GetAtomicNum() == self.func_grp.heavy_atomic_num:
                self.heavy_ids.append(atom.GetIdx())

    def _set_heavy_mol_orientation(self, start, end):
        """
        Rotates from `start` to `end`.
        
        Given two direction vectors, `start` and `end`, this method
        applies the rotation required to transform `start` to `end` on 
        the heavy molecule. The rotation occurs about the centroid
        of the heavy atoms. 
        
        For example, if the `start` and `end` vectors
        are 45 degrees apart, a 45 degree rotation will be applied to
        heavy molecule. The rotation will be along the appropriate axis.
        
        This method will likely have counterparts in derived classes of 
        ``StructUnit``. The counterparts will probably use a default
        `start` vector and will not be private.  This prevents 
        overwriting and means both versions of the function will be 
        available. For example, the ``Linker`` class will have the
        default `start` vector as the direction vector between the 2
        heavy atoms. This means that running the function on a 
        ``Linker`` instance will align the heavy atoms with the vector
        `end`.

        On the other hand, the ``BuildingBlock`` class will use the 
        normal to the plane formed by the heavy atoms as the `start`
        vector. This means that in ``BuildingBlock`` molecules the
        normal to the plane will be aligned with `end` when the method
        is run.
        
        As the above examples demonstrate, the great thing about this 
        method is that you as long as you can associate a geometric 
        feature of the molecule with a vector, then the molecule can be 
        roatated so that this vector is aligned with `end`. The defined 
        vector can be virtually anything. This means that any geomteric 
        feature of the molecule can be easily aligned with any arbitrary 
        direction.

        This method modifies the rdkit molecule in `heavy_mol`. It also
        returns of a copy of this rdkit molecule.
        
        Parameters
        ----------
        start : numpy.array
            A vector which is to be rotated so that it transforms to the
            `end` vector.
        
        end : numpy.array
            This array holds the directional vector along which the 
            heavy atoms in the linker should be placed.
            
        Modifies
        --------
        heavy_mol : rdkit.Chem.rdchem.Mol   
            The conformer in this rdkit instance is changed due to
            rotation of the molecule about the centroid of the heavy
            atoms.
        
        Returns
        -------
        rdkit.Chem.rdchem.Mol
            An rdkit molecule instance of the rotated molecule. This is 
            a copy of the rdkit molecule in `heavy_mol`.
        
        """
        
        # Normalize the input direction vectors.
        start = normalize_vector(start)
        end = normalize_vector(end)
        
        # Record the position of the molecule then translate the heavy
        # atom centroid to the origin. This is so that the rotation
        # occurs about this point.
        og_center = self.heavy_atom_centroid()
        self.set_heavy_atom_centroid(np.array([0,0,0])) 
        
        # Get the rotation matrix.
        rot_mat = rotation_matrix(start, end)
        
        # Apply the rotation matrix to the atomic positions to yield the
        # new atomic positions.
        new_pos_mat = np.dot(rot_mat, self.position_matrix('heavy'))

        # Set the positions in the heavy rdkit molecule.
        self.set_position_from_matrix('heavy', new_pos_mat)
        self.set_heavy_atom_centroid(og_center)

        return chem.Mol(self.heavy_mol)

    def centroid_centroid_dir_vector(self):
        """
        Returns the direction vector between the 2 molecular centroids.
        
        This method uses the substituted version of the molecule. The
        two centroids are the molecular centroid and the centroid of the
        heavy atoms only.
        
        Returns
        -------
        numpy.array
            The normalized direction vector running from the centroid of
            the heavy atoms to the molecular centroid.
        
        """
    
        return normalize_vector(self.centroid('heavy') - 
                                self.heavy_atom_centroid())

    def find_functional_group_atoms(self):
        """
        Returns a container of atom ids of atoms in functional groups.

        The ``StructUnit`` instance (`self`) represents a molecule. 
        This molecule is in turn represented in rdkit by a 
        ``rdkit.Chem.rdchem.Mol`` instance. This rdkit molecule instance 
        is held in the `prist_mol` attribute. The rdkit molecule 
        instance is made up of constitutent atoms which are
        ``rdkit.Chem.rdchem.Atom`` instances. Each atom instance has its
        own id. These are the ids contained in the tuple returned by 
        this function.   

        Returns
        -------
        tuple of tuples of ints
            The form of the returned tuple is:
            ((1,2,3), (4,5,6), (7,8,9)). This means that all atoms with
            ids 1 to 9 are in a functional group and that the atoms 1, 2
            and 3 all form one functional group together. So do 4, 5 and 
            5 and so on.

        """
        
        # Generate a ``rdkit.Chem.rdchem.Mol`` instance which represents
        # the functional group of the molecule.        
        func_grp_mol = chem.MolFromSmarts(self.func_grp.smarts_start)
        
        # Do a substructure search on the the molecule in `prist_mol`
        # to find which atoms match the functional group. Return the
        # atom ids of those atoms.
        return self.prist_mol.GetSubstructMatches(func_grp_mol)        

    def heavy_atom_centroid(self):
        """
        Returns the centroid of the heavy atoms.

        This is the centroid if only the heavy atoms are considered.

        Returns
        -------
        numpy.array
            A numpy array holding the midpoint of the heavy atoms.
        
        """

        centroid = sum(self.atom_coords('heavy', x) for x in 
                                                        self.heavy_ids) 
        return np.divide(centroid, len(self.heavy_ids))

    def set_heavy_atom_centroid(self, position):
        """
        Shift heavy molecule so heavy atom centroid is on `position`.         
        
        The entire molecule is shifted but it is the centroid of the
        heavy atoms that is placed on the target `position`.     
        
        This method changes the position of the rdkit molecule in 
        `heavy_mol`. It also returns a of copy of the shifted rdkit 
        molecule.
        
        Parameters
        ----------
        position : numpy.array
            A numpy array holding the desired the position. It holds the
            x, y and z coordinates, respectively.
            
        Modifies
        --------
        heavy_mol : rdkit.Chem.rdchem.Mol   
            The position of the molecule in this rdkit instance is
            changed, as described in this docstring.
            
        Returns
        -------
        rdkit.Chem.rdchem.Mol 
            A copy of the rdkit molecule after it has been shifted.
            
        """
        
        center = self.heavy_atom_centroid()
        shift = position - center
        new_conf = self.shift('heavy', shift).GetConformer()

        self.heavy_mol.RemoveAllConformers()
        self.heavy_mol.AddConformer(new_conf)
        
        return chem.Mol(self.heavy_mol)

    def similar_molecules(self, database):
        """
        Returns molecules from `database` ordered by similarity.
        
        This method uses the Morgan fingerprints of radius 4 to evaluate 
        how similar the molecules in `database` are to `self`.
        
        Parameters
        ----------
        database : str
            The full path of the database from which the molecules are
            checked for similarity.
            
        Returns
        -------
        list of tuples of form (float, str)
            The float is the similarity of a given molecule in the
            database to `self` while the str is the full path of the 
            .mol file of that molecule.
        
        """
        # First get the fingerprint of `self`.
        chem.GetSSSR(self.prist_mol)
        self.prist_mol.UpdatePropertyCache(strict=False)
        fp = ac.GetMorganFingerprint(self.prist_mol, 4)
        
        # For every .mol file in the database create an rdkit molecule.
        # Place these in a list.
        mols = []
        for file in os.listdir(database):
            path = os.path.join(database, file)
            # Ignore files which are not .mol and the .mol file of the
            # molecule itself.
            if not path.endswith('.mol') or path == self.prist_mol_file:
                continue
            mol = chem.MolFromMolFile(path, 
                                      sanitize=False, 
                                      removeHs=False)
            chem.GetSSSR(mol)
            mol.UpdatePropertyCache(strict=False)
            mol_fp = ac.GetMorganFingerprint(mol, 4)
            similarity = DataStructs.DiceSimilarity(fp, mol_fp)
            mols.append((similarity, path))
        
        return sorted(mols, reverse=True)

    def __eq__(self, other):
        return self.prist_mol_file == other.prist_mol_file
        
    def __lt__(self, other):
        return self.prist_mol_file < other.prist_mol_file
        
    def __hash__(self):
        return id(self)
    
    def __str__(self):
        return self.prist_mol_file
    
    def __repr__(self):
        repr_ =  "{0!r}".format(type(self))
        repr_ = repr_.replace(">", 
        ", prist_mol_file={0.prist_mol_file!r}>".format(self))
        repr_ = repr_.replace("class ", "class=")
        return repr_
  
class StructUnit2(StructUnit):
    """
    Represents building blocks with 2 functional groups.
    
    """
    
    def set_heavy_mol_orientation(self, end):
        """
        Rotate heavy molecule so heavy atoms lie on `end`.     
        
        The molecule is rotated about the centroid of the heavy atoms.
        It is rotated so that the direction vector running between the
        2 heavy atoms is aligned with the vector `end`.        
        
        Parameters
        ----------
        end : numpy.array
            The vector with which the molecule's heavy atoms should be
            aligned.
        
        Modifies
        --------        
        heavy_mol : rdkit.Chem.rdchem.Mol   
            The conformer in this rdkit instance is changed due to
            rotation of the molecule about the centroid of the heavy
            atoms.
        
        Returns
        -------
        rdkit.Chem.rdchem.Mol
            An rdkit molecule instance of the rotated molecule. This is 
            a copy of the rdkit molecule in `heavy_mol`.        
        
        """
        
        start = next(self.heavy_direction_vectors())
        return StructUnit._set_heavy_mol_orientation(self, start, end)
        
    def minimize_theta(self, vector, axis):
        """
        Rotates linker about `axis` to minimze theta with `vector`.
        
        The linker is iteratively rotated so that its heavy atom vector
        is as close as possible to `vector`.        
        
        Parameters
        ----------
        vector : numpy.array
            The vector to which the distance should be minimized.
            
        axis : numpy.array
            The direction vector along which the rotation happens.
        
        Returns
        -------
        None : NoneType        
        
        """
        
        vector = normalize_vector(vector)
        axis = normalize_vector(axis)
        
        # Size of iterative step in radians.
        step = 0.17
        
        theta = vector_theta(self.centroid_centroid_dir_vector(),
                             vector)
                             
        # First determine the direction in which iteration should occur.
        self.rotate2(step, axis)
        theta2 = vector_theta(self.centroid_centroid_dir_vector(),
                             vector)       
        if theta2 > theta:
            axis = np.multiply(axis, -1)
            
        prev_theta = theta2
        while True:
            self.rotate2(step, axis)
            theta = vector_theta(self.centroid_centroid_dir_vector(),
                             vector)
            
            if theta > prev_theta:
                axis = np.multiply(axis, -1)
                self.rotate2(step, axis)
                break
            
            prev_theta = theta

      
class StructUnit3(StructUnit):
    """
    Represents building blocks with 3 functional groups.
    
    """

    def heavy_plane(self):
        """
        Returns the coefficients of the plane formed by heavy atoms.
        
        A plane is defined by the scalar plane equation,
            
            ax + by + cz = d.
        
        This method returns the a, b, c and d coefficients of this 
        equation for the plane formed by the heavy atoms. The 
        coefficents a, b and c decribe the normal vector to the plane.
        The coefficent d is found by substituting these coefficients
        along with the x, y and z variables in the scalar equation and
        solving for d. The variables x, y and z are substituted by the
        coordinate of some point on the plane. For example, the position
        of one of the heavy atoms.
        
        Returns
        -------
        numpy.array
            This array has the form [a, b, c, d] and represents the 
            scalar equation of the plane formed by the heavy atoms.
        
        References
        ----------
        http://tutorial.math.lamar.edu/Classes/CalcIII/EqnsOfPlanes.aspx                
        
        """
        
        heavy_coord = self.atom_coords('heavy', self.heavy_ids[0])
        d = np.multiply(np.sum(np.multiply(self.heavy_plane_normal(), 
                                           heavy_coord)), -1)
        return np.append(self.heavy_plane_normal(), d)

    def heavy_plane_normal(self):
        """
        Returns the normal vector to the plane formed by heavy atoms.
        
        The normal of the plane is defined such that it goes in the
        direction toward the centroid of the building-block*.        
        
        Returns
        -------        
        numpy.array
            A unit vector which describes the normal to the plane of the
            heavy atoms.
        
        """
        
        if sum(1 for _ in self.heavy_direction_vectors()) < 2:
            raise ValueError(("StructUnit3 molecule "
                             "has fewer than 3 functional groups."))
        
        v1, v2 = itertools.islice(self.heavy_direction_vectors(), 2)
    
        normal_v = normalize_vector(np.cross(v1, v2))
        
        theta = vector_theta(normal_v, 
                             self.centroid_centroid_dir_vector())
                             
        if theta > np.pi/2:
            normal_v = np.multiply(normal_v, -1)
        
        return normal_v
        
    def set_heavy_mol_orientation(self, end):
        """
        Rotates heavy molecule so plane normal is aligned with `end`.

        Here ``plane normal`` referes to the normal of the plane formed
        by the heavy atoms in the substituted molecule. The molecule
        is rotated about the centroid of the heavy atoms. The rotation
        results in the normal of their plane being aligned with `end`.

        Parameters
        ----------
        end : numpy.array
            The vector with which the normal of plane of heavy atoms 
            shoould be aligned.
        
        Modifies
        --------
        heavy_mol : rdkit.Chem.rdchem.Mol   
            The conformer in this rdkit instance is changed due to
            rotation of the molecule about the centroid of the heavy
            atoms.        

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            An rdkit molecule instance of the rotated molecule. This is 
            a copy of the rdkit molecule in `heavy_mol`.
            
        """
        
        start = self.heavy_plane_normal()
        return StructUnit._set_heavy_mol_orientation(self, start, end)        

@total_ordering
class MacroMolecule(Molecule, metaclass=CachedMacroMol):
    """
    A class for MMEA assembled macromolecules.
    
    The goal of this class is to represent an individual used by the GA.
    As such, it holds attributes that are to be expected for this
    purpose. Mainly, it has a fitness value stored in its `fitness` 
    attribute and a genetic code - as defined by its `building_blocks` 
    and  `topology` attributes. If a change is made to either of these 
    attributes, they should describe a different macromolecule. On the
    other hand, the same attributes should always describe the same 
    macromolecule.
    
    Because of this, as well as the computational cost associated with
    macromolecule initialization, instances of this class are cached. 
    This means that providing the same arguments to the initializer will 
    not build a different instance with the same attribute values. It 
    will yield the original instance, retrieved from memory.
    
    To prevent bloating this class, any information that can be 
    categorized is. For example, storing information that concerns
    building blocks ``in a vacuum`` is stored in ``StructUnit``
    instances. Equally, manipulations of such data is also performed by
    those instances. Similarly, anything to do with topolgy should be 
    held by a ``Topology`` instance in the topology attribute. There is
    a notable exception to this however. This happens when retrieving 
    topological information directly from rdkit molecule instances 
    representing the macromolecule. Examples include the 
    information about atomic coordinates is stored in the rdkit molecule 
    instances, which are stored directly by this class.
    
    It should also be noted that only a single copy of each 
    ``StructUnit`` instance representing a specific building block needs
    to be held. How many of such building blocks are need to assemble
    the cage is the handled by the ``Topology`` class, which only needs
    a single copy of each building block to work with.    
    
    If new inormation associated with macromolecules, but not directly 
    concerning them as a whole, is to be added at some point in the 
    future, and that information can be grouped together in a logical 
    category, a new class should be created to store and manipulate this 
    data. It should not be given to the macromolecule directly. 
    Alternatively if more information to do with one of the already 
    categories, it should be added there. The attribute 
    `building_blocks` and its composing ``StructUnit`` instaces are an
    example of this approach.
    
    However information dealing with the cage as a whole can be added
    directly to attributes of this class. You can see examples of such 
    attributes below. Simple identifiers such as ``.mol`` files and 
    ``SMILES`` strings do not benefit from being grouped together. 
    (Unless they pertain to specific substructures within the cages such 
    as linkers and building-blocks* - as mentioned before.) Topology is 
    an exception to this because despite applying to the cage as a 
    whole, it a complex aspect with its own functions and data. 
    
    The goal is simplicity. Having too many categories causes unneeded
    complexity as does having too few.
    
    This class is not intended to be used directly but should be 
    inherited by subclasses representing specific macromolecules. The
    ``Cage`` and ``Polymer`` classes are examples of this. Any 
    information or methods that apply generally to all macromolecules
    should be defined within this class while specific non-general data
    should be included in derived classes.
    
    This class also supports comparison operations, these act on the 
    fitness value assiciated with a macromolecule. Comparison operations 
    not explicitly defined are included via the ``total_ordering`` 
    decorator. For other operations and methods supported by this class 
    examine the rest of the class definition.

    Finally, a word of caution. The equality operator ``==`` compares 
    fitness values. This means two macromolecules, made from different 
    building blocks, can compare equal if they happen to have the same 
    fitness. The operator is not to be used to check if one 
    macromolecule is the same structurally as another. To do this check 
    use the `same` method. This method may be overwritten in derived 
    classes, as necessary. In addition the ``is`` operator is 
    implemented as is default in Python. It compares whether two objects 
    are in the same location n memory. Because the ``MacroMolecule`` 
    class is cached the ``is`` operator could in principle be used 
    instead of the `same` method (including in derived classes). 
    However, this is not intended use and is not guuranteed to work in 
    future implementations. If caching stops being implemented such code 
    would break.
    
    Optimization of structures of ``MacroMolecule`` instances is not
    done by this class. This is because in order to run optimization
    functions in parallel, they cannot be defined as methods. As a
    result optimizations are implemented functionally in the
    ``optimization.py`` module.

    Attributes
    ----------
    building_blocks : list of ``StructUnit`` instances
        This attribute holds ``StructUnit`` instances which represent
        the monomers forming the macromolecule. Only one ``StructUnit``
        instance is needed per monomer, even if multiples of a monomer 
        join up to form the macromolecule

    topology : A child class of ``Topology``
        This instance represents the topology of the macromolecule. Any 
        information to do with how individual building blocks of the 
        macromolecule are organized and joined up in space is held by 
        this attribute. For more details about what information and 
        functions this entails see the docstring of the ``Topology`` 
        class and its derived classes.

    topology_args : list (default = [])
        This attribue holds the initializer arguments for the topology 
        instance. This is stored so that exceptions can print all
        values required to make an identical copy of a ``MacroMolecule``
        instance..

    prist_mol_file : str
        The full path of the ``.mol`` file holding the pristine version
        of the macromolecule.

    prist_mol : rdkit.Chem.rdchem.Mol
        An rdkit molecule instance representing the macromolecule.
        
    heavy_mol_file : str
        The full path of the ``.mol`` file holding the substituted
        version of the macromolecule.

    heavy_mol : rdkit.Chem.rdchem.Mol
        A rdkit molecule instance holding the substituted version of the
        macromolecule.

    optimized : bool (default = False)
        This is a flag to indicate if a molecule has been previously
        optimized. Optimization functions set this flag to ``True``
        after an optimization.
        
    energy : Energy
        This attribute handles information about the energy of the 
        instance. See the documentation of ``Energy`` to see what is
        available.        

    fitness : float (default = False)
        The fitness value of the macromolecule, as determined by the 
        chosen fitness function. This attribute is assigned by fitness
        functions and initialized with ``False``. If a scaling or
        normalization procedure is used, the fitness function will often
        assign to the attribute `unscaled_fitness` while normalization
        function will assign to `fitness`.
        
    unscaled_fitness = float or array (default = False)
        Fitness functions which couple with scaling or normalization 
        functions assign fitness values into this attribute.
        
    fitness_fail : bool (default = True)
        Used to indicate if the fitness function failed to calculate
        a value for the MacroMolecule. Not all fitness functions or
        normalization fucntions use this attribute but it can be useful.
        
    progress_params : list (default = False)
        Holds the fitness parameters which the GA should track to make
        progress plots. If the default ``False`` is used, the fitness
        value will be used.
        
    key : str
        The key used for caching the molecule. Necessary for 
        `update_cache` to work. This attribute is assigned by in the 
        `__call__` method of the ``CachedMacroMol`` class.
    
    """

    def __init__(self, building_blocks, topology, prist_mol_file, 
                 topology_args=None):
        """
        Initialize a ``MacroMolecule`` instance.
        
        The initialization is exectued inside a try block. This allows
        error handling for cases where cage initialization failed for
        some reason. In this case, when an exception occurs during
        initialization all the parameters which were provided to the 
        initializer are saved to a file ``failures.txt`` which should
        be located in the same directory as the ``output`` folder.
        
        Parameters
        ---------
        building_blocks : list of ``StructUnit`` instances
            A list of ``StructUnit`` instances which represent the 
            monomers forming the macromolecule.

        topology : A child class of ``Topology``
            The class which defines the topology of the macromolecule. 
            Such classes are defined in the topology module. The class 
            will be a child class which inherits ``Topology``.
        
        prist_mol_file : str
            The full path of the ``.mol`` file where the macromolecule
            will be stored.
            
        topology_args : list (default = None)
            Any additional arguments needed to initialize the topology
            class supplied in the `topology` argument.
            
        """
        
        try:
            self._std_init(building_blocks, topology, prist_mol_file, 
                                 topology_args)
            
        except Exception as ex:
            self.building_blocks = building_blocks
            self.topology = topology
            self.prist_mol = chem.Mol()
            self.prist_mol_file = prist_mol_file
            self.topology_args = topology_args
            MacroMolError(ex, self, 'During initialization.')

    def _std_init(self, building_blocks, topology, prist_mol_file, 
                 topology_args):
            
        if topology_args is None:
            topology_args = {}

        self.building_blocks = tuple(building_blocks)

        # A ``Topology`` subclass instance must be initialized with a 
        # copy of the cage it is describing.
        self.topology = topology(self, **topology_args)
        # The topology_args attribute is saved for error handling. See
        # MacroMolError class.
        self.topology_args = topology_args
        self.prist_mol_file = prist_mol_file

        # This generates the name of the heavy ``.mol`` file by adding
        # ``HEAVY_`` at the end of the pristine's ``.mol`` file's name. 
        heavy_mol_file = list(os.path.splitext(prist_mol_file))
        heavy_mol_file.insert(1,'HEAVY')        
        self.heavy_mol_file = '_'.join(heavy_mol_file) 
        
        self.energy = Energy(self)        
        
        # Ask the ``Topology`` instance to assemble/build the cage. This
        # creates the cage's ``.mol`` file all  the building blocks and
        # linkers joined up. Both the substituted and pristine versions.      
        self.topology.build()
    

        # Write the ``.mol`` files.
        self.write_mol_file('prist')
        self.write_mol_file('heavy')
                                             
        self.optimized = False

        # A numerical fitness is assigned by fitness functions evoked
        # by a ``Population`` instance's `GATools` attribute.
        self.fitness = False
        self.unscaled_fitness = False
        self.fitness_fail = True
        self.progress_params = False
        
        # Dump the molecule in case it causes an issue later on.
        self.dump(self.prist_mol_file.replace('.mol', '.dmp'))

    def same(self, other):
        """
        Check if the `other` instance has the same molecular structure.
        
        Parameters
        ----------
        other : MacroMolecule
            The ``MacroMolecule`` instance you are checking has the same 
            structure.
        
        Returns
        -------
        bool
            Returns ``True`` if the building-block*, linker and 
            topology of the cages are all the same.
        
        """
        
        # Compare the building blocks and topology making up the 
        # macromolecule. If these are the same then the cages have the 
        # same structure.
        return (self.building_blocks == other.building_blocks and 
                                    self.topology == other.topology)

    def update_cache(self):
        """
        Updates the caching dictionary so that it contains `self`.
        
        When an instance of ``MacroMolecule`` is first created it is
        cached. Using multiprocessing to perform optimizations returns
        modified copies of the cached molecules. In order to ensure that
        the cache points to the modified copies not to the originally
        initialized molecule, this method must be run.
        
        Returns
        -------
        None : NoneType
        
        """
        
        # The caching is done by the class. Access it and use its
        # ``_update_cache()`` method.
        cls = type(self)
        cls._update_cache(self)
    
    def __eq__(self, other):
        return self.fitness == other.fitness
        
    def __lt__(self, other):
        return self.fitness < other.fitness
    
    def __str__(self):
        return str({key: value for key, value in 
                                    self.__dict__.items() if 
                                    key in {'prist_mol_file', 
                                            'topology',
                                            'fitness',
                                            'optimized'}}) + "\n"
    
    def __repr__(self):
        return str(self)
        
    def __hash__(self):
        return id(self)

    """
    The following methods are inteded for convenience while 
    debugging or testing and should not be used during typical 
    execution of the program.
    
    """
    @classmethod
    def testing_init(cls, bb_str, lk_str, topology_str):
        key = (bb_str, lk_str, topology_str)
        if key in MacroMolecule._cache.keys():
            return MacroMolecule._cache[key]
        else:            
            cage = cls.__new__(cls)        
            cage.building_blocks = (bb_str, lk_str)
            cage.topology = topology_str
            cage.fitness = 3.14
            MacroMolecule._cache[key] = cage
            return cage


class Cage(MacroMolecule):
    """
    Used to represent molecular cages.
    
    Attributes
    ----------
    In addtion to those described in MacroMolecule, this class the
    following attributes.
    
    unscaled_fitness_vars : numpy.array
        An array holding the values of various contributions to the
        fitness of a cage, before they are scaled by the population
        averages. This attribute is created and modified the by the
        ``cage`` fitness function. See its documentation for more
        details.
    
    """

    @classmethod
    def init_fixed_bb(cls, bb_file, lk_db, topologies, prist_mol_file):
        bb = StructUnit3(bb_file)        
        
        while True:
            lk_file = np.random.choice(os.listdir(lk_db))
            if lk_file.endswith(".mol"):
                break
            
        lk_file = os.path.join(lk_db, lk_file)
        lk = StructUnit2(lk_file)
        
        topology = np.random.choice(topologies)        
        
        return cls((bb, lk), topology, prist_mol_file)
  
    @classmethod
    def init_random(cls, bb_db, lk_db, topologies, prist_mol_file):
        """
        Makes ``Cage`` from random building blocks and topology.
        
        Parameters
        ----------
        bb_db : str
        
        lk_db : str
        
        topologies : list of ``Topology`` child classes.
        
        prist_mol_file : str
        
        """
        
        while True:
            bb_file = np.random.choice(os.listdir(bb_db))
            if bb_file.endswith(".mol"):
                break
            
        bb_file = os.path.join(bb_db, bb_file)
        bb = StructUnit3(bb_file)
        
        while True:
            lk_file = np.random.choice(os.listdir(lk_db))
            if lk_file.endswith(".mol"):
                break
            
        lk_file = os.path.join(lk_db, lk_file)
        lk = StructUnit2(lk_file)
        
        topology = np.random.choice(topologies)
        
        return cls((bb, lk), topology, prist_mol_file)

class Polymer(MacroMolecule):
    """
    Used to represent polymers.
    
    """
    pass
