import numpy as np
from functools import partial
from multiprocessing import Pool
import itertools
import rdkit.Chem as chem

from .pyWindow import window_sizes
from .classes.exception import MacroMolError

def calc_fitness(func_data, population):
    """
    Calculates the fitness values of all members of a population.
    
    A fitness function should take a ``MacroMolecule`` instance and
    return a number representing its fitness. The assignement to the
    `fitness` attribute of a population member happens here, not by the
    fitness function.    
    
    Parameters
    ----------
    func_data : FunctionData
        A ``FunctionData`` instance representing the chosen fitness 
        function and any additional parameters it may require.
    
    population : Population
        The population whose members must have their fitness calculated.
        
    Returns
    -------
    None : NoneType    
    
    """

    # Get the fitness function object.
    func = globals()[func_data.name]

    # Apply the function to every member of the population.
    for macro_mol in population:
        try: 
            macro_mol.fitness = func(macro_mol, **func_data.params)
            print(macro_mol.fitness)
            
        except Exception as ex:
            MacroMolError(ex, macro_mol, 'During fitness calculation.')
            macro_mol.fitness = 0
            
def random_fitness(macro_mol):
    """
    Returns a random fitness value between 0 and 100.

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule to which a fitness value is to be assigned.
    
    Returns
    -------
    int
        An integer between 0 (including) and 100 (excluding).

    """

    return np.random.randint(0,100)

def cage(macro_mol, target_size):
    if isinstance(target_size, str):
        target_size = float(target_size)
        
    ideal_ratio = 1/(np.square(target_size))

    macro_mol.topology.windows = window_sizes(macro_mol.prist_mol_file)    
    
    if macro_mol.topology.windows != None:
        try:
            macro_mol.topology.max_window = max(macro_mol.topology.windows)
            if macro_mol.topology.max_window != 0:
                window_ratio = 1/(np.square(macro_mol.topology.max_window))
            else:
                window_ratio = 200
            window_ratio_tot = abs(ideal_ratio - window_ratio) * 100
            macro_mol.topology.window_ratio = window_ratio_tot
            
        except:
            macro_mol.topology.max_window = 0.5
            window_ratio = 1/(macro_mol.topology.max_window ** 2)
            window_ratio_tot = abs(ideal_ratio - window_ratio) * 100
            macro_mol.topology.window_ratio = window_ratio_tot
        
    if (macro_mol.topology.windows == None or 
        any(macro_mol.topology.windows) == 0 or 
        len(macro_mol.topology.windows) < macro_mol.topology.window_num):
        window_diff = 500
    else:    
        if len(macro_mol.topology.windows) > macro_mol.topology.window_num:
            macro_mol.topology.windows = macro_mol.topology.windows[:macro_mol.topology.window_num]

    
        window_diff = 0.0                   # Defining the initial window difference
        for i,j in itertools.combinations(range(len(macro_mol.topology.windows)), 2):
    
            window1 = macro_mol.topology.windows[i]
            window2 = macro_mol.topology.windows[j]
            window_diff += abs(window1 - window2)
        
        if len(macro_mol.topology.windows) > 1: 
            diff_num = 0
            i = len(macro_mol.topology.windows) - 1
        
            while i > 0:
                diff_num += i
                i -= 1
        else:
            diff_num = 1
        
        window_diff = (np.square(window_diff/diff_num)) * 100
    
    macro_mol.topology.window_diff = window_diff
        

    """
    Calculate the cage volume and divide it for the Cavity Size.
    The smaller the CageVolume/CavitySize is the better the structure is.
    The value will be substracted from the individual's val in order to obtain
    an improved fitness_value.
    """

    try:

        if target_size <= 10:
            try:
                val = (abs(target_size - GetCavitySize(macro_mol.topology.mol_file_location)) * 5) ** 2   
            except:
                val = 500
        else:
            try:
                val = abs(target_size - GetCavitySize(macro_mol.topology.mol_file_location)) * 50
            except:
                val = 500
            
        macro_mol.topology.val = val
        
        # Divide the final fitness function by 1000 in order to normalize the value of the fitness function over 1
        fitness_value  = (1000 - val - window_diff - window_ratio_tot) / 1000
        
    except:
        fitness_value = 0.0
    
    return fitness_value    
    
def GetCavitySize(mol_file_name):
    """
    gets the cavity size of a cage
    """
    
    #Atom vdW radii dictionary taken from www.ccdc.cam.ac.uk/Lists/ResourceFileList/Elemental_Radii.xlsx 13 Oct 2015
    #Excluding unstable radioisotopes, dummy atom denoted X (and D) have atomic vdW radii equal 1
    atom_vdw_radii = {
                      'Al': 2,    'Sb': 2,    'Ar': 1.88, 'As': 1.85, 'Ba': 2,    'Be': 2,    'Bi': 2, 
                      'B':  2,    'Br': 1.85, 'Cd': 1.58, 'Cs': 2,    'Ca': 2,    'C':  1.7,  'Ce': 2, 
                      'Cl': 1.75, 'Cr': 2,    'Co': 2,    'Cu': 1.4,  'Dy': 2,    'Er': 2,    'Eu': 2,
                      'F':  1.47, 'Gd': 2,    'Ga': 1.87, 'Ge': 2,    'Au': 1.66, 'Hf': 2,    'He': 1.4,
                      'Ho': 2,    'H':  1.09, 'In': 1.93, 'I':  1.98, 'Ir': 2,    'Fe': 2,    'Kr': 2.02,
                      'La': 2,    'Pb': 2.02, 'Li': 1.82, 'Lu': 2,    'Mg': 1.73, 'Mn': 2,    'Hg': 1.55,
                      'Mo': 2,    'Nd': 2,    'Ne': 1.54, 'Ni': 1.63, 'Nb': 2,    'N':  1.55, 'Os': 2,
                      'O':  1.52, 'Pd': 1.63, 'P':  1.8,  'Pt': 1.72, 'K':  2.75, 'Pr': 2,    'Pa': 2,
                      'Re': 2,    'Rh': 2,    'Rb': 2,    'Ru': 2,    'Sm': 2,    'Sc': 2,    'Se': 1.9,
                      'Si': 2.1,  'Ag': 1.72, 'Na': 2.27, 'Sr': 2,    'S':  1.8,  'Ta': 2,    'Te': 2.06,
                      'Tb': 2,    'Tl': 1.96, 'Th': 2,    'Tm': 2,    'Sn': 2.17, 'Ti': 2,    'W':  2,
                      'U':  1.86, 'V':  2,    'Xe': 2.16, 'Yb': 2,    'Y':  2,    'Zn': 1.29, 'Zr': 2,
                      'X':  1.0,  'D':  1.0
                     }
    
    mol = chem.MolFromMolFile(mol_file_name)
    try:
        conformer = mol.GetConformer(0)
    except:
        return 0
    
    atom_list = [x.GetIdx() for x in mol.GetAtoms()]
    centre_of_mass = GetAtomGroupCoM(mol, conformer, atom_list)
    shortest_distance = 999.9
    
    for atom_id in atom_list:
        
        cur_atom_x = conformer.GetAtomPosition(atom_id).x
        cur_atom_y = conformer.GetAtomPosition(atom_id).y
        cur_atom_z = conformer.GetAtomPosition(atom_id).z

        x_dist_sq = (cur_atom_x - centre_of_mass[0]) ** 2
        y_dist_sq = (cur_atom_y - centre_of_mass[1]) ** 2
        z_dist_sq = (cur_atom_z - centre_of_mass[2]) ** 2
        
        element = str(mol.GetAtomWithIdx(atom_id).GetSymbol())
        
        #Calculates the cage's radius        
        dist = float(np.sqrt(x_dist_sq + y_dist_sq + z_dist_sq)) - atom_vdw_radii[element]
        
        if dist < shortest_distance:                  
            shortest_distance = dist
    
    ##Returns the shortest diameter (radius * 2)
    return (shortest_distance * 2)

    
def GetAtomGroupCoM(molecule, conformer, atom_list):
    """
    gets the centre of mass of a group of atoms
    """
    
    mass_xposition_product_total = 0.0
    mass_yposition_product_total = 0.0
    mass_zposition_product_total = 0.0    
    mass_total = 0.0    
    
    for atom_id in atom_list:    
    
        atom_position = (conformer.GetAtomPosition(atom_id).x, conformer.GetAtomPosition(atom_id).y, 
                         conformer.GetAtomPosition(atom_id).z)
    
        atom_mass =  molecule.GetAtomWithIdx(atom_id).GetMass()
        
        mass_xposition_product_total += atom_mass * atom_position[0]        
        mass_yposition_product_total += atom_mass * atom_position[1]
        mass_zposition_product_total += atom_mass * atom_position[2]

        mass_total += atom_mass        
        
    x_coord = mass_xposition_product_total / mass_total
    y_coord = mass_yposition_product_total / mass_total
    z_coord = mass_zposition_product_total / mass_total
    
    return x_coord, y_coord,z_coord    
    
    
    
    
    