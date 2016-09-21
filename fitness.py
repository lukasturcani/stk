import numpy as np
import rdkit.Chem as chem

from .classes.exception import MacroMolError
from .classes.molecular import MacroMolecule, StructUnit
from .optimization import macromodel_opt

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
            
        except Exception as ex:
            MacroMolError(ex, macro_mol, 'During fitness calculation.')
            macro_mol.fitness = -np.inf
            macro_mol.topology.windows = None

        print(macro_mol.fitness, '-', macro_mol.prist_mol_file)            
            
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

def cage(macro_mol, target_size, coeffs=None, exponents=None):
    if macro_mol.fitness:
        print('Skipping {0}'.format(macro_mol.prist_mol_file))
        return macro_mol.fitness
    
    if macro_mol.topology.windows is None:
        return -np.inf

    if coeffs is None:
        coeffs = np.array([50,1,1])
        
    if exponents is None:
        exponents = np.array([1,1,1])

    if target_size <= 10:
        coeffs[0] = 5
        exponents[0] = 2        

    cavity_diff = abs(target_size - macro_mol.topology.cavity_size())

    target_window_area = np.square(target_size)
    window_area = np.square(max(macro_mol.topology.windows))
    window_area_diff = abs(target_window_area - window_area)
            
    fitness_value = np.array([
                             cavity_diff, 
                             window_area_diff,                                                          
                             macro_mol.topology.window_difference(500)
                             ])
    
    fitness_value = np.power(fitness_value, exponents)
    fitness_value = np.multiply(fitness_value, coeffs)    

    return 1/np.sum(fitness_value)
        
class MolInCage:
    
    def __call__(self, cage, target_mol_file, rotate=False):
        if cage.fitness:
            print('Skipping {0}'.format(cage.prist_mol_file))
            return cage.fitness            
        
        if not hasattr(self, 'target_mol'):
            target_mol = StructUnit(target_mol_file)
            macromodel_opt(target_mol)
            setattr(self, 'target_mol', target_mol)
        
        
        # This function creates a new molecule holding both the target
        # and the cage centered at the origin. It then calculates the 
        # energy of this complex and compares it to the energies of the
        # molecules when separate. The more stable the complex relative
        # to the individuals the higher the fitness.
        
        # Create rdkit instances of the target in the cage for each
        # rotation.        
        rdkit_complexes = self.generate_complexes(cage, rotate)
        
        # Place the rdkit instance(s) into a new .mol file and use that 
        # to make a new ``StructUnit`` instance of the complex.
        # ``StructUnit`` is the class of choice here because it can be
        # initialized from a .mol file. ``MacroMolecule`` requires
        # things like topology and building blocks for initialization.
        macromol_complexes = []        
        for i, complex_ in enumerate(rdkit_complexes):
            # First the data is loaded into a ``MacroMolecule`` instance
            # as this is the class which is able to write rdkit
            # instances to .mol files.
            mm_complex = MacroMolecule.__new__(MacroMolecule)
            mm_complex.prist_mol = complex_
            mm_complex.prist_mol_file = target_mol_file.replace('.mol', 
                                           '_COMPLEX_{0}.mol'.format(i))
            mm_complex.write_mol_file('prist')
            
            macromol_complex = StructUnit(mm_complex.prist_mol_file)
            macromodel_opt(macromol_complex)
            macromol_complexes.append(macromol_complex)
        
        # Calculate the energy of the complex and compare to the
        # individual energies. If more than complex was made, use the
        # most stable version.
        energy_separate = cage.energy + self.target_mol.energy
        
        energy_diff =  min(energy_separate - macromol_complex.energy for 
                                macromol_complex in macromol_complexes)
                                
        return np.exp(energy_diff)
   
    def generate_complexes(self, cage, rotate):
        
        # Define the rotations which are to be used on the target.
        rotations = [0, np.pi/2, np.pi, 3*np.pi/2]


        # First place both the target and cage at the origin.
        cage.set_prist_centroid([0,0,0])
        self.target_mol.set_prist_centroid([0,0,0])
        
        

        
molecule_in_cage = MolInCage()


class molecule_in_cage:
    """
    Fitness function for finding cages for holding `target_mol`.
    
    Parameters
    ----------
    cage : Cage
        The cage whose fitness must be calculated.
        
    target_mol_file : str
        The full path of the .mol file holding the target molecule which 
        is to be placed inside the cage.
    
    Returns
    -------
    float
        The fitness value of the cage.
    
    """
    
    # This function is implemented a closure. The idea is that the user
    # will only input the .mol file location into the MMEA input script.
    # From this .mol file a ``MacroMolecule`` instance will be generated 
    # of the target molecule. However it is obviously a complete waste 
    # of resources if the target ``MacroMolecule`` is generated each
    # time the fitness function is called. A closure allows the creation
    # of it only once - the first time it is called.

    
    
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    