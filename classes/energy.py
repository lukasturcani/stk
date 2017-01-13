import os
import rdkit.Chem as chem
import rdkit.Chem.AllChem as ac
import subprocess as sp
import shutil
import uuid

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
        
        # To prevent conflicts when running this function in parallel,
        # a temporary copy of the molecular structure file is made and
        # used for macromodel calculations.
        
        # Unique file name is generated by inserting a random int into 
        # the prist_mol_file path.
        tmp_file = os.path.split(self.molecule.prist_mol_file)[-1]
        tmp_file, ext = os.path.splitext(tmp_file)
        r_int = uuid.uuid4().int
        tmp_file = tmp_file + str(r_int) + ext
        tmp_file = os.path.join(os.getcwd(), tmp_file)
        
        # Create the file.
        shutil.copy(self.molecule.prist_mol_file, tmp_file)
        
        file_root, ext = os.path.splitext(tmp_file)
        convrt_app = os.path.join(macromodel_path, 'utilities', 
                                                    'structconvert')
        convrt_cmd = [convrt_app, 
                     tmp_file, file_root+'.mae']
        sp.call(convrt_cmd, stdout=sp.PIPE, stderr=sp.PIPE)

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
        
        # Clean up temporary files.
        for filename in os.listdir(os.path.split(tmp_file)[0]):
            if str(r_int) in filename:
                os.remove(filename)
        
        return eng