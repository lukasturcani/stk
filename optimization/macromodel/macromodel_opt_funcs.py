import os
import subprocess as sp
import time
import rdkit.Chem as chem
import rdkit.Chem.AllChem as ac
import warnings
import psutil
from multiprocessing import Pool
from functools import partial

# More imports at bottom.

class ConversionError(Exception):
    def __init__(self, message):
        self.message = message
class PathError(Exception):
    def __init__(self, message):
        self.message = message
class LicenseError(Exception):
    def __init__(self, message):
        print('License not found.')
        self.message = message
class ForceFieldError(Exception):
    def __init__(self, message):
        self.message = message
class OptimizationError(Exception):
    def __init__(self, message):
        self.message = message
class ConformerIdentificationError(Exception):
    def __init__(self, message):
        self.message = message
class ConformerExtractionError(Exception):
    def __init__(self, message):
        self.message = message

def macromodel_opt(macro_mol, force_field=16,
                 macromodel_path=r"C:\Program Files\Schrodinger2016-2",
                 no_fix=False, md=False):
    """
    Optimizes the molecule using MacroModel.

    This function runs a restricted optimization. The structures of the
    building blocks are frozen and only the new bonds formed between
    building blocks during assembly are optimized.    
    
    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule who's structure must be optimized.
        
    macromodel_path : str
        The full path of the ``Schrodinger`` suite within the user's 
        machine. For example, in a default Microsoft installation the 
        folder will probably be something like
        ``C:\Program Files\Schrodinger2016-2``.

    force_field : int (default = 16)
        The number of the force field, within MacroModel, to be used in 
        the optimization.
        
    no_fix : bool (default = False)
        When ``True`` the molecular parameters will not be fixed during
        the optimization.
        
    md : bool (default = False)
        If ``True`` then a macromodel optimization using MD is carried
        out to sample different conformtions.
    
    Modifies
    --------
    macro_mol.prist_mol
        The rdkit molecule held in this attribute is replaced by an 
        rdkit molecule with an optimized structure.
    
    macro_mol.prist_mol_file's content
        The content of the ``.mol`` file located at 
        `macro_mol.prist_mol_file`, is changed so that it holds the
        structure of the optimized rdkit molecule.
    
    macro_mol.optimized
        After the optimization, this attribute is set to ``True``.
    
    Returns
    -------
    macro_mol : MacroMolecule
        The macromolecule that was passed as an argument and modified
        by the optimization function. Returned to accomodate 
        parallelization. See ``optimize_all`` function documentation for
        more details.               
    
    """
    
    # If the molecule is already optimized, return.
    if macro_mol.optimized:
        print('Skipping {0}.'.format(macro_mol.prist_mol_file))       
        return macro_mol
    
    print('\nOptimizing {0}.'.format(macro_mol.prist_mol_file))    
    try:
        # MacroModel requires a ``.mae`` file as input. This creates a 
        # ``.mae`` file holding the molding the pristine molecule.    
        convert_mol_to_mae(macro_mol, macromodel_path)        
        # generate the ``.com`` file for the MacroModel run.
        generate_com(macro_mol, force_field, no_fix)        
        # Run the optimization.
        run_bmin(macro_mol, macromodel_path)
        # Get the ``.maegz`` file output from the optimization and 
        # convert it to a ``.mol2`` file.
        convert_maegz_to_mol2(macro_mol, macromodel_path)
        update_prist_attrs_from_mol2(macro_mol) 

        if md:
            macromodel_md_opt(macro_mol, macromodel_path)
            
        macro_mol.optimized = True       
        return macro_mol 

    except ConversionError as ex:
        MacroMolError(ex, macro_mol, '`structconvert` failed.')
        return macro_mol        

    except PathError as ex:
        MacroMolError(ex, macro_mol, 'Wrong MacroModel path.')
        return macro_mol

    except OptimizationError as ex:
        MacroMolError(ex, macro_mol, 'Optimization by `bmin` failed.')
        return macro_mol

    except ForceFieldError as ex:        
        # If OPLS_2005 has been tried already - record an exception.
        if force_field == 14:
            MacroMolError(Exception(), macro_mol, 
                          'Both force fields failed.')
            return macro_mol
            
        # If OPLSE_2005 has not been tried - try it.
        print(('Minimization with OPLS3 failed. Trying OPLS_2005. '
               '- {0}').format(macro_mol.prist_mol_file))
        return macromodel_opt(macro_mol, force_field=14, 
                              macromodel_path=macromodel_path,
                              no_fix=no_fix, md=md)   

    except Exception as ex:
        MacroMolError(ex, macro_mol, ('Uncategorized '
                      'exception during `macromodel_opt`.'))
        return macro_mol
       
def macromodel_md_opt(macro_mol, macromodel_path, 
                      timeout=True, force_field=16, 
                      temp=300, confs=50, eq_time=10, sim_time=200):  

    print('\nRunning MD on {0}.'.format(macro_mol.prist_mol_file))    
    try:
        # MacroModel requires a ``.mae`` file as input. This creates a 
        # ``.mae`` file holding the molding the pristine molecule.    
        convert_mol_to_mae(macro_mol, macromodel_path)        
        # Generate the ``.com`` file for the MacroModel MD run.
        generate_md_com(macro_mol, force_field=force_field, temp=temp, 
                        confs=confs, eq_time=eq_time, sim_time=sim_time)
        # Run the optimization.
        run_bmin(macro_mol, macromodel_path, timeout)        
        # Find the number of the lowest energy conformer.
        conf_num = low_energy_conf(macro_mol)
        # Extract the lowest energy conformer into its own .mae file.
        extract_conformer(macro_mol, conf_num, macromodel_path)
        # Convert the .mae file of the lowest energy conformer to a .mol
        # file.
        convert_mae_to_mol2(macro_mol, macromodel_path)
        update_prist_attrs_from_mol2(macro_mol) 

    except ConversionError as ex:
        MacroMolError(ex, macro_mol, '`structconvert` failed.')
        return macro_mol

    except PathError as ex:
        MacroMolError(ex, macro_mol, 'Wrong MacroModel path.')
        return macro_mol

    except OptimizationError as ex:
        MacroMolError(ex, macro_mol, 'Optimization by `bmin` failed.')
        return macro_mol

    except ConformerIdentificationError as ex:
        MacroMolError(ex, macro_mol, 'Lowest energy conformer not found.')
        return macro_mol

    except ConformerExtractionError as ex:
        MacroMolError(ex, macro_mol, 'Conformer extraction failed.')
        return macro_mol

    except ForceFieldError as ex:        
        # If OPLS_2005 has been tried already - record an exception.
        if force_field == 14:
            MacroMolError(Exception(), macro_mol, 
                          'Both force fields failed.')
            return macro_mol
            
        # If OPLSE_2005 has not been tried - try it.
        print(('Minimization with OPLS3 failed. Trying OPLS_2005. '
               '- {0}').format(macro_mol.prist_mol_file))
        return macromodel_md_opt(macro_mol, macromodel_path, 
                                 timeout=timeout, force_field=14, 
                                 temp=temp, confs=confs, 
                                 eq_time=eq_time, sim_time=sim_time) 
        
    except Exception as ex:
        MacroMolError(ex, macro_mol, ('Uncategorized'
                       ' exception during `macromodel_md_opt`.'))
        return macro_mol
    
def macromodel_cage_opt(macro_mol, force_field=16,
                 macromodel_path=r"C:\Program Files\Schrodinger2016-2",
                 no_fix=False, md=False):
    """
    Optimizes the molecule using MacroModel.

    This function runs a restricted optimization. The structures of the
    building blocks are frozen and only the new bonds formed between
    building blocks during assembly are optimized.    
    
    This function differes from `macromodel_opt` in that it checks the
    number of windows the `macro_mol` has before running the MD. The MD
    is only run all windows are found  (i.e. the cage is not collapsed).    
    
    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule who's structure must be optimized.
        
    macromodel_path : str
        The full path of the ``Schrodinger`` suite within the user's 
        machine. For example, in a default Microsoft installation the 
        folder will probably be something like
        ``C:\Program Files\Schrodinger2016-2``.

    force_field : int (default = 16)
        The number of the force field, within MacroModel, to be used in 
        the optimization.
        
    no_fix : bool (default = False)
        When ``True`` the molecular parameters will not be fixed during
        the optimization.
        
    md : bool (default = False)
        If ``True`` then a macromodel optimization using MD is carried
        out to sample different conformtions.
    
    Modifies
    --------
    macro_mol.prist_mol
        The rdkit molecule held in this attribute is replaced by an 
        rdkit molecule with an optimized structure.
    
    macro_mol.prist_mol_file's content
        The content of the ``.mol`` file located at 
        `macro_mol.prist_mol_file`, is changed so that it holds the
        structure of the optimized rdkit molecule.
    
    macro_mol.optimized
        After a successful optimization, this attribute is set to 
        ``True``.
    
    Returns
    -------
    macro_mol : MacroMolecule
        The macromolecule that was passed as an argument and modified
        by the optimization function. Returned to accomodate 
        parallelization. See ``optimize_all`` function documentation for
        more details.               
    
    """
    
    # If the molecule is already optimized, return.
    if macro_mol.optimized:
        print('Skipping {0}.'.format(macro_mol.prist_mol_file))       
        return macro_mol
    
    print('\nOptimizing {0}.'.format(macro_mol.prist_mol_file))    
    try:    
        # MacroModel requires a ``.mae`` file as input. This creates a 
        # ``.mae`` file holding the molding the pristine molecule.    
        convert_mol_to_mae(macro_mol, macromodel_path)        
        # generate the ``.com`` file for the MacroModel run.
        generate_com(macro_mol, force_field, no_fix)
        # Run the optimization.
        run_bmin(macro_mol, macromodel_path)
        # Get the ``.maegz`` file output from the optimization and 
        # convert it to a ``.mol2`` file.
        convert_maegz_to_mol2(macro_mol, macromodel_path)
        update_prist_attrs_from_mol2(macro_mol) 

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if macro_mol.topology.windows is not None:
                all_windows = (len(macro_mol.topology.windows) == 
                                       macro_mol.topology.n_windows)

                if md and all_windows:
                    macromodel_md_opt(macro_mol, macromodel_path)

        macro_mol.optimized = True       
        return macro_mol

    except ConversionError as ex:
        MacroMolError(ex, macro_mol, '`structconvert` failed.')
        return macro_mol        

    except PathError as ex:
        MacroMolError(ex, macro_mol, 'Wrong MacroModel path.')
        return macro_mol

    except OptimizationError as ex:
        MacroMolError(ex, macro_mol, 'Optimization by `bmin` failed.')
        return macro_mol

    except ForceFieldError as ex:        
        # If OPLS_2005 has been tried already - record an exception.
        if force_field==14:
            MacroMolError(Exception(), macro_mol, 
                          'Both force fields failed.')
            return macro_mol
            
        # If OPLSE_2005 has not been tried - try it.
        print(('Minimization with OPLS3 failed. Trying OPLS_2005. '
       '- {0}').format(macro_mol.prist_mol_file))
        return macromodel_cage_opt(macro_mol, force_field=14, 
                              macromodel_path=macromodel_path,
                              no_fix=no_fix, md=md)    

    except Exception as ex:
        MacroMolError(ex, macro_mol, ('Uncategorized '
                      'exception during `macromodel_cage_opt`.'))
        return macro_mol

def run_bmin(macro_mol, macromodel_path, timeout=True):

    print("", time.ctime(time.time()),
    'Running bmin - {0}.'.format(macro_mol.prist_mol_file), sep='\n')
    
    # To run MacroModel a command is issued to to the console via
    # ``subprocess.Popen``. The command is the full path of the ``bmin``
    # program. ``bmin`` is located in the Schrodinger installation
    # folder.
    file_root = macro_mol.prist_mol_file.replace(".mol", "")
    opt_app = os.path.join(macromodel_path, "bmin")
    # The first member of the list is the command, the following ones
    # are any additional arguments.
    opt_cmd = [opt_app, file_root] 

    opt_proc = psutil.Popen(opt_cmd, stdout=sp.PIPE, 
                            stderr=sp.STDOUT, 
                            universal_newlines=True)
    try:
        if timeout:
            proc_out, _ = opt_proc.communicate(timeout=600)
        else:
            proc_out, _ = opt_proc.communicate()   
        opt_proc.wait()
    except sp.TimeoutExpired:
        print(('\nMinimization took too long and was terminated '
               'by force - {}\n').format(macro_mol.prist_mol_file))
        kill_bmin()
        proc_out = ""

    # If optimization fails because a wrong Schrodinger path was given,
    # raise.
    if 'The system cannot find the path specified' in proc_out:
        raise PathError(('Wrong Schrodinger path supplied to'
                              ' `macromodel_opt` function.'))

    # If optimization fails because the license is not found, rerun the
    # function.
    if not license_found(macro_mol, proc_out):
        return run_bmin(macro_mol, macromodel_path)

    # Make sure the .log and .maegz files which should be created by
    # the optimization are present.
    log_file = macro_mol.prist_mol_file.replace('.mol', '.log')
    wait_for_file(log_file)
    maegz = macro_mol.prist_mol_file.replace('.mol', '-out.maegz')
    wait_for_file(maegz)
    if not os.path.exists(log_file) or not os.path.exists(maegz):
        raise OptimizationError(('The .log and/or .maegz '
                     'files were not created by the optimization.'))   

def kill_bmin():
    bmins = []
    for i, process in enumerate(psutil.process_iter()):
        if process.name() == 'bmin.exe':
            process_time = (process.cpu_times().user + 
                            process.cpu_times().system)
    
            bmins.append((process_time, i, process))

    bmins.sort(reverse=True)
    for x in bmins:
        try:
            x[-1].kill()
            break
        except:
            print('``kill_bmin`` excepted.')
            continue

def license_found(macro_mol, output):
    """
    Checks to see if minimization failed due to a missing license.

    The user can be notified of this in one of two ways. Sometimes the
    output of the submission contains the message informing that the 
    license was not found and in other cases it will be the log file.
    This function checks both of these sources for this message.

    Parameters
    ----------
    macro_mol : MacroMolecule or any other type
        The macromolecule being optimized. If the .log file is not to
        be checked, a non MacroMolecule object should be placed here
        instead. An empty string would do nicely.
    
    output : str
        The outout from submitting the minimization of the structure
        to the ``bmin`` program.
        
    Returns
    -------
    bool
        ``True`` if the license was found. ``False`` if the minimization
        did not occur due to a missing license.
    
    """

    if 'Could not check out a license for mmlibs' in output:
        return False
    if not isinstance(macro_mol, MacroMolecule):
        return True

    # To check if the log file mentions a missing license file open the
    # the log file and scan for the apporpriate string.
    
    # Check if the file exists first. If not, this is often means the
    # calculation must be redone so return False anyway.
    log_file_path = macro_mol.prist_mol_file.replace('mol', 'log')
    if not os.path.exists(log_file_path):    
        return False
    with open(log_file_path, 'r') as log_file:
        log_file_content = log_file.read()
        
    if 'Could not check out a license for mmlibs' in log_file_content:
        return False

    return True
 
def generate_com(macro_mol, force_field=16, no_fix=False):
    """
    Create a ``.com`` file for a MacroModel optimization.

    The created ``.com`` file fixes all bond parameters which were not
    added during assembly. This means all bond distances, bond angles
    and torsional angles are fixed, except for cases where it involves
    a bond added during assembly of the macromolecule.
    
    This fixing is implemented by creating a ``.com`` file with various
    ``FX`` commands written within its body.
    
    This function is called by ``macromodel_opt``. It is private because
    it should probably not be used outside of this context.

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule which is to be optimized.
        
    force_field : int (default = 16)
        The number of the force field, within MacroModel, to be used in 
        the optimization.
        
    no_fix : bool (default = False)
        When ``True`` the generated .com file will not contain commands
        which fix the molecular parameters during optimization.

    Modifies
    --------
    This function creates a new ``.com`` file holding the instructions
    for optimizing the pristine macromolecule using MacroModel.

    Returns
    -------
    None : NoneType    
    
    """

    print('Creating .com file - {0}.'.format(macro_mol.prist_mol_file))
    
    # This is the body of the ``.com`` file. The line that begins and
    # ends with exclamation lines is replaced with the various commands
    # that fix bond distances and angles.
    main_string= (" MMOD       0      1      0      0     0.0000     "
    "0.0000     0.0000     0.0000\n"
" DEBG      55      0      0      0     0.0000     0.0000     "
"0.0000     0.0000\n"
" FFLD{0:8}      1      0      0     1.0000     0.0000     "
"0.0000     0.0000\n"
" BDCO       0      0      0      0    41.5692 99999.0000     "
"0.0000     0.0000\n"
" CRMS       0      0      0      0     0.0000     0.5000     "
"0.0000     0.0000\n"
" BGIN       0      0      0      0     0.0000     0.0000     "
"0.0000     0.0000\n"
" READ       0      0      0      0     0.0000     0.0000     "
"0.0000     0.0000\n"
"!!!BLOCK_OF_FIXED_PARAMETERS_COMES_HERE!!!\n"
" CONV       2      0      0      0     0.0500     0.0000     "
"0.0000     0.0000\n"
" MINI       1      0   2500      0     0.0000     0.0000     "
"0.0000     0.0000\n"
" END        0      0      0      0     0.0000     0.0000     "
"0.0000     0.0000\n"
" ").format(force_field)

    # Create a path for the ``.com`` file. It is the same as that of the
    # ``.mol`` file but with a ``.com`` extension. Get the path of the
    # ``.mae`` file and the output file in the same way. 
    com_file = macro_mol.prist_mol_file.replace(".mol", ".com")
    mae = macro_mol.prist_mol_file.replace(".mol", ".mae")
    output = macro_mol.prist_mol_file.replace(".mol", "-out.maegz")
    
    # This function adds all the lines which fix bond distances and 
    # angles into ``main_string``.
    main_string = fix_params_in_com_file(macro_mol, 
                                          main_string, no_fix)
    
    # Writes the ``.com`` file.
    with open(com_file, "w") as com:
        # The first line hold the ``.mae`` file containing the molecule
        # to be optimized.
        com.write(str(mae + "\n"))
        # The second line holds the name of the output file of the 
        # optimization.
        com.write(str(output + "\n"))
        # Next is the body of the ``.com`` file, held in
        # ``main_string``.
        com.write(main_string)

def generate_md_com(macro_mol, force_field=16, temp=300, confs=50, eq_time=10, sim_time=200):

    print('Creating .com file - {0}.'.format(macro_mol.prist_mol_file))

    # Defining the string to be printed in the COM file - uses OPLS3 (FFLD = 16)
    # run a 200 ns MD, at 300K and optimize 50 random conformations generated during the trajectory
    
    main_string= """ MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000
 FFLD{force_field:8}      1      0      0     1.0000     0.0000     0.0000     0.0000
 BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000
 READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000
 CONV       2      0      0      0     0.0500     0.0000     0.0000     0.0000
 MINI       1      0   2500      0     0.0000     0.0000     0.0000     0.0000
 MDIT       0      0      0      0   300.0000     0.0000     0.0000     0.0000
 MDYN       0      0      0      0     1.5000{eq_time:6}.0000{temp:6}.0000     0.0000
 MDSA{confs:8}      0      0      0     0.0000     0.0000     1.0000     0.0000
 MDYN       1      0      0      0     1.5000{sim_time:6}.0000{temp:6}.0000     0.0000
 WRIT       0      0      0      0     0.0000     0.0000     0.0000     0.0000
 RWND       0      1      0      0     0.0000     0.0000     0.0000     0.0000
 BGIN       0      0      0      0     0.0000     0.0000     0.0000     0.0000
 READ      -2      0      0      0     0.0000     0.0000     0.0000     0.0000
 CONV       2      0      0      0     0.0500     0.0000     0.0000     0.0000
 MINI       1      0   2500      0     0.0000     0.0000     0.0000     0.0000
 END        0      0      0      0     0.0000     0.0000     0.0000     0.0000"""

    main_string = main_string.format(force_field=force_field, temp=temp, 
                        confs=confs, eq_time=eq_time, sim_time=sim_time)

    com_file = macro_mol.prist_mol_file.replace(".mol", ".com")
    mae = macro_mol.prist_mol_file.replace(".mol", ".mae")
    output = macro_mol.prist_mol_file.replace(".mol", "-out.maegz")

    # Generate the com file containing the info for the run
    with open(com_file, "w") as com:
        # name of the macromodel file
        com.write(str(mae + "\n"))
        # name of the output file
        com.write(str(output + "\n"))
        # details of the macromodel run
        com.write(main_string)

def convert_mol_to_mae(macro_mol, macromodel_path):
    """
    Creates the ``.mae`` file holding the molecule to be optimized.    

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule which is to be optimized. Its ``.mol`` file is
        converted to a ``.mae`` file. The original ``.mol`` file is also
        kept.
        
    macromodel_path : str
        The full path of the installation directory of the Schrodinger 
        suite. By default on a Windows machine it should be something
        like: "C:\Program Files\Schrodinger2016-2".

    Modifies
    --------
    This function creates a new ``.mae`` file from the ``.mol`` file in
    `macro_mol.prist_mol_file`. This new file is placed in the same
    folder as the ``.mol`` file and has the same name. Only the 
    extensions are different.

    Returns
    -------
    str
        The full path of the newly created ``.mae`` file.     
    
    """
 
    print('Converting .mol to .mae - {}.'.format(
                                              macro_mol.prist_mol_file))
   
    # Create the name of the new ``.mae`` file. It is the same as the
    # ``.mol`` file, including the same path. Only the extensions are
    # different.
    mae_file = macro_mol.prist_mol_file.replace('.mol', '.mae')  
    structconvert(macro_mol.prist_mol_file, mae_file, macromodel_path)
    return mae_file

def convert_maegz_to_mol2(macro_mol, macromodel_path):
    """
    Converts a ``.mae`` file to a ``.mol2`` file.
    
    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule being optimized. The ``.mae`` file holding its
        optimized structure is converted to a ``.mol2`` file. Both
        versions are kept.
        
    macromodel_path : str
        The full path of the installation directory of the Schrodinger 
        suite. By default on a Windows machine it should be something
        like: "C:\Program Files\Schrodinger2016-2".   
    
    Modifies
    --------    
    This function creates a new ``.mol2`` file from the optimized 
    ``.mae`` file. This new file is placed in the same folder as the 
    ``.mae`` file.
    
    Returns
    -------
    None : NoneType

    Raises
    ------
    ForceFieldError
        If the OPLS3 force field failed to optimize the molecule. If
        this happens the conversion function is unable to convert the 
        output of the optimization function and as a result this error
        is raised.
    
    """

    print('Converting .maegz to .mol2 - {}.'.format(
                                            macro_mol.prist_mol_file))

    # ``out`` is the full path of the optimized ``.mae`` file.
    maegz = macro_mol.prist_mol_file.replace(".mol", "-out.maegz")      
    # Replace extensions to get the names of the various files.
    mol2 = macro_mol.prist_mol_file.replace(".mol", ".mol2")
    structconvert(maegz, mol2, macromodel_path)

def convert_mae_to_mol2(macro_mol, macromodel_path):
    """
    Converts a ``.mae`` file to a ``.mol2`` file.
    
    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule being optimized. The ``.mae`` file holding its
        optimized structure is converted to a ``.mol2`` file. Both
        versions are kept.
        
    macromodel_path : str
        The full path of the installation directory of the Schrodinger 
        suite. By default on a Windows machine it should be something
        like: "C:\Program Files\Schrodinger2016-2".   
    
    Modifies
    --------    
    This function creates a new ``.mol2`` file from the optimized 
    ``.mae`` file. This new file is placed in the same folder as the 
    ``.mae`` file.
    
    Returns
    -------
    None : NoneType

    Raises
    ------
    ForceFieldError
        If the OPLS3 force field failed to optimize the molecule. If
        this happens the conversion function is unable to convert the 
        output of the optimization function and as a result this error
        is raised.
    
    """

    print('Converting .mae to .mol2 - {}.'.format(
                                            macro_mol.prist_mol_file))
    
    # ``mae`` is the full path of the optimized ``.mae`` file.
    mae = macro_mol.prist_mol_file.replace(".mol", ".mae")        
    # Replace extensions to get the names of the various files.
    mol2 = macro_mol.prist_mol_file.replace(".mol", ".mol2")
    structconvert(mae, mol2, macromodel_path)    
    
def structconvert(iname, oname, macromodel_path):
  
    convrt_app = os.path.join(macromodel_path, 'utilities', 
                                                     'structconvert')
    convrt_cmd = [convrt_app, iname, oname]

    # Execute the file conversion.
    convrt_return = sp.run(convrt_cmd, stdout=sp.PIPE, 
                           stderr=sp.STDOUT, universal_newlines=True) 

    # If conversion fails because a wrong Schrodinger path was given,
    # raise.
    if 'The system cannot find the path specif' in convrt_return.stdout:
        raise PathError(('Wrong Schrodinger path supplied to'
                              ' `structconvert` function.'))

    # If no license if found, keep re-running the function until it is.
    if not license_found('', convrt_return.stdout):
        return structconvert(iname, oname, macromodel_path)    

    # If force field failed, raise.
    if 'number 1' in convrt_return.stdout:
        raise ForceFieldError(convrt_return.stdout)    

    wait_for_file(oname)
    if not os.path.exists(oname):
        raise ConversionError(('Conversion output file {} was not found.'
                ' Console output was {}.').format(oname, 
                                                  convrt_return.stdout))

    return convrt_return

def fix_params_in_com_file(macro_mol, main_string, no_fix=False):
    """
    Adds lines to the ``.com`` body fixing bond distances and angles.
    
    For each bond distance, bond angle and torisional angle that does
    not involve a bond created during assembly a ``FX`` command is added
    to the string holding holding the body of the ``.com`` file.
    
    These lines replace the filler line in the main string.
    
    This function is called by ``macromodel_opt``. It is private because
    it should probably not be used outside of this context.    
    
    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule which is to be optimized.
        
    main_string : str
        The body of the ``.com`` file which is to have fix commands
        added.
    
    no_fix : bool (default = False)
        When ``True`` the block containing instructions to fix molecular
        parameters is not added to the .com file.
        
    Returns
    -------
    str
        A string holding the body of the ``.com`` file with instructions
        to fix the various bond distances and angles as described in the
        docstring.
    
    """
    
    # Make a string to hold all of the ``FX`` lines.
    fix_block = "" 

    # If no_fix is ``True`` do not add a fix block.
    if no_fix:
        return main_string.replace(("!!!BLOCK_OF_FIXED_PARAMETERS_"
                                    "COMES_HERE!!!\n"), fix_block)
    # Add lines that fix the bond distance.
    fix_block = fix_distance_in_com_file(macro_mol, fix_block)  
    # Add lines that fix the bond angles.                          
    fix_block = fix_bond_angle_in_com_file(macro_mol, fix_block)
    # Add lines that fix the torsional angles.
    fix_block = fix_torsional_angle_in_com_file(macro_mol, fix_block)
    
    return main_string.replace(("!!!BLOCK_OF_FIXED_PARAMETERS_"
                                "COMES_HERE!!!\n"), fix_block)

def fix_distance_in_com_file(macro_mol, fix_block):
    """
    Adds lines fixing bond distances to ``.com`` body string.
    
    Only bond distances which do not involve bonds created during
    assembly are fixed.

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule to be optimized.
        
    fix_block : str
        The string holding all the lines containing fix commands for the
        ``.com`` file.
        
    Returns
    -------
    str
        A string holding lines containg fix commands for the ``.com`` 
        file. The string has the lines fixing bond distances added to it
        by this function.
    
    """
    
    # This line holds the format for a line fixing the bond distance
    # between two atoms. The first two ``{...}`` are replaced with the
    # ids of atoms to be fixed. The last ``{...}`` is replaced with the
    # bond distance. Note that in the ``.mae`` files the indices of
    # atoms start at 1 while in rdkit they start at 0. As far as I can 
    # tell this corresponds to a shift of one for each atom index, with
    # the ordering being the same.
    fix_distance = (" FXDI {0:>7}{1:>7}      0      0"
                    "   100.0000 {2:>10.4f}     0.0000     0.0000")
   
    # Go through all the bonds in the heavy rdkit molecule. If the bond
    # is not between heavy atoms get its distance. Add a fix line using
    # the bond distance and atomic indices to the ``fix_block``. If the
    # bond does invovle two heavy atoms go to the next bond. This is
    # because a bond between 2 heavy atoms was added during assembly and
    # should therefore not be fixed.  
    for bond in macro_mol.heavy_mol.GetBonds():
        atom1 = bond.GetBeginAtom() 
        atom2 = bond.GetEndAtom()
        
        if (atom1.GetAtomicNum() in FGInfo.heavy_atomic_nums and
            atom2.GetAtomicNum() in FGInfo.heavy_atomic_nums):
            continue
        
        atom1_id = atom1.GetIdx() 
        atom2_id = atom2.GetIdx() 
        
        bond_len = macro_mol.atom_distance('prist', atom1_id, atom2_id)
        
        # Make sure that the indices are increased by 1 in the ``.mae``
        # file from their rdkit value.
        fix_block += (fix_distance.format(atom1_id+1, atom2_id+1,
                                         bond_len) + "\n")

    return fix_block
    
def fix_bond_angle_in_com_file(macro_mol, fix_block):
    """
    Adds lines fixing bond angles to ``.com`` body string.
    
    Only bond angles which do not involve bonds created during
    assembly are fixed. The exception to this is for bond angles next
    to bonds added during assembly. For example, consider:
        
        A-B-C-D=E
        
    if the bond between D and E (``=``) represents the bond added during
    assembly, the bond angle A-B-C will be fixed but B-C-D will not be.
    This is an artifact of the implementation but is not expected to
    play a significant role as the vast majority of bond angles which
    should be fixed, will be. The bond angle C-D=E will also not be
    fixed as that is the purpose of this function.         

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule to be optimized.
        
    fix_block : str
        The string holding all the lines containing fix commands for the
        ``.com`` file.
        
    Returns
    -------
    str
        A string holding lines containg fix commands for the ``.com`` 
        file. The string has the lines fixing bond angles added to it by
        this function.
    
    """
    # This line holds the format for a line fixing the bond angles
    # between 3 atoms. The first 3 ``{...}`` are replaced with the ids
    # of atoms. The last ``{...}`` is replaced with the bond angle. Note 
    # that in the ``.mae`` files the indices of atoms start at 1 while 
    # in rdkit they start at 0. As far as I can tell this corresponds to 
    # a shift of 1 for each atom index, with the ordering being the 
    # same.    
    fix_ba = (" FXBA {0:>7}{1:>7}{2:>7}      0   100.0000 "
              "{3:>10.4f}     0.0000     0.0000")
    
    # Create a substructure consisting of 3 dummy atoms bonded with 3
    # dummy bonds. This substructure will match with any 3 atoms which
    # are bonded together with any combination of bonds. These 3 atoms
    # will therefore have a bond angle.          
    ba_mol = chem.MolFromSmarts('[*]~[*]~[*]')
    
    # Get the indices of all atoms which have a bond angle. ``ba_atoms``
    # is a tuple of tuples of the form ((1,2,3), (4,5,6), (7,8,9), ...).
    # Each inner tuple holds the indicies of the atoms which form a bond
    # angle.
    ba_atoms = macro_mol.heavy_mol.GetSubstructMatches(ba_mol)
    
    # Get the conformer holding the atomic positions.
    conf = macro_mol.heavy_mol.GetConformer()    
    
    # For each bond angle check if a heavy atom is involved in forming
    # it. If no, a line fixing the bond angle is added to ``fix_block``.
    # If any atom of the 3 is a heavy atom the bond angle is not fixed.
    # This means that there will be some bond angles which consist of 2
    # bonds not added during assembly which will not be fixed. However,
    # it is assumed that the effect of this will be minimal.
    for atom1_id, atom2_id, atom3_id in ba_atoms:
        atom1 = macro_mol.heavy_mol.GetAtomWithIdx(atom1_id)
        atom2 = macro_mol.heavy_mol.GetAtomWithIdx(atom2_id)
        atom3 = macro_mol.heavy_mol.GetAtomWithIdx(atom3_id)
        
        if (atom1.GetAtomicNum() in FGInfo.heavy_atomic_nums or
            atom2.GetAtomicNum() in FGInfo.heavy_atomic_nums or
            atom3.GetAtomicNum() in FGInfo.heavy_atomic_nums):
            continue
        
        ba = ac.GetAngleDeg(conf, atom1_id, atom2_id, atom3_id)
        
        fix_block += (fix_ba.format(atom1_id+1, atom2_id+1, 
                                    atom3_id+1, ba) + "\n")
    
    return fix_block
    
def fix_torsional_angle_in_com_file(macro_mol, fix_block):
    """
    Adds lines fixing torsional bond angles to ``.com`` body string.
    
    Only torsional angles which do not involve bonds created during
    assembly are fixed. The exception to this is for torsional angles 
    next to bonds added during assembly. For example, consider:
        
        A-B-C-D-E=F
        
    if the bond between E and F (``=``) represents the bond added during
    assembly, the torsional angle A-B-C-D will be fixed but B-C-D-E will 
    not be. This is an artifact of the implementation but is not 
    expected to play a significant role as the vast majority of bond 
    angles which should be fixed, will be. The bond angle C-D-E=F will 
    also not be fixed as that is the purpose of this function.

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule to be optimized.
        
    fix_block : str
        The string holding all the lines containing fix commands for the
        ``.com`` file.
        
    Returns
    -------
    str
        A string holding lines containg fix commands for the ``.com`` 
        file. The string has the lines fixing bond angles added to it by
        this function.
    
    """

    # This line holds the format for a line fixing the torsional angles
    # between 4 atoms. The first 4 ``{...}`` are replaced with the ids
    # of atoms. The last ``{...}`` is replaced with the torsional angle. 
    # Note that in the ``.mae`` files the indices of atoms start at 1 
    # while in rdkit they start at 0. As far as I can tell this 
    # corresponds to a shift of 1 for each atom index, with the ordering 
    # being the same.    
    fix_ta = (" FXTA {0:>7}{1:>7}{2:>7}{3:>7}   100.0000 "
                "{4:>10.4f}     0.0000     0.0000")        

    # Create a substructure consisting of 4 dummy atoms bonded with 4
    # dummy bonds. This substructure will match with any 4 atoms which
    # are bonded together with any combination of bonds. These 4 atoms
    # will therefore have a torsinal angle.              
    ta_mol = chem.MolFromSmarts('[*]~[*]~[*]~[*]')
    
    # Get the indices of all atoms which have a torsional angle. 
    # ``ta_atoms`` is a tuple of tuples of the form ((1,2,3,4), 
    # (4,5,6,7), ...). Each inner tuple holds the indicies of the atoms 
    # which form a torsional angle.
    ta_atoms = macro_mol.heavy_mol.GetSubstructMatches(ta_mol)
    # Get the conformer holding the atomic positions.
    conf = macro_mol.heavy_mol.GetConformer()
    
    # For each torsional angle check if a heavy atom is involved in 
    # forming it. If no, a line fixing the torsional angle is added to 
    # ``fix_block``. If any atom of the 4 is a heavy atom the bond angle 
    # is not fixed. This means that there will be some bond angles which 
    # consist of 3 bonds not added during assembly which will not be 
    # fixed. However, it is assumed that the effect of this will be 
    # minimal.    
    for atom1_id, atom2_id, atom3_id, atom4_id in ta_atoms:
        atom1 = macro_mol.heavy_mol.GetAtomWithIdx(atom1_id)
        atom2 = macro_mol.heavy_mol.GetAtomWithIdx(atom2_id)
        atom3 = macro_mol.heavy_mol.GetAtomWithIdx(atom3_id)
        atom4 = macro_mol.heavy_mol.GetAtomWithIdx(atom4_id)
        
        if (atom1.GetAtomicNum() in FGInfo.heavy_atomic_nums or
            atom2.GetAtomicNum() in FGInfo.heavy_atomic_nums or
            atom3.GetAtomicNum() in FGInfo.heavy_atomic_nums or
            atom4.GetAtomicNum() in FGInfo.heavy_atomic_nums):
            continue
        
        ta = ac.GetDihedralDeg(conf, atom1_id, atom2_id, 
                                     atom3_id, atom4_id)
        
        fix_block += (fix_ta.format(atom1_id+1, atom2_id+1, 
                                atom3_id+1, atom4_id+1, ta) + "\n")

    return fix_block

def low_energy_conf(macro_mol):
    """
    Return the lowest energy conformer's number.
    
    An MD conformer seach will produce a .log file which holds the
    energies of the various conformers found. This function extracts
    the number of the lowest energy conformer from the .log file.
    
    Parameters
    ----------
    macro_mol : MacroMolecule
        The macro_mol which had an MD conformer search run on it.
        
    Returns
    -------
    int
        The number of the lowest energy conformer.

    """

    log_name = macro_mol.prist_mol_file.replace('.mol', '.log')    
    
    with open(log_name, 'r') as log:
        conformers = []
        for line in log:
            if "Conf" in line and "kJ/mol" in line and "****" not in line:
                conf_num = int(line.split()[1])
                conf_en = float(line.split()[4])
                conformers.append((conf_en, conf_num))

    if not conformers:
        raise ConformerIdentificationError(('No conformers'
                                            ' found in .log file.'))
                                            
    return min(conformers)[1] - 1

def extract_conformer(macro_mol, conf_num, macromodel_path):
    """
    Creates an individual .mae file for a conformer in a .maegz file.
    
    An MD conformer search produces a .maegz file holding the all the
    found conformers together. This function places the conformer with
    the number `conf_num` in its own .mae file. The .mae file reflects
    the name of `macro_mol`.
    
    Parameters
    ----------
    macro_mol : MacroMolecule
        The mocromolecule which had an MD conformer search run on it.
        
    conf_num : int
        The number of the conformer to be extracted. Likely the lowest
        energy conformer.
        
    macromodel_path : str
        The full path of the Schrodinger directory.
        
    Modifies
    --------
    This function creates a .mae file which has the name of 
    `macro_mol.prist_mol_file` except with .mae instead of .mol at the
    end. The .mae file holds the conformer of macro_mol held in a
    .maegz file. The .maegz file was likely produced during an MD
    conformer search.

    Returns
    -------
    None : NoneType
    
    """
    
    print('Extracting conformer - {}.'.format(macro_mol.prist_mol_file)) 
    
    # If the conformer number is 0, the original structure had the
    # lowest energy and no extraction needs to take place.
    if conf_num == 0:
        return
    
    # The names of the input and output files.
    maegz = macro_mol.prist_mol_file.replace('.mol', '-out.maegz')
    mae =  macro_mol.prist_mol_file.replace('.mol', '.mae' )  
    
    # The full path of the application doing the extraction.
    extract_app = os.path.join(macromodel_path, 'utilities', 'maesubset')

    # The command needed to run the application and extract the
    # conformer via bash or whatever.             
    extract_cmd = [extract_app, maegz, "-n", str(conf_num), "-o", mae]
  
    # Execute the extraction.
    extract_return = sp.run(extract_cmd, stdout=sp.PIPE, 
                           stderr=sp.STDOUT, universal_newlines=True)

    # If no license if found, keep re-running the function until it is.
    if not license_found('', extract_return.stdout):   
        return extract_conformer(macro_mol, conf_num, macromodel_path)        

    # Make sure the .mae file is in fact created.
    wait_for_file(mae)
    if not os.path.exists(mae):
        raise ConformerExtractionError(('Conformer extraction failed.'
        ' Console output was {}.').format(extract_return))

def wait_for_file(file_name, timeout=20):
    """
    Stalls until a given file exists or `timeout` expires.
    
    Parameters
    ----------
    file_name : str
        The full path of the file which should be waited for.
        
    timeout : int or float
        The number of seconds before the function stops waiting and
        returns.
        
    Returns
    --------
    None : NoneType
    
    """
    
    t_start = time.time()
    tick = 0
    while True:
        time_taken = time.time() - t_start
        if divmod(time_taken, 5)[0] == tick + 1:
            print('Waiting for {0}.'.format(file_name))
            tick += 1
        
        if os.path.exists(file_name) or time_taken > timeout:
            break 

def kill_macromodel():
    """
    Kills any applications left open as a result running MacroModel.    
    
    Applications that are typically left open are 
    ``jserver-watcher.exe`` and ``jservergo.exe``.    
    
    Returns
    -------
    None : NoneType    
    
    """
    
    if os.name == 'nt':
        # In Windows, use the ``Taskkill`` command to force a close on
        # the applications.           
        sp.run(["Taskkill", "/IM", "jserver-watcher.exe", "/F"])
        sp.run(["Taskkill", "/IM", "jservergo.exe", "/F"])

def optimize_folder(path, macromodel_path):
    
    # First make a list holding all macromolecule objects to be 
    # optimzied. Because the objects need to be initialized from a .mol
    # file a StructUnit instance not a MacroMolecule instance is used.
    # minimal = True, because it is all that is needed to run 
    # optimziations, plus you want to avoid doing functional group
    # substitutions.
    names = [os.path.join(path, file_name) for file_name in 
             os.listdir(path) if file_name.endswith(".mol")]    
    macro_mols = [StructUnit(file_path, minimal=True) for file_path in 
                                                                  names]

    md_opt = partial(macromodel_md_opt, macromodel_path=macromodel_path, 
        timeout=False, temp=1000, confs=500, eq_time=100, sim_time=5000)    
    
    with Pool() as p:
        p.map(md_opt, macro_mols)
    
    

from ...classes.exception import MacroMolError       
from ...classes.molecular import FGInfo, StructUnit, MacroMolecule        
from ..optimization import update_prist_attrs_from_mol2        