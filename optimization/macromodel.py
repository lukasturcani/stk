"""
Defines optimization functions which use MacroModel.

"""

import os
import subprocess as sp
import time
import rdkit.Chem as chem
import rdkit.Chem.AllChem as ac
import warnings
import psutil
import re

from ..classes.exception import MolError
from ..convenience_tools import MAEExtractor

class _ConversionError(Exception):
    def __init__(self, message):
        self.message = message
class _PathError(Exception):
    def __init__(self, message):
        self.message = message
class _LicenseError(Exception):
    def __init__(self, message):
        print('License not found.')
        self.message = message
class _ForceFieldError(Exception):
    def __init__(self, message):
        self.message = message
class _OptimizationError(Exception):
    def __init__(self, message):
        self.message = message
class _LewisStructureError(Exception):
    def __init__(self, message):
        self.message = message

def macromodel_opt(macro_mol, force_field=16,
                 macromodel_path=r"C:\Program Files\Schrodinger2016-2",
                 no_fix=False, md=False, lewis_fixed=False):
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
    macro_mol.mol
        The rdkit molecule held in this attribute is replaced by an 
        rdkit molecule with an optimized structure.
    
    macro_mol.file's content
        The content of the structure file located at `macro_mol.file`, 
        is changed so that it holds the structure of the optimized 
        molecule.
    
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
        print('Skipping {0}.'.format(macro_mol.file))       
        return macro_mol
    
    print('\nOptimizing {0}.'.format(macro_mol.file))    
    try:
        # MacroModel requires a ``.mae`` file as input. This creates a 
        # ``.mae`` file holding the molecule.    
        _create_mae(macro_mol, macromodel_path)        
        # generate the ``.com`` file for the MacroModel run.
        _generate_com(macro_mol, force_field, no_fix)        
        # Run the optimization.
        _run_bmin(macro_mol, macromodel_path)
        # Get the ``.maegz`` file output from the optimization and 
        # convert it to a ``.mae`` file.
        _convert_maegz_to_mae(macro_mol, macromodel_path)
        macro_mol.update_from_mae() 

        if md:
            macromodel_md_opt(macro_mol, macromodel_path)
            
        macro_mol.optimized = True       
        return macro_mol 

    except _ConversionError as ex:
        MolError(ex, macro_mol, '`structconvert` failed.')
        return macro_mol        

    except _PathError as ex:
        MolError(ex, macro_mol, 'Wrong MacroModel path.')
        return macro_mol

    except _OptimizationError as ex:
        MolError(ex, macro_mol, 'Optimization by `bmin` failed.')
        return macro_mol

    except _ForceFieldError as ex:        
        # If OPLS_2005 has been tried already - record an exception.
        if force_field == 14:
            MolError(Exception(), macro_mol, 
                          'Both force fields failed.')
            return macro_mol
            
        # If OPLSE_2005 has not been tried - try it.
        print(('Minimization with OPLS3 failed. Trying OPLS_2005. '
               '- {0}').format(macro_mol.file))
        return macromodel_opt(macro_mol, force_field=14,
                              lewis_fixed=lewis_fixed,
                              macromodel_path=macromodel_path,
                              no_fix=no_fix, md=md)

    except _LewisStructureError as ex:
        if not lewis_fixed:
            _run_applyhtreat(macro_mol, macromodel_path)
            return macromodel_opt(macro_mol, force_field=force_field,
                              lewis_fixed=True,
                              macromodel_path=macromodel_path,
                              no_fix=no_fix, md=md) 
        else:
            MolError(ex, macro_mol, ('A viable Lewis'
                                       ' structure was not generated.'))
            return macro_mol

    except Exception as ex:
        MolError(ex, macro_mol, ('Uncategorized '
                      'exception during `macromodel_opt`.'))
        return macro_mol
       
def macromodel_md_opt(macro_mol, macromodel_path, lewis_fixed=False,
                      timeout=True, force_field=16, 
                      temp=300, confs=50, eq_time=10, sim_time=200):  

    print('\nRunning MD on {0}.'.format(macro_mol.file))    
    try:
        # MacroModel requires a ``.mae`` file as input. This creates a 
        # ``.mae`` file holding the molecule.    
        _create_mae(macro_mol, macromodel_path)     
        # Generate the ``.com`` file for the MacroModel MD run.
        _generate_md_com(macro_mol, force_field=force_field, temp=temp, 
                        confs=confs, eq_time=eq_time, sim_time=sim_time)
        # Run the optimization.
        _run_bmin(macro_mol, macromodel_path, timeout)        
        # Extract the lowest energy conformer into its own .mae file.        
        conformer_mae = MAEExtractor(macro_mol).path
        macro_mol.update_from_mae(conformer_mae) 

    except _ConversionError as ex:
        MolError(ex, macro_mol, '`structconvert` failed.')
        return macro_mol

    except _PathError as ex:
        MolError(ex, macro_mol, 'Wrong MacroModel path.')
        return macro_mol

    except _OptimizationError as ex:
        MolError(ex, macro_mol, 'Optimization by `bmin` failed.')
        return macro_mol

    except _ForceFieldError as ex:        
        # If OPLS_2005 has been tried already - record an exception.
        if force_field == 14:
            MolError(Exception(), macro_mol, 
                          'Both force fields failed.')
            return macro_mol
            
        # If OPLSE_2005 has not been tried - try it.
        print(('Minimization with OPLS3 failed. Trying OPLS_2005. '
               '- {0}').format(macro_mol.file))
        return macromodel_md_opt(macro_mol, macromodel_path, 
                                 timeout=timeout, force_field=14,
                                 lewis_fixed=lewis_fixed,
                                 temp=temp, confs=confs, 
                                 eq_time=eq_time, sim_time=sim_time) 

    except _LewisStructureError as ex:
        if not lewis_fixed:
            _run_applyhtreat(macro_mol, macromodel_path)
            return macromodel_md_opt(macro_mol, macromodel_path, 
                                 timeout=timeout, force_field=force_field,
                                 lewis_fixed=True,
                                 temp=temp, confs=confs, 
                                 eq_time=eq_time, sim_time=sim_time) 
        else:
            MolError(ex, macro_mol, ('A viable Lewis'
                                       ' structure was not generated.'))
            return macro_mol
        
    except Exception as ex:
        MolError(ex, macro_mol, ('Uncategorized'
                       ' exception during `macromodel_md_opt`.'))
        return macro_mol
    
def macromodel_cage_opt(macro_mol, force_field=16,
                 macromodel_path=r"C:\Program Files\Schrodinger2016-2",
                 no_fix=False, md=False, lewis_fixed=False):
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
    macro_mol.mol
        The rdkit molecule held in this attribute is replaced by a 
        rdkit molecule with an optimized structure.
    
    macro_mol.file's content
        The content of the structure file located at `macro_mol.file`, 
        is changed so that it holds the structure of the optimized 
        molecule.
    
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
        print('Skipping {0}.'.format(macro_mol.file))       
        return macro_mol
    
    print('\nOptimizing {0}.'.format(macro_mol.file))    
    try:    
        # MacroModel requires a ``.mae`` file as input. This creates a 
        # ``.mae`` file holding the molecule.    
        _create_mae(macro_mol, macromodel_path)        
        # generate the ``.com`` file for the MacroModel run.
        _generate_com(macro_mol, force_field, no_fix)
        # Run the optimization.
        _run_bmin(macro_mol, macromodel_path)
        # Get the ``.maegz`` file output from the optimization and 
        # convert it to a ``.mae`` file.
        _convert_maegz_to_mae(macro_mol, macromodel_path)
        macro_mol.update_from_mae() 

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if macro_mol.topology.windows is not None:
                all_windows = (len(macro_mol.topology.windows) == 
                                       macro_mol.topology.n_windows)

                if md and all_windows:
                    macromodel_md_opt(macro_mol, macromodel_path)

        macro_mol.optimized = True       
        return macro_mol

    except _ConversionError as ex:
        MolError(ex, macro_mol, '`structconvert` failed.')
        return macro_mol        

    except _PathError as ex:
        MolError(ex, macro_mol, 'Wrong MacroModel path.')
        return macro_mol

    except _OptimizationError as ex:
        MolError(ex, macro_mol, 'Optimization by `bmin` failed.')
        return macro_mol

    except _ForceFieldError as ex:        
        # If OPLS_2005 has been tried already - record an exception.
        if force_field==14:
            MolError(Exception(), macro_mol, 
                          'Both force fields failed.')
            return macro_mol
            
        # If OPLSE_2005 has not been tried - try it.
        print(('Minimization with OPLS3 failed. Trying OPLS_2005. '
       '- {0}').format(macro_mol.file))
        return macromodel_cage_opt(macro_mol, force_field=14, 
                              macromodel_path=macromodel_path,
                              no_fix=no_fix, md=md, 
                              lewis_fixed=lewis_fixed)    

    except _LewisStructureError as ex:
        if not lewis_fixed:
            _run_applyhtreat(macro_mol, macromodel_path)
            return macromodel_cage_opt(macro_mol, force_field=force_field,
                              lewis_fixed=True,
                              macromodel_path=macromodel_path,
                              no_fix=no_fix, md=md) 
        else:
            MolError(ex, macro_mol, ('A viable Lewis'
                                       ' structure was not generated.'))
            return macro_mol

    except Exception as ex:
        MolError(ex, macro_mol, ('Uncategorized '
                      'exception during `macromodel_cage_opt`.'))
        return macro_mol

def _run_bmin(macro_mol, macromodel_path, timeout=True):

    print("", time.ctime(time.time()),
    'Running bmin - {0}.'.format(macro_mol.file), sep='\n')
    
    # To run MacroModel a command is issued to to the console via
    # ``subprocess.Popen``. The command is the full path of the ``bmin``
    # program. ``bmin`` is located in the Schrodinger installation
    # folder.
    file_root, ext = os.path.splitext(macro_mol.file)
    log_file = file_root + '.log'
    opt_app = os.path.join(macromodel_path, "bmin")
    # The first member of the list is the command, the following ones
    # are any additional arguments.

    opt_cmd = [opt_app, file_root, "-WAIT", "-LOCAL"]
    opt_proc = psutil.Popen(opt_cmd, stdout=sp.PIPE, 
                            stderr=sp.STDOUT, 
                            universal_newlines=True)
    try:
        if timeout:
            proc_out, _ = opt_proc.communicate(timeout=600)
        else:
            proc_out, _ = opt_proc.communicate()    
    
    except sp.TimeoutExpired:
        print(('\nMinimization took too long and was terminated '
               'by force - {}\n').format(macro_mol.file))
        _kill_bmin(macro_mol, macromodel_path)
        proc_out = ""

    with open(log_file, 'r') as log: 
        log_content = log.read()

    # Check the log for error reports.
    if "termination due to error condition           21-" in log_content:
        raise _OptimizationError(("`bmin` crashed due to"
                                " an error condition. See .log file."))
    if "FATAL do_nosort_typing: NO MATCH found for atom " in log_content:
        raise _ForceFieldError('The log implies the force field failed.')
    if (("FATAL gen_lewis_structure(): could not find best Lewis"
         " structure") in log_content and ("skipping input structure  "
         "due to forcefield interaction errors") in log_content):
        raise _LewisStructureError(
                '`bmin` failed due to poor Lewis structure.')
        

    # If optimization fails because a wrong Schrodinger path was given,
    # raise.
    if 'The system cannot find the path specified' in proc_out:
        raise _PathError(('Wrong Schrodinger path supplied to'
                              ' `macromodel_opt` function.'))

    # If optimization fails because the license is not found, rerun the
    # function.
    if not _license_found(proc_out, macro_mol):
        return _run_bmin(macro_mol, macromodel_path)

    # Make sure the .maegz file created by the optimization is present.
    maegz = file_root +  '-out.maegz'
    _wait_for_file(maegz)
    if not os.path.exists(log_file) or not os.path.exists(maegz):
        raise _OptimizationError(('The .log and/or .maegz '
                     'files were not created by the optimization.'))
        
def _kill_bmin(macro_mol, macromodel_path):
    name, ext = os.path.splitext(macro_mol.file)
    name = re.split(r'\\|/', name)[-1]
    app = os.path.join(macromodel_path, 'jobcontrol')
    cmd = [app, '-stop', name]
    out = sp.run(cmd, stdout=sp.PIPE, 
                 stderr=sp.STDOUT, universal_newlines=True)
 
    # If no license if found, keep re-running the function until it is.
    if not _license_found(out.stdout):
        return _kill_bmin(macro_mol, macromodel_path)  
   
   
   # This loop causes the function to wait until the job has been killed
   # via job control. This means the output files will have been written
   # by the time the function exits. Essentially the loop continues
   # until the job is no longer found by "./jobcontrol -list"
    cmd = [app, '-list']
    output = name
    start = time.time()
    while name in output:
        output = sp.run(cmd, stdout=sp.PIPE, 
                 stderr=sp.STDOUT, universal_newlines=True).stdout
        if time.time() - start > 600:
            break
                 
def _run_applyhtreat(macro_mol, macromodel_path):
    name, ext = os.path.splitext(macro_mol.file)
    mae = name + '.mae'
    mae_out = name + '_htreated.mae'
    _create_mae(macro_mol, macromodel_path)
    
    app = os.path.join(macromodel_path, 'utilities', 'applyhtreat')
    cmd = [app, mae, mae_out]
    out = sp.run(cmd, stdout=sp.PIPE, 
                 stderr=sp.STDOUT, universal_newlines=True)

    # If no license if found, keep re-running the function until it is.
    if not _license_found(out.stdout):
        return _run_applyhtreat(macro_mol, macromodel_path)     
    
    macro_mol.update_from_mae(mae_out)    
    
def _license_found(output, macro_mol=None):
    """
    Checks to see if minimization failed due to a missing license.

    The user can be notified of this in one of two ways. Sometimes the
    output of the submission contains the message informing that the 
    license was not found and in other cases it will be the log file.
    This function checks both of these sources for this message.

    Parameters
    ----------
    output : str
        The outout from submitting the minimization of the structure
        to the ``bmin`` program.

    macro_mol : MacroMolecule (default=None)
        The macromolecule being optimized. If the .log file is not to
        be checked, the default ``None`` should be used.

    Returns
    -------
    bool
        ``True`` if the license was found. ``False`` if the minimization
        did not occur due to a missing license.
    
    """

    if 'Could not check out a license for mmlibs' in output:
        return False
    if macro_mol is None:
        return True

    # To check if the log file mentions a missing license file open the
    # the log file and scan for the apporpriate string.
    
    # Check if the file exists first. If not, this is often means the
    # calculation must be redone so return False anyway.
    log_file_path = macro_mol.file.replace('mol', 'log')
    with open(log_file_path, 'r') as log_file:
        log_file_content = log_file.read()
        
    if 'Could not check out a license for mmlibs' in log_file_content:
        return False

    return True
 
def _generate_com(macro_mol, force_field=16, no_fix=False):
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
    for optimizing the macromolecule using MacroModel.

    Returns
    -------
    None : NoneType    
    
    """

    print('Creating .com file - {}.'.format(macro_mol.file))
    
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
    # structure file but with a ``.com`` extension. Get the path of the
    # ``.mae`` file and the output file in the same way.
    name, ext = os.path.splitext(macro_mol.file)
    com_file = name + '.com'
    mae = name + '.mae'
    output = name + '-out.maegz'
    
    # This function adds all the lines which fix bond distances and 
    # angles into ``main_string``.
    main_string = _fix_params_in_com_file(macro_mol, 
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

def _generate_md_com(macro_mol, force_field=16, temp=300, confs=50, eq_time=10, sim_time=200):

    print('Creating .com file - {0}.'.format(macro_mol.file))

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

    name, ext = os.path.splitext(macro_mol.file)
    com_file = name + '.com'
    mae = name + '.mae'
    output = name + '-out.maegz'

    # Generate the com file containing the info for the run
    with open(com_file, "w") as com:
        # name of the macromodel file
        com.write(str(mae + "\n"))
        # name of the output file
        com.write(str(output + "\n"))
        # details of the macromodel run
        com.write(main_string)

def _create_mae(macro_mol, macromodel_path):
    """
    Creates the ``.mae`` file holding the molecule to be optimized.    

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule which is to be optimized. Its molecular 
        structure file is converted to a ``.mae`` file. The original 
        file is also kept.
        
    macromodel_path : str
        The full path of the installation directory of the Schrodinger 
        suite. By default on a Windows machine it should be something
        like: "C:\Program Files\Schrodinger2016-2".

    Modifies
    --------
    This function creates a new ``.mae`` file from the structure file in
    `macro_mol.file`. This new file is placed in the same
    folder as the original file and has the same name. Only the 
    extensions are different.

    Returns
    -------
    str
        The full path of the newly created ``.mae`` file.     
    
    """
 
    _, ext = os.path.splitext(macro_mol.file)
    
    print('Converting {} to .mae - {}.'.format(ext,
                                              macro_mol.file))
   
    # Create the name of the new ``.mae`` file. It is the same as the
    # original structure file, including the same path. Only the 
    # extensions are different.
    mae_file = macro_mol.file.replace(ext, '.mae')  
    _structconvert(macro_mol.file, mae_file, macromodel_path)
    return mae_file

def _convert_maegz_to_mae(macro_mol, macromodel_path):
    """
    Converts a ``.maegz`` file to a ``.mae`` file.
    
    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule being optimized. The ``.maegz`` file holding 
        its optimized structure is converted to a ``.mae`` file. Both
        versions are kept.
        
    macromodel_path : str
        The full path of the installation directory of the Schrodinger 
        suite. By default on a Windows machine it should be something
        like: "C:\Program Files\Schrodinger2016-2".   
    
    Modifies
    --------    
    This function creates a new ``.mae`` file from a ``.maegz`` file. 
    This new file is placed in the same folder as the ``.maegz`` file.
    
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

    print('Converting .maegz to .mae - {}.'.format(macro_mol.file))

    name, ext = os.path.splitext(macro_mol.file)
    # ``out`` is the full path of the optimized ``.mae`` file.
    maegz = name + '-out.maegz'      
    # Replace extensions to get the names of the various files.
    mae = name + '.mae'
    return _structconvert(maegz, mae, macromodel_path)
    
def _structconvert(iname, oname, macromodel_path):
  
    convrt_app = os.path.join(macromodel_path, 'utilities', 
                                                     'structconvert')
    convrt_cmd = [convrt_app, iname, oname]

    # Execute the file conversion.
    try:
        convrt_return = sp.run(convrt_cmd, stdout=sp.PIPE, 
                           stderr=sp.STDOUT, universal_newlines=True) 

    # If conversion fails because a wrong Schrodinger path was given,
    # raise.
    except FileNotFoundError:
        raise _PathError(('Wrong Schrodinger path supplied to'
                            ' `structconvert` function.'))

    # If no license if found, keep re-running the function until it is.
    if not _license_found(convrt_return.stdout):
        return _structconvert(iname, oname, macromodel_path)    

    # If force field failed, raise.
    if 'number 1' in convrt_return.stdout:
        raise _ForceFieldError(convrt_return.stdout)    

    _wait_for_file(oname)
    if not os.path.exists(oname):
        raise _ConversionError(('Conversion output file {} was not found.'
                ' Console output was {}.').format(oname, 
                                                  convrt_return.stdout))

    return convrt_return

def _fix_params_in_com_file(macro_mol, main_string, no_fix=False):
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
    fix_block = _fix_distance_in_com_file(macro_mol, fix_block)  
    # Add lines that fix the bond angles.                          
    fix_block = _fix_bond_angle_in_com_file(macro_mol, fix_block)
    # Add lines that fix the torsional angles.
    fix_block = _fix_torsional_angle_in_com_file(macro_mol, fix_block)
    
    return main_string.replace(("!!!BLOCK_OF_FIXED_PARAMETERS_"
                                "COMES_HERE!!!\n"), fix_block)

def _fix_distance_in_com_file(macro_mol, fix_block):
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
   
    # Go through all the bonds in the rdkit molecule. If the bond
    # is not between bonder atoms get its distance. Add a fix line using
    # the bond distance and atomic indices to the ``fix_block``. If the
    # bond does invovle two bonder atoms go to the next bond. This is
    # because a bond between 2 bonder atoms was added during assembly 
    # and should therefore not be fixed.  
    for bond in macro_mol.mol.GetBonds():
        atom1 = bond.GetBeginAtom() 
        atom2 = bond.GetEndAtom()
        
        if atom1.HasProp('bonder') and atom2.HasProp('bonder'):
            continue
        
        atom1_id = atom1.GetIdx() 
        atom2_id = atom2.GetIdx() 
        
        bond_len = macro_mol.atom_distance(atom1_id, atom2_id)
        
        # Make sure that the indices are increased by 1 in the ``.mae``
        # file from their rdkit value.
        fix_block += (fix_distance.format(atom1_id+1, atom2_id+1,
                                         bond_len) + "\n")

    return fix_block
    
def _fix_bond_angle_in_com_file(macro_mol, fix_block):
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
    ba_atoms = macro_mol.mol.GetSubstructMatches(ba_mol)
    
    # Get the conformer holding the atomic positions.
    conf = macro_mol.mol.GetConformer()    
    
    # For each bond angle check if a bonder atom is involved in forming
    # it. If no, a line fixing the bond angle is added to ``fix_block``.
    # If any atom of the 3 is a bonder atom the bond angle is not fixed.
    # This means that there will be some bond angles which consist of 2
    # bonds not added during assembly which will not be fixed. However,
    # it is assumed that the effect of this will be minimal.
    for atom1_id, atom2_id, atom3_id in ba_atoms:
        atom1 = macro_mol.mol.GetAtomWithIdx(atom1_id)
        atom2 = macro_mol.mol.GetAtomWithIdx(atom2_id)
        atom3 = macro_mol.mol.GetAtomWithIdx(atom3_id)
        
        if (atom1.HasProp('bonder') or atom2.HasProp('bonder') or
            atom3.HasProp('bonder')):
            continue
        
        ba = ac.GetAngleDeg(conf, atom1_id, atom2_id, atom3_id)
        
        fix_block += (fix_ba.format(atom1_id+1, atom2_id+1, 
                                    atom3_id+1, ba) + "\n")
    
    return fix_block
    
def _fix_torsional_angle_in_com_file(macro_mol, fix_block):
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
    ta_atoms = macro_mol.mol.GetSubstructMatches(ta_mol)
    # Get the conformer holding the atomic positions.
    conf = macro_mol.mol.GetConformer()
    
    # For each torsional angle check if a bonder atom is involved in 
    # forming it. If no, a line fixing the torsional angle is added to 
    # ``fix_block``. If any atom of the 4 is a bonder atom the bond 
    # angle is not fixed. This means that there will be some bond angles 
    # which consist of 3 bonds not added during assembly which will not 
    # be fixed. However, it is assumed that the effect of this will be 
    # minimal.    
    for atom1_id, atom2_id, atom3_id, atom4_id in ta_atoms:
        atom1 = macro_mol.mol.GetAtomWithIdx(atom1_id)
        atom2 = macro_mol.mol.GetAtomWithIdx(atom2_id)
        atom3 = macro_mol.mol.GetAtomWithIdx(atom3_id)
        atom4 = macro_mol.mol.GetAtomWithIdx(atom4_id)
        
        if (atom1.HasProp('bonder') or atom2.HasProp('bonder') or
            atom3.HasProp('bonder') or atom4.HasProp('bonder')):
            continue
        
        ta = ac.GetDihedralDeg(conf, atom1_id, atom2_id, 
                                     atom3_id, atom4_id)
        
        fix_block += (fix_ta.format(atom1_id+1, atom2_id+1, 
                                atom3_id+1, atom4_id+1, ta) + "\n")

    return fix_block

def _wait_for_file(file_name, timeout=10):
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