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
from uuid import uuid4

from ...convenience_tools import MAEExtractor

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

def macromodel_opt(macro_mol, macromodel_path, settings={}, md={}):
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

    settings : dict (default = {})
        A dictionary which maps the names of optimization parameters to
        their values. Valid values are:

            'restricted' : bool (default = True)
                If False then all bonds are optimized, not just the
                ones created during macromolecular assembly.

            'timeout' : float (default = 0)
                The amount in seconds the optimization is allowed to
                run before being terminated. 0 means there is no
                timeout.

            'force_field' : int (default = 16)
                The number of the force field to be used.

            'max_iter' : int (default = 2500)
                The maximum number of iterations done during the
                optimization.

            'gradient' : float (default = 0.05)
                The gradient at which optimization is stopped.

            'md' : bool (default=False)
                Toggles whether a MD conformer search should be
                performed.

        Only values which need to be changed from the default need to
        be specified. For exmaple:

            macromodel_opt(mol, 'path', {'max_iter' : 10})


    md : dict (default = {})
        A dictionary holding settings for the MD conformer search.
        This parameter is used in the same way as `settings` except
        the values effect the MD only. See docstring of
        macromodel_md_opt() for valid values.

    Modifies
    --------
    macro_mol.mol
        The rdkit molecule held in this attribute is replaced by an
        rdkit molecule with an optimized structure.

    Returns
    -------
    None : NoneType

    """

    vals = {
             'restricted' : True,
             'timeout' : 0,
             'force_field' : 16,
             'max_iter' : 2500,
             'gradient' : 0.05,
             'md' : False
            }
    vals.update(settings)

    # Default initiailize some parameters. These are for the internal
    # use of the function (to prevent infinite recursion) not for the
    # user. Also see comments below.
    if 'lewis_fixed' not in vals:
        vals['lewis_fixed'] = False

    try:
        macro_mol._file = getattr(macro_mol, '_file',
                                 '{}.mol'.format(uuid4().int))
        # First write a .mol file of the molecule.
        macro_mol.write(macro_mol._file)
        # MacroModel requires a ``.mae`` file as input. This creates a
        # ``.mae`` file holding the molecule.
        _create_mae(macro_mol, macromodel_path)
        # generate the ``.com`` file for the MacroModel run.
        _generate_com(macro_mol, vals)
        # Run the optimization.
        _run_bmin(macro_mol, macromodel_path)
        # Get the ``.maegz`` file output from the optimization and
        # convert it to a ``.mae`` file.
        _convert_maegz_to_mae(macro_mol, macromodel_path)
        macro_mol.update_from_mae(
                            macro_mol._file.replace('.mol', '.mae'))

        if vals['md']:
            macromodel_md_opt(macro_mol, macromodel_path, md)

    except _ForceFieldError as ex:
        # If OPLS_2005 has been tried already - record an exception.
        if vals['force_field'] == 14:
            raise ex

        # If OPLSE_2005 has not been tried - try it.
        print(('Minimization with OPLS3 failed. Trying OPLS_2005. '
               '- {}').format(macro_mol.name))

        vals['force_field'] = 14
        return macromodel_opt(macro_mol,
                              macromodel_path, vals, md)

    # The 'lewis_fixed' parameter should not be used by the user.
    # Sometimes the Lewis structure of `macro_mol` is wrong. If this is
    # the case the function tries to fix it and then runs itself again.
    # The `lewis_fixed` parameter indicates if the fix has already been
    # tried to prevent infinite recursion.
    except _LewisStructureError as ex:
        print('Attempting to fix Lewis structure.')
        if not vals['lewis_fixed']:
            _run_applyhtreat(macro_mol, macromodel_path)
            vals['lewis_fixed'] = True
            return macromodel_opt(macro_mol,
                                  macromodel_path, vals, md)
        else:
            raise ex

def macromodel_md_opt(macro_mol, macromodel_path, settings={}):
    """
    Runs a MD conformer search on `macro_mol`.

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule who's structure must be optimized.

    macromodel_path : str
        The full path of the ``Schrodinger`` suite within the user's
        machine. For example, in a default Microsoft installation the
        folder will probably be something like
        ``C:\Program Files\Schrodinger2016-2``.

    settings : dict (default = {})
        Each key is a string describing an MD parameter. The value is
        the corresponding value. Only values which need to be changed
        from the default need to be specified. The allowed parameters
        are:

            'timeout' : float (default = 0)
                The amount in seconds the MD is allowed to run before
                being terminated. 0 means there is no timeout.

            'force_field' : int (default = 16)
                The number of the force field to be used.

            'temp' : float (default = 300)
                The temperature in Kelvin at which the MD is run.

            'confs' : int (default = 50)
                The number of conformers sampled and optimized from the
                MD.

            'time_step' : float (default = 1.5)
                The time step for the MD.

            'eq_time' : float (default = 10)
                The equilibriation time before the MD is run.

            'sim_time' : float (default = 200)
                The simulation time of the MD.

            'max_iter' : int (default = 2500)
                The maximum number of iterations done during the
                optimization.

            'gradient' : float (default = 0.05)
                The gradient at which optimization is stopped.

    Modifies
    --------
    macro_mol.mol
        The rdkit molecule held in this attribute is replaced by an
        rdkit molecule with an optimized structure.

    Returns
    -------
    None : NoneType

    """

    vals = {
               'timeout' : 0,
               'force_field' : 16,
               'temp' : 300,
               'confs' : 50,
               'time_step' : 1.5,
               'eq_time' : 10,
               'sim_time' : 200,
               'max_iter' : 2500,
               'gradient' : 0.05
              }

    vals.update(settings)

    # Default initiailize some parameters. These are for the internal
    # use of the function (to prevent infinite recursion) not for the
    # user. Also see comments below.
    if 'lewis_fixed' not in vals:
        vals['lewis_fixed'] = False

    print('\nRunning MD on {}.'.format(macro_mol.name))
    try:
        macro_mol._file = getattr(macro_mol, '_file',
                                 '{}.mol'.format(uuid4().int))
        # First write a .mol file of the molecule.
        macro_mol.write(macro_mol._file)
        # MacroModel requires a ``.mae`` file as input. This creates a
        # ``.mae`` file holding the molecule.
        _create_mae(macro_mol, macromodel_path)
        # Generate the ``.com`` file for the MacroModel MD run.
        _generate_md_com(macro_mol, vals)
        # Run the optimization.
        _run_bmin(macro_mol, macromodel_path, vals['timeout'])
        # Extract the lowest energy conformer into its own .mae file.
        conformer_mae = MAEExtractor(macro_mol._file).path
        macro_mol.update_from_mae(conformer_mae)

    except _ForceFieldError as ex:
        # If OPLS_2005 has been tried already - record an exception.
        if vals['force_field'] == 14:
            raise ex
        # If OPLSE_2005 has not been tried - try it.
        print(('Minimization with OPLS3 failed. Trying OPLS_2005. '
               '- {}').format(macro_mol.name))

        vals['force_field'] = 14
        return macromodel_md_opt(macro_mol, macromodel_path, vals)

    # The 'lewis_fixed' parameter should not be used by the user.
    # Sometimes the Lewis structure of `macro_mol` is wrong. If this is
    # the case the function tries to fix it and then runs itself again.
    # The `lewis_fixed` parameter indicates if the fix has already been
    # tried to prevent infinite recursion.
    except _LewisStructureError as ex:
        print('Attempting to fix Lewis structure.')
        if not vals['lewis_fixed']:
            vals['lewis_fixed'] = True
            _run_applyhtreat(macro_mol, macromodel_path)
            return macromodel_md_opt(macro_mol,
                                     macromodel_path, vals)
        else:
            raise ex

def macromodel_cage_opt(macro_mol,
                        macromodel_path, settings={}, md={}):
    """
    Optimizes the molecule using MacroModel.

    This function runs a restricted optimization. The structures of the
    building blocks are frozen and only the new bonds formed between
    building blocks during assembly are optimized.

    This function differes from `macromodel_opt` in that it checks the
    number of windows the `macro_mol` has before running the MD. The MD
    is only run all windows are found  (i.e. the cage is not
    collapsed).

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule who's structure must be optimized.

    macromodel_path : str
        The full path of the ``Schrodinger`` suite within the user's
        machine. For example, in a default Microsoft installation the
        folder will probably be something like
        ``C:\Program Files\Schrodinger2016-2``.

    settings : dict (default = {})
        A dictionary which maps the names of optimization parameters to
        their values. Valid values are:

            'restricted' : bool (default = True)
                If False then all bonds are optimized, not just the
                ones created during macromolecular assembly.

            'timeout' : float (default = 0)
                The amount in seconds the optimization is allowed to
                run before being terminated. 0 means there is no
                timeout.

            'force_field' : int (default = 16)
                The number of the force field to be used.

            'max_iter' : int (default = 2500)
                The maximum number of iterations done during the
                optimization.

            'gradient' : float (default = 0.05)
                The gradient at which optimization is stopped.

            'md' : bool (default=False)
                Toggles whether a MD conformer search should be
                performed.

        Only values which need to be changed from the default need to
        be specified. For exmaple:

            macromodel_opt(mol, 'path', {'max_iter' : 10})


    md : dict (default = {})
        A dictionary holding settings for the MD conformer search.
        This parameter is used in the same way as `settings` except
        the values effect the MD only. See docstring of
        macromodel_md_opt() for valid values.

    Modifies
    --------
    macro_mol.mol
        The rdkit molecule held in this attribute is replaced by an
        rdkit molecule with an optimized structure.

    Returns
    -------
    None : NoneType

    """

    vals = {
             'restricted' : True,
            'timeout' : 0,
             'force_field' : 16,
             'max_iter' : 2500,
             'gradient' : 0.05,
             'md' : False
            }
    vals.update(settings)

    # Default initialize some parameters. These are for the internal
    # use of the function (to prevent infinite recursion) not for the
    # user. Also see comments below.
    if 'lewis_fixed' not in vals:
        vals['lewis_fixed'] = False

    try:
        macro_mol._file = getattr(macro_mol, '_file',
                                 '{}.mol'.format(uuid4().int))
        # First write a .mol file of the molecule.
        macro_mol.write(macro_mol._file)
        # MacroModel requires a ``.mae`` file as input. This creates a
        # ``.mae`` file holding the molecule.
        _create_mae(macro_mol, macromodel_path)
        # generate the ``.com`` file for the MacroModel run.
        _generate_com(macro_mol, vals)
        # Run the optimization.
        _run_bmin(macro_mol, macromodel_path)
        # Get the ``.maegz`` file output from the optimization and
        # convert it to a ``.mae`` file.
        _convert_maegz_to_mae(macro_mol, macromodel_path)
        macro_mol.update_from_mae(
                             macro_mol._file.replace('.mol', '.mae'))

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if macro_mol.windows is not None:
                all_windows = (len(macro_mol.windows) ==
                                       macro_mol.topology.n_windows)

                if vals['md'] and all_windows:
                    macromodel_md_opt(macro_mol, macromodel_path, md)

    except _ForceFieldError as ex:
        # If OPLS_2005 has been tried already - record an exception.
        if vals['force_field'] == 14:
            raise ex

        # If OPLSE_2005 has not been tried - try it.
        print(('Minimization with OPLS3 failed. Trying OPLS_2005. '
        '- {}').format(macro_mol.name))
        vals['force_field'] = 14
        return macromodel_cage_opt(macro_mol, macromodel_path, vals)

    # The 'lewis_fixed' parameter should not be used by the user.
    # Sometimes the Lewis structure of `macro_mol` is wrong. If this is
    # the case the function tries to fix it and then runs itself again.
    # The `lewis_fixed` parameter indicates if the fix has already been
    # tried to prevent infinite recursion.
    except _LewisStructureError as ex:
        print('Attempting to fix Lewis structure.')
        if not vals['lewis_fixed']:
            _run_applyhtreat(macro_mol, macromodel_path)
            vals['lewis_fixed'] = True
            return macromodel_cage_opt(macro_mol,
                                       macromodel_path, vals)
        else:
            raise ex

def _run_bmin(macro_mol, macromodel_path, timeout=0):

    print("", time.ctime(time.time()),
    'Running bmin - {}.'.format(macro_mol.name), sep='\n')

    # To run MacroModel a command is issued to the console via
    # ``subprocess.Popen``. The command is the full path of the
    # ``bmin`` program. ``bmin`` is located in the Schrodinger
    # installation folder.
    file_root, ext = os.path.splitext(macro_mol._file)
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
            proc_out, _ = opt_proc.communicate(timeout=timeout)
        else:
            proc_out, _ = opt_proc.communicate()


    except sp.TimeoutExpired:
        print(('\nMinimization took too long and was terminated '
               'by force - {}\n').format(macro_mol.name))
        _kill_bmin(macro_mol, macromodel_path)
        proc_out = ""

    with open(log_file, 'r') as log:
        log_content = log.read()

    # Check the log for error reports.
    if ("termination due to error condition           21-" in
                                                         log_content):
        raise _OptimizationError(("`bmin` crashed due to"
                            " an error condition. See .log file."))

    if ("FATAL do_nosort_typing: NO MATCH found for atom " in
                                                         log_content):
        raise _ForceFieldError(
                        'The log implies the force field failed.')

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
    name, ext = os.path.splitext(macro_mol._file)
    name = re.split(r'\\|/', name)[-1]
    app = os.path.join(macromodel_path, 'jobcontrol')
    cmd = [app, '-stop', name]
    out = sp.run(cmd, stdout=sp.PIPE,
                 stderr=sp.STDOUT, universal_newlines=True)

    # If no license if found, keep re-running the function until it is.
    if not _license_found(out.stdout):
        return _kill_bmin(macro_mol, macromodel_path)


   # This loop causes the function to wait until the job has been
   # killed via job control. This means the output files will have been
   # written by the time the function exits. Essentially the loop
   # continues until the job is no longer found by "./jobcontrol -list"
    cmd = [app, '-list']
    output = name
    start = time.time()
    while name in output:
        output = sp.run(cmd, stdout=sp.PIPE,
                 stderr=sp.STDOUT, universal_newlines=True).stdout
        if time.time() - start > 600:
            break

def _run_applyhtreat(macro_mol, macromodel_path):
    name, ext = os.path.splitext(macro_mol._file)
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
        ``True`` if the license was found. ``False`` if the
        minimization did not occur due to a missing license.

    """

    if 'Could not check out a license for mmlibs' in output:
        return False
    if macro_mol is None:
        return True

    # To check if the log file mentions a missing license file open the
    # the log file and scan for the apporpriate string.

    # Check if the file exists first. If not, this is often means the
    # calculation must be redone so return False anyway.
    log_file_path = macro_mol._file.replace('mol', 'log')
    with open(log_file_path, 'r') as log_file:
        log_file_content = log_file.read()

    if 'Could not check out a license for mmlibs' in log_file_content:
        return False

    return True

def _com_line(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9):
    return (" {:<5}{:>7}{:>7}{:>7}{:>7}{:>11.4f}{:>11.4f}"
            "{:>11.4f}{:>11.4f}").format(arg1, arg2, arg3, arg4,
                                         arg5, arg6, arg7, arg8, arg9)

def _generate_com(macro_mol, settings):
    """
    Create a ``.com`` file for a MacroModel optimization.

    The created ``.com`` file fixes all bond parameters which were not
    added during assembly. This means all bond distances, bond angles
    and torsional angles are fixed, except for cases where it involves
    a bond added during assembly of the macromolecule.

    This fixing is implemented by creating a ``.com`` file with various
    ``FX`` commands written within its body.

    This function is called by ``macromodel_opt``. It is private
    because it should probably not be used outside of this context.

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule which is to be optimized.

    settings : dict
        A dictionary of settings for the optimization. See
        macromodel_opt() documentation.

    Modifies
    --------
    This function creates a new ``.com`` file holding the instructions
    for optimizing the macromolecule using MacroModel.

    Returns
    -------
    None : NoneType

    """

    print('Creating .com file - {}.'.format(macro_mol.name))

    # This is the body of the ``.com`` file. The line that begins and
    # ends with exclamation lines is replaced with the various commands
    # that fix bond distances and angles.
    main_string = "\n".join([
        _com_line('MMOD', 0, 1, 0, 0, 0, 0, 0, 0),

        _com_line('FFLD',
                  settings['force_field'], 1, 0, 0, 1, 0, 0, 0),

        _com_line('BGIN', 0, 0, 0, 0, 0, 0, 0, 0),
        _com_line('READ', 0, 0, 0, 0, 0, 0, 0, 0),
        "!!!BLOCK_OF_FIXED_PARAMETERS_COMES_HERE!!!",
        _com_line('CONV', 2, 0, 0, 0, settings['gradient'], 0, 0, 0),
        _com_line('MINI', 1, 0, settings['max_iter'], 0, 0, 0, 0, 0),
        _com_line('END', 0, 1, 0, 0, 0, 0, 0, 0)])

    # Create a path for the ``.com`` file. It is the same as that of
    # the structure file but with a ``.com`` extension. Get the path of
    # the ``.mae`` file and the output file in the same way.
    name, ext = os.path.splitext(macro_mol._file)
    com_file = name + '.com'
    mae = name + '.mae'
    output = name + '-out.maegz'

    # This function adds all the lines which fix bond distances and
    # angles into ``main_string``.
    main_string = _fix_params_in_com_file(macro_mol, main_string,
                                          settings['restricted'])

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

def _generate_md_com(macro_mol, settings):
    print('Creating .com file - {}.'.format(macro_mol.name))

    main_string = "\n".join([
        _com_line('MMOD', 0, 1, 0, 0, 0, 0, 0, 0),

        _com_line('FFLD',
                  settings['force_field'], 1, 0, 0, 1, 0, 0, 0),

        _com_line('READ', 0, 0, 0, 0, 0, 0, 0, 0),
        _com_line('MDIT', 0, 0, 0, 0, settings['temp'], 0, 0, 0),

        _com_line('MDYN', 0, 0, 0, 0, settings['time_step'],
                    settings['eq_time'], settings['temp'], 0),

        _com_line('MDSA', settings['confs'], 0, 0, 0, 0, 0, 1, 0),

        _com_line('MDYN', 1, 0, 0, 0, settings['time_step'],
                   settings['sim_time'], settings['temp'], 0),

        _com_line('WRIT', 0, 0, 0, 0, 0, 0, 0, 0),
        _com_line('RWND', 0, 1, 0, 0, 0, 0, 0, 0),
        _com_line('BGIN', 0, 0, 0, 0, 0, 0, 0, 0),
        _com_line('READ', -2, 0, 0, 0, 0, 0, 0, 0),
        _com_line('CONV', 2, 0, 0, 0, settings['gradient'], 0, 0, 0),
        _com_line('MINI', 1, 0, settings['max_iter'], 0, 0, 0, 0, 0),
        _com_line('END', 0, 1, 0, 0, 0, 0, 0, 0)])

    name, ext = os.path.splitext(macro_mol._file)
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
    This function creates a new ``.mae`` file from the structure file
    in `macro_mol._file`. This new file is placed in the same
    folder as the original file and has the same name. Only the
    extensions are different.

    Returns
    -------
    str
        The full path of the newly created ``.mae`` file.

    """

    _, ext = os.path.splitext(macro_mol._file)

    print('Converting {} to .mae - {}.\n'.format(ext, macro_mol.name))

    # Create the name of the new ``.mae`` file. It is the same as the
    # original structure file, including the same path. Only the
    # extensions are different.
    mae_file = macro_mol._file.replace(ext, '.mae')
    _structconvert(macro_mol._file, mae_file, macromodel_path)
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

    print('Converting .maegz to .mae - {}.'.format(macro_mol.name))

    name, ext = os.path.splitext(macro_mol._file)
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
        raise _ConversionError(
        ('Conversion output file {} was not found.'
        ' Console output was {}.').format(oname, convrt_return.stdout))

    return convrt_return

def _fix_params_in_com_file(macro_mol, main_string, restricted):
    """
    Adds lines to the ``.com`` body fixing bond distances and angles.

    For each bond distance, bond angle and torisional angle that does
    not involve a bond created during assembly a ``FX`` command is
    added to the string holding holding the body of the ``.com`` file.

    These lines replace the filler line in the main string.

    This function is called by ``macromodel_opt``. It is private
    because it should probably not be used outside of this context.

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule which is to be optimized.

    main_string : str
        The body of the ``.com`` file which is to have fix commands
        added.

    restricted : bool (default = True)
        When ``False`` the block containing instructions to fix
        molecular parameters is not added to the .com file.

    Returns
    -------
    str
        A string holding the body of the ``.com`` file with
        instructions to fix the various bond distances and angles as
        described in the docstring.

    """

    # Make a string to hold all of the ``FX`` lines.
    fix_block = ""

    # If `restricted` is ``False`` do not add a fix block.
    if not restricted:
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
        The string holding all the lines containing fix commands for
        the ``.com`` file.

    Returns
    -------
    str
        A string holding lines containg fix commands for the ``.com``
        file. The string has the lines fixing bond distances added to
        it by this function.

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
    # is not between bonder atoms get its distance. Add a fix line
    # using the bond distance and atomic indices to the ``fix_block``.
    # If the bond does invovle two bonder atoms go to the next bond.
    # This is because a bond between 2 bonder atoms was added during
    # assembly and should therefore not be fixed.
    for bond in macro_mol.mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()

        if (atom1.GetIdx() in macro_mol.bonder_ids and
            atom2.GetIdx() in macro_mol.bonder_ids):
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

    if the bond between D and E (``=``) represents the bond added
    during assembly, the bond angle A-B-C will be fixed but B-C-D will
    not be. This is an artifact of the implementation but is not
    expected to play a significant role as the vast majority of bond
    angles which should be fixed, will be. The bond angle C-D=E will
    also not be fixed as that is the purpose of this function.

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule to be optimized.

    fix_block : str
        The string holding all the lines containing fix commands for
        the ``.com`` file.

    Returns
    -------
    str
        A string holding lines containg fix commands for the ``.com``
        file. The string has the lines fixing bond angles added to it
        by this function.

    """
    # This line holds the format for a line fixing the bond angles
    # between 3 atoms. The first 3 ``{...}`` are replaced with the ids
    # of atoms. The last ``{...}`` is replaced with the bond angle.
    # Note that in the ``.mae`` files the indices of atoms start at 1
    # while in rdkit they start at 0. As far as I can tell this
    # corresponds to a shift of 1 for each atom index, with the
    # ordering being the same.
    fix_ba = (" FXBA {0:>7}{1:>7}{2:>7}      0   100.0000 "
              "{3:>10.4f}     0.0000     0.0000")

    # Create a substructure consisting of 3 dummy atoms bonded with 3
    # dummy bonds. This substructure will match with any 3 atoms which
    # are bonded together with any combination of bonds. These 3 atoms
    # will therefore have a bond angle.
    ba_mol = chem.MolFromSmarts('[*]~[*]~[*]')

    # Get the indices of all atoms which have a bond angle.
    # ``ba_atoms`` is a tuple of tuples of the form ((1,2,3), (4,5,6),
    # (7,8,9), ...). Each inner tuple holds the indicies of the atoms
    # which form a bond angle.
    ba_atoms = macro_mol.mol.GetSubstructMatches(ba_mol)

    # Get the conformer holding the atomic positions.
    conf = macro_mol.mol.GetConformer()

    # For each bond angle check if a bonder atom is involved in forming
    # it. If no, a line fixing the bond angle is added to
    # ``fix_block``. If any atom of the 3 is a bonder atom the bond
    # angle is not fixed. This means that there will be some bond
    # angles which consist of 2 bonds not added during assembly which
    # will not be fixed. However, it is assumed that the effect of this
    # will be minimal.
    for atom1_id, atom2_id, atom3_id in ba_atoms:
        if (atom1_id in macro_mol.bonder_ids or
            atom2_id in macro_mol.bonder_ids or
            atom3_id in macro_mol.bonder_ids):
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

    if the bond between E and F (``=``) represents the bond added
    during assembly, the torsional angle A-B-C-D will be fixed but
    B-C-D-E will not be. This is an artifact of the implementation but
    is not expected to play a significant role as the vast majority of
    bond angles which should be fixed, will be. The bond angle C-D-E=F
    will also not be fixed as that is the purpose of this function.

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule to be optimized.

    fix_block : str
        The string holding all the lines containing fix commands for
        the ``.com`` file.

    Returns
    -------
    str
        A string holding lines containg fix commands for the ``.com``
        file. The string has the lines fixing bond angles added to it
        by this function.

    """

    # This line holds the format for a line fixing the torsional angles
    # between 4 atoms. The first 4 ``{...}`` are replaced with the ids
    # of atoms. The last ``{...}`` is replaced with the torsional
    # angle. Note that in the ``.mae`` files the indices of atoms start
    # at 1 while in rdkit they start at 0. As far as I can tell this
    # corresponds to a shift of 1 for each atom index, with the
    # ordering being the same.
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
    # angle is not fixed. This means that there will be some bond
    # angles which consist of 3 bonds not added during assembly which
    # will not be fixed. However, it is assumed that the effect of this
    # will be minimal.
    for atom1_id, atom2_id, atom3_id, atom4_id in ta_atoms:
        if (atom1_id in macro_mol.bonder_ids or
            atom2_id in macro_mol.bonder_ids or
            atom3_id in macro_mol.bonder_ids or
            atom4_id in macro_mol.bonder_ids):
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
            print('Waiting for {}.'.format(file_name))
            tick += 1

        if os.path.exists(file_name) or time_taken > timeout:
            break
