"""
Defines optimization functions which use MacroModel.

"""

import os
import subprocess as sp
import time
import rdkit.Chem.AllChem as rdkit
import warnings
import psutil
import re
from uuid import uuid4
import logging

from ..convenience_tools import MAEExtractor


logger = logging.getLogger(__name__)


class _ConversionError(Exception):
    def __init__(self, message):
        self.message = message


class _PathError(Exception):
    def __init__(self, message):
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


def macromodel_opt(macro_mol, macromodel_path,
                   settings={}, md={}, logger=logger):
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

    logger : FakeLogger or logging.Logger, optional
        Used for logging.

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
             'restricted': True,
             'timeout': 0,
             'force_field': 16,
             'max_iter': 2500,
             'gradient': 0.05,
             'md': False
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
        _create_mae(macro_mol, macromodel_path, logger)
        # generate the ``.com`` file for the MacroModel run.
        _generate_com(macro_mol, vals, logger)
        # Run the optimization.
        _run_bmin(macro_mol, macromodel_path, logger, vals['timeout'])
        # Get the ``.maegz`` file output from the optimization and
        # convert it to a ``.mae`` file.
        _convert_maegz_to_mae(macro_mol, macromodel_path, logger)
        macro_mol.update_from_mae(macro_mol._file.replace('.mol',
                                                          '.mae'))

        if vals['md']:
            _macromodel_md_opt(macro_mol, macromodel_path, logger, md)

    except _ForceFieldError as ex:
        # If OPLS_2005 has been tried already - record an exception.
        if vals['force_field'] == 14:
            raise ex

        # If OPLSE_2005 has not been tried - try it.
        logger.warning(('Minimization with OPLS3 failed on "{}". '
                        'Trying OPLS_2005.').format(macro_mol.name))

        vals['force_field'] = 14
        return macromodel_opt(macro_mol, macromodel_path,
                              vals, md, logger)

    # The 'lewis_fixed' parameter should not be used by the user.
    # Sometimes the Lewis structure of `macro_mol` is wrong. If this is
    # the case the function tries to fix it and then runs itself again.
    # The `lewis_fixed` parameter indicates if the fix has already been
    # tried to prevent infinite recursion.
    except _LewisStructureError as ex:
        logger.warning(('Attempting to fix Lewis '
                        'structure of "{}".'.format(macro_mol.name)))
        if not vals['lewis_fixed']:
            _run_applyhtreat(macro_mol, macromodel_path, logger)
            vals['lewis_fixed'] = True
            return macromodel_opt(macro_mol, macromodel_path,
                                  vals, md, logger)
        else:
            raise ex


def macromodel_cage_opt(macro_mol, macromodel_path,
                        settings={}, md={}, logger=logger):
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

    logger : FakeLogger or logging.Logger, optional
        Used for logging.

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
             'restricted': True,
             'timeout': 0,
             'force_field': 16,
             'max_iter': 2500,
             'gradient': 0.05,
             'md': False
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
        _create_mae(macro_mol, macromodel_path, logger)
        # generate the ``.com`` file for the MacroModel run.
        _generate_com(macro_mol, vals, logger)
        # Run the optimization.
        _run_bmin(macro_mol, macromodel_path, logger, vals['timeout'])
        # Get the ``.maegz`` file output from the optimization and
        # convert it to a ``.mae`` file.
        _convert_maegz_to_mae(macro_mol, macromodel_path, logger)
        macro_mol.update_from_mae(macro_mol._file.replace('.mol',
                                                          '.mae'))

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if macro_mol.windows is not None:
                all_windows = (len(macro_mol.windows) ==
                               macro_mol.topology.n_windows)

                if vals['md'] and all_windows:
                    _macromodel_md_opt(macro_mol, macromodel_path,
                                       logger, md)

    except _ForceFieldError as ex:
        # If OPLS_2005 has been tried already - record an exception.
        if vals['force_field'] == 14:
            raise ex

        # If OPLSE_2005 has not been tried - try it.
        logger.warning(('Minimization with OPLS3 failed on "{}". '
                        'Trying OPLS_2005.').format(macro_mol.name))

        vals['force_field'] = 14
        return macromodel_cage_opt(macro_mol, macromodel_path,
                                   vals, logger)

    # The 'lewis_fixed' parameter should not be used by the user.
    # Sometimes the Lewis structure of `macro_mol` is wrong. If this is
    # the case the function tries to fix it and then runs itself again.
    # The `lewis_fixed` parameter indicates if the fix has already been
    # tried to prevent infinite recursion.
    except _LewisStructureError as ex:
        logger.warning(('Attempting to fix Lewis '
                        'structure of "{}".'.format(macro_mol.name)))
        if not vals['lewis_fixed']:
            _run_applyhtreat(macro_mol, macromodel_path, logger)
            vals['lewis_fixed'] = True
            return macromodel_cage_opt(macro_mol,
                                       macromodel_path, vals, logger)
        else:
            raise ex


def _macromodel_md_opt(macro_mol, macromodel_path,
                       logger, settings={}):
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

    logger : FakeLogger or logging.Logger
        Used for logging.

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

            'time_step' : float (default = 1.0)
                The time step in fs for the MD.

            'eq_time' : float (default = 10)
                The equilibriation time in ps before the MD is run.

            'sim_time' : float (default = 200)
                The simulation time in ps of the MD.

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
               'timeout': 0,
               'force_field': 16,
               'temp': 300,
               'confs': 50,
               'time_step': 1.0,
               'eq_time': 10,
               'sim_time': 200,
               'max_iter': 2500,
               'gradient': 0.05
              }

    vals.update(settings)

    # Default initiailize some parameters. These are for the internal
    # use of the function (to prevent infinite recursion) not for the
    # user. Also see comments below.
    if 'lewis_fixed' not in vals:
        vals['lewis_fixed'] = False

    logger.info('Running MD on "{}".'.format(macro_mol.name))
    try:
        macro_mol._file = getattr(macro_mol, '_file',
                                  '{}.mol'.format(uuid4().int))
        # First write a .mol file of the molecule.
        macro_mol.write(macro_mol._file)
        # MacroModel requires a ``.mae`` file as input. This creates a
        # ``.mae`` file holding the molecule.
        _create_mae(macro_mol, macromodel_path, logger)
        # Generate the ``.com`` file for the MacroModel MD run.
        _generate_md_com(macro_mol, vals, logger)
        # Run the optimization.
        _run_bmin(macro_mol, macromodel_path, logger, vals['timeout'])
        # Extract the lowest energy conformer into its own .mae file.
        conformer_mae = MAEExtractor(macro_mol._file).path
        macro_mol.update_from_mae(conformer_mae)

    except _ForceFieldError as ex:
        # If OPLS_2005 has been tried already - record an exception.
        if vals['force_field'] == 14:
            raise ex
        # If OPLSE_2005 has not been tried - try it.
        logger.warning(('Minimization with OPLS3 failed on "{}". '
                        'Trying OPLS_2005.').format(macro_mol.name))

        vals['force_field'] = 14
        return _macromodel_md_opt(macro_mol, macromodel_path, logger,
                                  vals)

    # The 'lewis_fixed' parameter should not be used by the user.
    # Sometimes the Lewis structure of `macro_mol` is wrong. If this is
    # the case the function tries to fix it and then runs itself again.
    # The `lewis_fixed` parameter indicates if the fix has already been
    # tried to prevent infinite recursion.
    except _LewisStructureError as ex:
        logger.warning(('Attempting to fix Lewis '
                        'structure of "{}".'.format(macro_mol.name)))
        if not vals['lewis_fixed']:
            vals['lewis_fixed'] = True
            _run_applyhtreat(macro_mol, macromodel_path, logger)
            return _macromodel_md_opt(macro_mol,
                                      macromodel_path, logger, vals)
        else:
            raise ex


def _run_bmin(macro_mol, macromodel_path, logger, timeout=0):

    logger.info('Running bmin on "{}".'.format(macro_mol.name))

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
        logger.warning(('Minimization took too long'
                        ' and was terminated '
                        'by force on "{}".').format(macro_mol.name))
        _kill_bmin(macro_mol, macromodel_path)
        proc_out = ""

    logger.debug('Output of bmin on "{}" was: {}.'.format(
                                             macro_mol.name, proc_out))

    with open(log_file, 'r') as log:
        log_content = log.read()

    # Check the log for error reports.
    if ("termination due to error condition           21-" in
       log_content):
        raise _OptimizationError(("`bmin` crashed due to"
                                  " an error condition. "
                                  "See .log file."))

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
        return _run_bmin(macro_mol, macromodel_path, logger, timeout)

    # Make sure the .maegz file created by the optimization is present.
    maegz = file_root + '-out.maegz'
    _wait_for_file(maegz, logger)
    if not os.path.exists(log_file) or not os.path.exists(maegz):
        raise _OptimizationError(('The .log and/or .maegz '
                                  'files were not created by '
                                  'the optimization.'))


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
        output = sp.run(cmd, stdout=sp.PIPE, stderr=sp.STDOUT,
                        universal_newlines=True).stdout
        if time.time() - start > 600:
            break


def _run_applyhtreat(macro_mol, macromodel_path, logger):
    name, ext = os.path.splitext(macro_mol._file)
    mae = name + '.mae'
    mae_out = name + '_htreated.mae'
    _create_mae(macro_mol, macromodel_path, logger)

    app = os.path.join(macromodel_path, 'utilities', 'applyhtreat')
    cmd = [app, mae, mae_out]
    out = sp.run(cmd, stdout=sp.PIPE,
                 stderr=sp.STDOUT, universal_newlines=True)

    # If no license if found, keep re-running the function until it is.
    if not _license_found(out.stdout):
        return _run_applyhtreat(macro_mol, macromodel_path, logger)

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


def _generate_com(macro_mol, settings, logger):
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

    logger : FakeLogger or logging.Logger
        Used for logging.

    Modifies
    --------
    This function creates a new ``.com`` file holding the instructions
    for optimizing the macromolecule using MacroModel.

    Returns
    -------
    None : NoneType

    """

    logger.debug('Creating .com file for "{}".'.format(macro_mol.name))

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


def _generate_md_com(macro_mol, settings, logger):
    logger.debug('Creating .com file for "{}".'.format(macro_mol.name))

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


def _create_mae(macro_mol, macromodel_path, logger):
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

    logger : FakeLogger or logging.Logger
        Used for logging.

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

    logger.debug('Converting {} of "{}" to .mae.'.format(
                                                  ext, macro_mol.name))

    # Create the name of the new ``.mae`` file. It is the same as the
    # original structure file, including the same path. Only the
    # extensions are different.
    mae_file = macro_mol._file.replace(ext, '.mae')
    _structconvert(macro_mol._file, mae_file, macromodel_path, logger)
    return mae_file


def _convert_maegz_to_mae(macro_mol, macromodel_path, logger):
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

    logger : FakeLogger or logging.Logger
        Used for logging.

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

    logger.debug('Converting .maegz of "{}" to .mae.'.format(
                                                       macro_mol.name))
    name, ext = os.path.splitext(macro_mol._file)
    # ``out`` is the full path of the optimized ``.mae`` file.
    maegz = name + '-out.maegz'
    # Replace extensions to get the names of the various files.
    mae = name + '.mae'
    return _structconvert(maegz, mae, macromodel_path, logger)


def _structconvert(iname, oname, macromodel_path, logger):

    convrt_app = os.path.join(macromodel_path, 'utilities',
                              'structconvert')
    convrt_cmd = [convrt_app, iname, oname]

    # Execute the file conversion.
    try:
        convrt_return = sp.run(convrt_cmd,
                               stdout=sp.PIPE,
                               stderr=sp.STDOUT,
                               universal_newlines=True)

    # If conversion fails because a wrong Schrodinger path was given,
    # raise.
    except FileNotFoundError:
        raise _PathError(('Wrong Schrodinger path supplied to'
                          ' `structconvert` function.'))

    # If no license if found, keep re-running the function until it is.
    if not _license_found(convrt_return.stdout):
        return _structconvert(iname, oname, macromodel_path, logger)

    # If force field failed, raise.
    if 'number 1' in convrt_return.stdout:
        raise _ForceFieldError(convrt_return.stdout)

    _wait_for_file(oname, logger)
    if not os.path.exists(oname):
        raise _ConversionError(
         ('Conversion output file {} was not found.'
          ' Console output was {}.').format(oname,
                                            convrt_return.stdout))

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

    # Go through all the bonds in the rdkit molecule. If the bond
    # is not between bonder atoms add a fix line to the ``fix_block``.
    # If the bond does invovle two bonder atoms go to the next bond.
    # This is because a bond between 2 bonder atoms was added during
    # assembly and should therefore not be fixed.
    for bond in macro_mol.mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()

        if (atom1.GetIdx() in macro_mol.bonder_ids and
           atom2.GetIdx() in macro_mol.bonder_ids):
            continue

        # Make sure that the indices are increased by 1 in the .mae
        # file.
        atom1_id = atom1.GetIdx() + 1
        atom2_id = atom2.GetIdx() + 1
        fix_block += (_com_line('FXDI', atom1_id, atom2_id, 0, 0,
                                99999, 0, 0, 0) + '\n')

    return fix_block


def _fix_bond_angle_in_com_file(macro_mol, fix_block):
    """
    Adds lines fixing bond angles to ``.com`` body string.

    All bond angles of the molecule are fixed.

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

    # Create a substructure consisting of 3 dummy atoms bonded with 3
    # dummy bonds. This substructure will match with any 3 atoms which
    # are bonded together with any combination of bonds. These 3 atoms
    # will therefore have a bond angle.
    ba_mol = rdkit.MolFromSmarts('[*]~[*]~[*]')

    # Get the indices of all atoms which have a bond angle.
    # ``ba_atoms`` is a tuple of tuples of the form ((1,2,3), (4,5,6),
    # (7,8,9), ...). Each inner tuple holds the indicies of the atoms
    # which form a bond angle.
    ba_atoms = macro_mol.mol.GetSubstructMatches(ba_mol)

    for atom_ids in ba_atoms:
        atom_ids = [i+1 for i in atom_ids]
        fix_block += (_com_line('FXBA', *atom_ids,
                                99999, 0, 0, 0, 0) + '\n')

    return fix_block


def _fix_torsional_angle_in_com_file(macro_mol, fix_block):
    """
    Adds lines fixing torsional bond angles to ``.com`` body string.

    All torsional angles of the macromolecule are fixed.

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

    # Create a substructure consisting of 4 dummy atoms bonded with 4
    # dummy bonds. This substructure will match with any 4 atoms which
    # are bonded together with any combination of bonds. These 4 atoms
    # will therefore have a torsinal angle.
    ta_mol = rdkit.MolFromSmarts('[*]~[*]~[*]~[*]')

    # Get the indices of all atoms which have a torsional angle.
    # ``ta_atoms`` as a tuple of tuples of the form ((1,2,3,4),
    # (4,5,6,7), ...). Each inner tuple holds the indicies of the atoms
    # which have a torsional angle.
    ta_atoms = macro_mol.mol.GetSubstructMatches(ta_mol)

    # Apply the fix.
    for atom_ids in ta_atoms:
        atom_ids = [i+1 for i in atom_ids]
        fix_block += (_com_line('FXTA', *atom_ids,
                                99999, 361, 0, 0) + '\n')

    return fix_block


def _wait_for_file(file_name, logger, timeout=10):
    """
    Stalls until a given file exists or `timeout` expires.

    Parameters
    ----------
    file_name : str
        The full path of the file which should be waited for.

    logger : FakeLogger or logging.Logger
        Used for logging.

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
            logger.warning('Waiting for "{}".'.format(file_name))
            tick += 1

        if os.path.exists(file_name) or time_taken > timeout:
            break
