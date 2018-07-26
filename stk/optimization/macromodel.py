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
import gzip

from ..utilities import MAEExtractor, flatten


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


def macromodel_opt(mol,
                   macromodel_path,
                   settings=None,
                   md=None,
                   conformer=-1):
    """
    Optimizes the molecule using MacroModel.

    This function runs a restricted optimization. The structures of the
    building blocks are frozen and only the new bonds formed between
    building blocks during assembly are optimized.

    Parameters
    ----------
    mol : :class:`.Molecule`
        The molecule who's structure must be optimized.

    macromodel_path : :class:`str`
        The full path of the Schrodinger suite within the user's
        machine. For example, on a Linux machine this may be something
        like ``'/opt/schrodinger2017-2'``.

    settings : :class:`dict`, optional
        A dictionary which maps the names of optimization parameters to
        their values. Valid values are:

            'restricted' : :class:`bool` (default = ``True``)
                If ``False`` then all bonds are optimized, not just the
                ones created during macromolecular assembly. If
                ``True`` then an optimization is performed only on the
                bonds added during molecular assembly. If
                ``'both'`` then a restricted optimization is performed
                first, followed by a regular optimization.

            'timeout' : :class:`float` (default = ``None``)
                The amount in seconds the optimization is allowed to
                run before being terminated. ``None`` means there is no
                timeout.

            'force_field' : :class:`int` (default = ``16``)
                The number of the force field to be used.

            'max_iter' : :class:`int` (default = ``2500``)
                The maximum number of iterations done during the
                optimization.

            'gradient' : :class:`float` (default = ``0.05``)
                The gradient at which optimization is stopped.

            'md' : :class:`bool` (default = ``False``)
                Toggles whether a MD conformer search should be
                performed.

        Only values which need to be changed from the default need to
        be specified. For exmaple,

        .. code-block:: python

            macromodel_opt(mol, 'path', {'max_iter' : 10})


    md : :class:`dict`, optional
        A dictionary holding settings for the MD conformer search.
        This parameter is used in the same way as `settings` except
        the values effect the MD only. See docstring of
        :func:`_macromodel_md_opt` for valid values.

    conformer : :class:`int`, optional
        The id of the conformer to be optimized.

    Returns
    -------
    None : :class:`NoneType`

    """

    if settings is None:
        settings = {}
    if md is None:
        md = {}

    vals = {
             'restricted': True,
             'timeout': None,
             'force_field': 16,
             'max_iter': 2500,
             'gradient': 0.05,
             'md': False
            }
    vals.update(settings)

    try:
        mol._file = '{}.mol'.format(uuid4().int)
        # First write a .mol file of the molecule.
        mol.write(mol._file, conformer)
        # MacroModel requires a ``.mae`` file as input. This creates a
        # ``.mae`` file holding the molecule.
        _create_mae(mol, macromodel_path)
        # generate the ``.com`` file for the MacroModel run.
        _generate_com(mol, vals)
        # Run the optimization.
        _run_bmin(mol, macromodel_path, vals['timeout'])
        # Get the ``.maegz`` file output from the optimization and
        # convert it to a ``.mae`` file.
        _convert_maegz_to_mae(mol)
        mol.update_from_mae(mol._file.replace('.mol', '.mae'),
                            conformer)

        if vals['restricted'] == 'both':
            new_vals = dict(vals)
            new_vals['md'] = False
            new_vals['restricted'] = False
            macromodel_opt(mol=mol,
                           macromodel_path=macromodel_path,
                           settings=new_vals,
                           md={})

        if vals['md']:
            _macromodel_md_opt(mol,
                               macromodel_path,
                               md,
                               conformer)

    except _ForceFieldError as ex:
        # If OPLS_2005 has been tried already - record an exception.
        if vals['force_field'] == 14:
            raise ex

        # If OPLSE_2005 has not been tried - try it.
        logger.warning(('Minimization with OPLS3 failed on "{}". '
                        'Trying OPLS_2005.').format(mol.name))

        vals['force_field'] = 14
        return macromodel_opt(mol,
                              macromodel_path,
                              vals,
                              md,
                              conformer)


def macromodel_cage_opt(mol,
                        macromodel_path,
                        settings=None,
                        md=None,
                        conformer=-1):
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
    mol : :class:`.Molecule`
        The molecule who's structure must be optimized.

    macromodel_path : :class:`str`
        The full path of the Schrodinger suite within the user's
        machine. For example, on a Linux machine this may be something
        like ``'/opt/schrodinger2017-2'``.

    settings : :class:`dict`, optional
        A dictionary which maps the names of optimization parameters to
        their values. Valid values are:

            'restricted' : :class:`bool` (default = ``True``)
                If ``False`` then all bonds are optimized, not just the
                ones created during macromolecular assembly. If
                ``True`` then an optimization is performed only on the
                bonds added during molecular assembly. If
                ``'both'`` then a restricted optimization is performed
                first, followed by a regular optimization.

            'timeout' : :class:`float` (default = ``None``)
                The amount in seconds the optimization is allowed to
                run before being terminated. ``None`` means there is no
                timeout.

            'force_field' : :class:`int` (default = ``16``)
                The number of the force field to be used.

            'max_iter' : :class:`int` (default = ``2500``)
                The maximum number of iterations done during the
                optimization.

            'gradient' : :class:`float` (default = ``0.05``)
                The gradient at which optimization is stopped.

            'md' : :class:`bool` (default = ``False``)
                Toggles whether a MD conformer search should be
                performed.

        Only values which need to be changed from the default need to
        be specified. For exmaple,

        .. code-block:: python

            macromodel_opt(mol, 'path', {'max_iter' : 10})


    md : :class:`dict`, optional
        A dictionary holding settings for the MD conformer search.
        This parameter is used in the same way as `settings` except
        the values effect the MD only. See docstring of
        :func:`_macromodel_md_opt` for valid values.

    conformer : :class:`int`, optional
        The id of the conformer to be optimized.

    Returns
    -------
    None : :class:`NoneType`

    """

    if settings is None:
        settings = {}
    if md is None:
        md = {}

    vals = {
             'restricted': True,
             'timeout': None,
             'force_field': 16,
             'max_iter': 2500,
             'gradient': 0.05,
             'md': False
            }
    vals.update(settings)

    try:
        mol._file = '{}.mol'.format(uuid4().int)
        # First write a .mol file of the molecule.
        mol.write(mol._file, conformer)
        # MacroModel requires a ``.mae`` file as input. This creates a
        # ``.mae`` file holding the molecule.
        _create_mae(mol, macromodel_path)
        # generate the ``.com`` file for the MacroModel run.
        _generate_com(mol, vals)
        # Run the optimization.
        _run_bmin(mol, macromodel_path, vals['timeout'])
        # Get the ``.maegz`` file output from the optimization and
        # convert it to a ``.mae`` file.
        _convert_maegz_to_mae(mol)
        mol.update_from_mae(mol._file.replace('.mol', '.mae'),
                            conformer)

        if vals['restricted'] == 'both':
            new_vals = dict(vals)
            new_vals['md'] = False
            new_vals['restricted'] = False
            macromodel_opt(mol=mol,
                           macromodel_path=macromodel_path,
                           settings=new_vals,
                           md={})

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            windows = mol.windows(conformer)
            if windows is not None:
                all_windows = (len(windows) ==
                               mol.topology.n_windows)

                if vals['md'] and all_windows:
                    _macromodel_md_opt(mol,
                                       macromodel_path,
                                       md,
                                       conformer)

    except _ForceFieldError as ex:
        # If OPLS_2005 has been tried already - record an exception.
        if vals['force_field'] == 14:
            raise ex

        # If OPLSE_2005 has not been tried - try it.
        logger.warning(('Minimization with OPLS3 failed on "{}". '
                        'Trying OPLS_2005.').format(mol.name))

        vals['force_field'] = 14
        return macromodel_cage_opt(mol,
                                   macromodel_path,
                                   vals,
                                   md,
                                   conformer)


def _macromodel_md_opt(mol,
                       macromodel_path,
                       settings=None,
                       conformer=-1):
    """
    Runs a MD conformer search on `mol`.

    Parameters
    ----------
    mol : :class:`.Molecule`
        The molecule who's structure must be optimized.

    macromodel_path : :class:`str`
        The full path of the Schrodinger suite within the user's
        machine. For example, on a Linux machine this may be something
        like ``'/opt/schrodinger2017-2'``.

    settings : :class:`dict`, optional
        Each key is a string describing an MD parameter. The value is
        the corresponding value. Only values which need to be changed
        from the default need to be specified. The allowed parameters
        are:

            'timeout' : :class:`float` (default = ``None``)
                The amount in seconds the MD is allowed to run before
                being terminated. ``None`` means there is no timeout.

            'force_field' : :class:`int` (default = ``16``)
                The number of the force field to be used.

            'temp' : :class:`float` (default = ``300``)
                The temperature in Kelvin at which the MD is run.

            'confs' : :class:`int` (default = ``50``)
                The number of conformers sampled and optimized from the
                MD.

            'time_step' : :class:`float` (default = ``1.0``)
                The time step in ``fs`` for the MD.

            'eq_time' : :class:`float` (default = ``10``)
                The equilibriation time in ``ps`` before the MD is run.

            'sim_time' : :class:`float` (default = ``200``)
                The simulation time in ``ps`` of the MD.

            'max_iter' : :class:`int` (default = ``2500``)
                The maximum number of iterations done during the
                optimization.

            'gradient' : float (default = ``0.05``)
                The gradient at which optimization is stopped.

    conformer : :class:`int`, optional
        The id of the conformer to be optimized.

    Returns
    -------
    None : :class:`NoneType`

    """

    if settings is None:
        settings = {}

    vals = {
               'timeout': None,
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

    logger.info('Running MD on "{}".'.format(mol.name))
    try:
        mol._file = '{}.mol'.format(uuid4().int)
        # First write a .mol file of the molecule.
        mol.write(mol._file, conformer)
        # MacroModel requires a ``.mae`` file as input. This creates a
        # ``.mae`` file holding the molecule.
        _create_mae(mol, macromodel_path)
        # Generate the ``.com`` file for the MacroModel MD run.
        _generate_md_com(mol, vals)
        # Run the optimization.
        _run_bmin(mol, macromodel_path, vals['timeout'])
        # Extract the lowest energy conformer into its own .mae file.
        conformer_mae = MAEExtractor(mol._file).path
        mol.update_from_mae(conformer_mae, conformer)

    except _ForceFieldError as ex:
        # If OPLS_2005 has been tried already - record an exception.
        if vals['force_field'] == 14:
            raise ex
        # If OPLSE_2005 has not been tried - try it.
        logger.warning(('Minimization with OPLS3 failed on "{}". '
                        'Trying OPLS_2005.').format(mol.name))

        vals['force_field'] = 14
        return _macromodel_md_opt(mol,
                                  macromodel_path,
                                  vals,
                                  conformer)


def _run_bmin(macro_mol, macromodel_path, timeout):

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

    incomplete = True
    while incomplete:
        opt_proc = psutil.Popen(opt_cmd,
                                stdout=sp.PIPE,
                                stderr=sp.STDOUT,
                                universal_newlines=True)
        try:
            proc_out, _ = opt_proc.communicate(timeout=timeout)

        except sp.TimeoutExpired:
            logger.warning(('Minimization took too long'
                            ' and was terminated '
                            f'by force on "{macro_mol.name}".'))
            _kill_bmin(macro_mol, macromodel_path)
            proc_out = ""

        logger.debug(
            f'Output of bmin on "{macro_mol.name}" was: {proc_out}.')

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
             " structure") in log_content and
            ("skipping input structure  due to "
             "forcefield interaction errors") in log_content):
            raise _LewisStructureError(
                    '`bmin` failed due to poor Lewis structure.')

        # If optimization fails because a wrong Schrodinger path was
        # given, raise.
        if 'The system cannot find the path specified' in proc_out:
            raise _PathError(('Wrong Schrodinger path supplied to'
                              ' `macromodel_opt` function.'))

        # If optimization fails because the license is not found, rerun
        # the function.
        if _license_found(proc_out, macro_mol):
            incomplete = False

    # Make sure the .maegz file created by the optimization is present.
    maegz = file_root + '-out.maegz'
    _wait_for_file(maegz)
    if not os.path.exists(log_file) or not os.path.exists(maegz):
        raise _OptimizationError(('The .log and/or .maegz '
                                  'files were not created by '
                                  'the optimization.'))


def _kill_bmin(macro_mol, macromodel_path):
    name, ext = os.path.splitext(macro_mol._file)
    name = re.split(r'\\|/', name)[-1]
    app = os.path.join(macromodel_path, 'jobcontrol')
    cmd = [app, '-stop', name]

    incomplete = True
    while incomplete:
        out = sp.run(cmd, stdout=sp.PIPE,
                     stderr=sp.STDOUT, universal_newlines=True)

        # If no license if found, keep re-running the function until it
        # is.
        if _license_found(out.stdout):
            incomplete = False

    # This loop causes the function to wait until the job has been
    # killed via job control. This means the output files will have
    # been written by the time the function exits. Essentially the loop
    # continues until the job is no longer found by
    # "./jobcontrol -list"
    cmd = [app, '-list']
    output = name
    start = time.time()
    while name in output:
        output = sp.run(cmd, stdout=sp.PIPE, stderr=sp.STDOUT,
                        universal_newlines=True).stdout
        if time.time() - start > 600:
            break


def _license_found(output, mol=None):
    """
    Checks to see if minimization failed due to a missing license.

    The user can be notified of this in one of two ways. Sometimes the
    output of the submission contains the message informing that the
    license was not found and in other cases it will be the log file.
    This function checks both of these sources for this message.

    Parameters
    ----------
    output : :class:`str`
        The output from submitting the minimization of the structure
        to the ``bmin`` program.

    mol : :class:`.Molecule`, optional
        The molecule being optimized. If the ``.log`` file is not to
        be checked, the default ``None`` should be used.

    Returns
    -------
    :class:`bool`
        ``True`` if the license was found. ``False`` if the
        minimization did not occur due to a missing license.

    """

    if 'Could not check out a license for mmlibs' in output:
        return False
    if mol is None:
        return True

    # To check if the log file mentions a missing license file open the
    # the log file and scan for the apporpriate string.

    # Check if the file exists first. If not, this is often means the
    # calculation must be redone so return False anyway.
    log_file_path = mol._file.replace('mol', 'log')
    with open(log_file_path, 'r') as log_file:
        log_file_content = log_file.read()

    if 'Could not check out a license for mmlibs' in log_file_content:
        return False

    return True


def _com_line(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9):
    return (" {:<5}{:>7}{:>7}{:>7}{:>7}{:>11.4f}{:>11.4f}"
            "{:>11.4f}{:>11.4f}").format(arg1, arg2, arg3, arg4,
                                         arg5, arg6, arg7, arg8, arg9)


def _generate_com(mol, settings):
    """
    Create a ``.com`` file for a MacroModel optimization.

    The created ``.com`` file fixes all bond parameters which were not
    added during assembly. This means all bond distances, bond angles
    and torsional angles are fixed, except for cases where it involves
    a bond added during assembly of the macromolecule.

    This fixing is implemented by creating a ``.com`` file with various
    "FX" commands written within its body.

    Parameters
    ----------
    mol : :class:`.Molecule`
        The molecule which is to be optimized.

    settings : :class:`dict`
        A dictionary of settings for the optimization. See
        :func:`macromodel_opt` documentation.

    Returns
    -------
    None : :class:`NoneType`

    """

    logger.debug('Creating .com file for "{}".'.format(mol.name))

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
    name, ext = os.path.splitext(mol._file)
    com_file = name + '.com'
    mae = name + '.mae'
    output = name + '-out.maegz'

    # This function adds all the lines which fix bond distances and
    # angles into ``main_string``.
    main_string = _fix_params_in_com_file(mol, main_string,
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


def _create_mae(mol, macromodel_path):
    """
    Creates the ``.mae`` file holding the molecule to be optimized.

    Parameters
    ----------
    mol : :class:`.Molecule`
        The molecule which is to be optimized. Its molecular
        structure file is converted to a ``.mae`` file. The original
        file is also kept.

    macromodel_path : :class:`str`
        The full path of the Schrodinger suite within the user's
        machine. For example, on a Linux machine this may be something
        like ``'/opt/schrodinger2017-2'``.

    Returns
    -------
    :class:`str`
        The full path of the newly created ``.mae`` file.

    """

    _, ext = os.path.splitext(mol._file)

    logger.debug(f'Converting {ext} of "{mol.name}" to .mae.')

    # Create the name of the new ``.mae`` file. It is the same as the
    # original structure file, including the same path. Only the
    # extensions are different.
    mae_file = mol._file.replace(ext, '.mae')
    _structconvert(mol._file, mae_file, macromodel_path)
    return mae_file


def _convert_maegz_to_mae(mol):
    """
    Converts a ``.maegz`` file to a ``.mae`` file.

    Parameters
    ----------
    mol : :class:`.Molecule`
        The molecule being optimized. The ``.maegz`` file holding
        its optimized structure is converted to a ``.mae`` file. Both
        versions are kept.

    Returns
    -------
    None : :class:`NoneType`

    """

    logger.debug(f'Converting .maegz of "{mol.name}" to .mae.')
    name, ext = os.path.splitext(mol._file)
    # ``out`` is the full path of the optimized ``.mae`` file.
    maegz = name + '-out.maegz'
    # Replace extensions to get the names of the various files.
    mae = name + '.mae'

    gz_file = gzip.open(maegz)
    with open(mae, 'wb') as f:
        f.write(gz_file.read())
    gz_file.close()


def _structconvert(iname, oname, macromodel_path):

    convrt_app = os.path.join(macromodel_path,
                              'utilities',
                              'structconvert')
    convrt_cmd = [convrt_app, iname, oname]

    incomplete = True
    while incomplete:

        # Execute the file conversion.
        try:
            convrt_return = sp.run(convrt_cmd,
                                   stdout=sp.PIPE,
                                   stderr=sp.STDOUT,
                                   universal_newlines=True)

        # If conversion fails because a wrong Schrodinger path was
        # given, raise.
        except FileNotFoundError:
            raise _PathError(('Wrong Schrodinger path supplied to'
                              ' `structconvert` function.'))

        if 'File does not exist' in convrt_return.stdout:
            raise _ConversionError(
                    (f'structconvert input file, {iname}, missing. '
                     f'Console output was {convrt_return.stdout}'))

        # If no license if found, keep re-running the function until it
        # is.
        if _license_found(convrt_return.stdout):
            incomplete = False

    # If force field failed, raise.
    if 'number 1' in convrt_return.stdout:
        raise _ForceFieldError(convrt_return.stdout)

    _wait_for_file(oname)
    if not os.path.exists(oname):
        raise _ConversionError(
         (f'Conversion output file {oname} was not found.'
          f' Console output was {convrt_return.stdout}.'))

    return convrt_return


def _fix_params_in_com_file(mol, main_string, restricted):
    """
    Adds lines to the ``.com`` body fixing bond distances and angles.

    For each bond distance, bond angle and torisional angle that does
    not involve a bond created during assembly a "FX" command is
    added to the string holding holding the body of the ``.com`` file.

    These lines replace the filler line in the main string.

    Parameters
    ----------
    mol : :class:`.Molecule`
        The molecule which is to be optimized.

    main_string : :class:`str`
        The body of the ``.com`` file which is to have fix commands
        added.

    restricted : :class:`bool`
        When ``False`` the block containing instructions to fix
        molecular parameters is not added to the ``.com`` file.

    Returns
    -------
    :class:`str`
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
    fix_block = _fix_distance_in_com_file(mol, fix_block)
    # Add lines that fix the bond angles.
    fix_block = _fix_bond_angle_in_com_file(mol, fix_block)
    # Add lines that fix the torsional angles.
    fix_block = _fix_torsional_angle_in_com_file(mol, fix_block)

    return main_string.replace(("!!!BLOCK_OF_FIXED_PARAMETERS_"
                                "COMES_HERE!!!\n"), fix_block)


def _fix_distance_in_com_file(mol, fix_block):
    """
    Adds lines fixing bond distances to ``.com`` body string.

    Only bond distances which do not involve bonds created during
    assembly are fixed.

    Parameters
    ----------
    mol : :class:`.Molecule`
        The molecule to be optimized.

    fix_block : :class:`str`
        The string holding all the lines containing fix commands for
        the ``.com`` file.

    Returns
    -------
    :class:`str`
        A string holding lines containg fix commands for the ``.com``
        file. The string has the lines fixing bond distances added to
        it by this function.

    """

    bonder_ids = set(flatten(mol.bonder_ids))

    # Go through all the bonds in the rdkit molecule. If the bond
    # is not between bonder atoms add a fix line to the ``fix_block``.
    # If the bond does invovle two bonder atoms go to the next bond.
    # This is because a bond between 2 bonder atoms was added during
    # assembly and should therefore not be fixed.
    for bond in mol.mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()

        if (atom1.GetIdx() in bonder_ids and
           atom2.GetIdx() in bonder_ids):
            continue

        # Make sure that the indices are increased by 1 in the .mae
        # file.
        atom1_id = atom1.GetIdx() + 1
        atom2_id = atom2.GetIdx() + 1
        fix_block += (_com_line('FXDI', atom1_id, atom2_id, 0, 0,
                                99999, 0, 0, 0) + '\n')

    return fix_block


def _fix_bond_angle_in_com_file(mol, fix_block):
    """
    Adds lines fixing bond angles to ``.com`` body string.

    All bond angles of the molecule are fixed.

    Parameters
    ----------
    mol : :class:`.Molecule`
        The molecule to be optimized.

    fix_block : :class:`str`
        The string holding all the lines containing fix commands for
        the ``.com`` file.

    Returns
    -------
    :class:`str`
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
    ba_atoms = mol.mol.GetSubstructMatches(ba_mol)

    for atom_ids in ba_atoms:
        atom_ids = [i+1 for i in atom_ids]
        fix_block += (_com_line('FXBA', *atom_ids,
                                99999, 0, 0, 0, 0) + '\n')

    return fix_block


def _fix_torsional_angle_in_com_file(mol, fix_block):
    """
    Adds lines fixing torsional bond angles to ``.com`` body string.

    All torsional angles of the molecule are fixed.

    Parameters
    ----------
    macro_mol : :class:`.Molecule`
        The molecule to be optimized.

    fix_block : :class:`str`
        The string holding all the lines containing fix commands for
        the ``.com`` file.

    Returns
    -------
    :class:`str`
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
    ta_atoms = mol.mol.GetSubstructMatches(ta_mol)

    # Apply the fix.
    for atom_ids in ta_atoms:
        atom_ids = [i+1 for i in atom_ids]
        fix_block += (_com_line('FXTA', *atom_ids,
                                99999, 361, 0, 0) + '\n')

    return fix_block


def _wait_for_file(file_name, timeout=10):
    """
    Stalls until a given file exists or `timeout` expires.

    Parameters
    ----------
    file_name : :class:`str`
        The full path of the file which should be waited for.

    timeout : :class:`int` or :class:`float`, optional
        The number of seconds before the function stops waiting and
        returns.

    Returns
    --------
    None : :class:`NoneType`

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
