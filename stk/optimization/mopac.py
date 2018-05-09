"""
Defines optimization functions which use MOPAC.

"""

import os
import subprocess as sp
import psutil
import time
import logging
import rdkit.Chem.AllChem as rdkit
from uuid import uuid4

logger = logging.getLogger(__name__)


def mopac_opt(mol, mopac_path, settings=None):
    """
    Optimizes the molecule using MOPAC.

    This function runs an optimization. It is possible to provide
    different options, which correspond to the input keywords from
    MOPAC:

    http://openmopac.net/Manual/index.html

    Parameters
    ----------
    mol : :class:`.Molecule`
        The molecule to be optimized.

    mopac_path : :class:`str`
        The full path to the MOPAC suite on the user's
        machine. For example, in a default MacOS installation the
        folder will probably be something like
        ``/opt/mopac/MOPAC2016.exe``.

    settings: :class:`dict`, optional
        A dictionary which maps the names of the optimization
        parameters to their values. Valid values are:

            'hamiltonian' : :class:`str` (default = ``'PM7'``)
                A series of different methods can be selected:
                PM7, PM6, AM1, CIS (CISD, CISDT), MNDO, RM1, etc..

                PM7 is the latest version of the reparametrization of
                NDDO theory, where all the atomic and diatomic
                parameters were re-optimized / updated from PM6 [#]_.

            'method' : :class:`str` (default = ``'OPT'``)
                The default calculation consists of a geometry
                optimization. You can run single point calculations
                (SCF) or transition search algorithms (TS). Refer to
                the MOPAC website for specific keywords.

            'gradient' : :class:`float` (default = ``0.01``)
                The gradient at which the geometry optimization reaches
                the convergence criteria (``kcal`` / ``mol`` /
                ``Angstrom``). For small system and high precision
                work, ``0.01`` is recommended.

            'eps' : :class:`float` (default = ``80.1``)
                Sets the dielectric constant for the solvent. Presence
                of this keyword will cause the COSMO (Conductor-like
                Screening Model) method to be used to approximate the
                effect of a solvent model surrounding the molecule.
                Solvents with a low dielectric constant are not likely
                to work well with this model. ``0`` means that the
                dielectric constant is not included in the calculation.
                ``80.1`` can be used to model a water environment at
                room temperature.

            'charge' : :class:`float` (default = ``0``)
                The charge of the system.

            'fileout' : :class:`str` (default = ``'PDBOUT'``)
                Determines the output file type.

            'timeout' : :class:`float` (default = ``172800``)
                The amount in seconds the optimization is allowed to
                run before being terminated. The default value is ``2``
                days or ``172,800`` seconds. ``None`` means there
                is no timeout.

    Returns
    -------
    None : :class:`NoneType`

    References
    ----------
    .. [#] http://openmopac.net/PM7_accuracy/PM7_accuracy.html

    """

    if settings is None:
        settings = {}

    vals = {
            'hamiltonian': 'PM7',
            'method': 'OPT',
            'gradient': 0.01,
            'eps': 80.1,
            'charge': 0,
            'fileout': 'PDBOUT',
            'timeout': 172800,
            }
    vals.update(settings)

    mol._file = '{}.mol'.format(uuid4().int)

    # First write a .mol file of the molecule.
    mol.write(mol._file)
    # MOPAC requires a ``.mop`` file as input. This creates a ``.mop``
    # file holding the molecule.
    _create_mop(mol, vals)
    # Run the optimization
    _run_mopac(mol, mopac_path, settings)
    # Update the rdkit mol info with the ``.pdb`` file generated from
    # the MOPAC run
    _convert_mopout_to_mol(mol)


def _run_mopac(mol, mopac_path, settings, timeout=7200):

    name, ext = os.path.splitext(mol._file)
    mop_file = name + '.mop'

    print("", time.ctime(time.time()),
          'Running MOPAC - {}.'.format(mol.name), sep='\n')

    # To run MOPAC a command is issued to the console via
    # ``subprocess.Popen``. The command is the full path of the
    # ``mopac`` program.
    file_root, ext = os.path.splitext(mop_file)
    opt_cmd = [mopac_path, file_root]
    opt_proc = psutil.Popen(opt_cmd, stdout=sp.PIPE,
                            stderr=sp.STDOUT,
                            universal_newlines=True)

    try:
        proc_out, _ = opt_proc.communicate(timeout=timeout)
    except sp.TimeoutExpired:
        logger.warning('\nMinimization took too long and was terminated '
                       'by force - {}\n'.format(mol.name))
        _kill_mopac(mol)

    return


def _kill_mopac(mol):
    """
    Kills an in-progress MOPAC run.

    To kill a MOPAC run, a file with the molecule's name and a ``.end``
    extension is written.

    Parameters
    ----------
    mol : :class:`.Molecule`
        The molecule being optimized.

    Returns
    -------
    None : :class:`NoneType`

    """
    name, ext = os.path.splitext(mol._file)
    end_file = name + '.end'

    with open(end_file, 'w') as end:
        end.write('SHUT')


def _mop_line(settings):
    """
    Formats `settings` into a MOPAC input string.

    Parameters
    ----------
    settings : :class:`dict`
        Dictionary defined in :func:`mopac_opt`, where all the run
        details are defined.

    Returns
    -------
    :class:`str`
        String containing all the MOPAC keywords correctly formatted
        for the input file.

    """

    # Generate an empty string
    mopac_run_str = ""

    # Add Hamiltonian info
    mopac_run_str = mopac_run_str + settings['hamiltonian']
    # Add method and GNORM if 'OPT' otherwise just print the method
    if settings['method'] == 'OPT':
        gnorm_info = ' GNORM={} '.format(settings['gradient'])
        mopac_run_str = mopac_run_str + " " + settings['method'] + gnorm_info
    else:
        mopac_run_str = mopac_run_str + settings['method']
    # Add EPS info
    eps_info = ' EPS={} '.format(settings['eps'])
    mopac_run_str = mopac_run_str + eps_info
    # Add Charge info
    charge_info = ' CHARGE={} '.format(settings['charge'])
    mopac_run_str = mopac_run_str + charge_info
    # Add fileout info
    mopac_run_str = mopac_run_str + " " + settings['fileout']

    # Add the let keyword avoiding the crash of MOPAC
    mopac_run_str = mopac_run_str + " LET "

    return mopac_run_str


def _create_mop(mol, settings):
    """
    Creates the ``.mop`` file holding the molecule to be optimized.

    Parameters
    ----------
    mol : :class:`.Molecule`
        The molecule which is to be optimized. Its molecular
        structure file is converted to a ``.mop`` file. The original
        file is also kept.

    settings : :class:`dict`
        Dictionary defined in :func:`mopac_opt`, where all the run
        details are defined.

    Returns
    -------
    :class:`str`
        The full path of the newly created ``.mop`` file.

    """
    name, ext = os.path.splitext(mol._file)
    mop_file = name + '.mop'
    mol = mol.mol

    logger.info('\nCreating .mop file - {}.'.format(mol.name))

    # Generate the mop file containing the MOPAC run info
    with open(mop_file, 'w') as mop:
        # line for the run info
        mop.write(_mop_line(settings) + "\n")
        # line with the name of the molecule
        mop.write(name + "\n\n")

        # print the structural info
        for atom in mol.GetAtoms():
            atom_id = atom.GetIdx()
            atom_symbol = atom.GetSymbol()
            x, y, z = mol.GetConformer().GetAtomPosition(atom_id)
            atom_info = f"{atom_symbol}   {x}   +1  {y}   +1  {z}   +1 \n"
            mop.write(atom_info)

    return mop_file


def _convert_mopout_to_mol(mol):
    """
    Updates the molecular structure if the optimization is successful.

    Takes the ``.pdb`` file of the neutral file generated from the
    MOPAC run and initializes a new ``rdkit`` molecule with those
    coordinates. `mol` is then updated to hold the new molecule.

    Parameters
    ----------
    mol : :class:`.Molecule`
        The macromolecule being optimized. The ``.pdb`` file holding
        its optimized structure is converted to a rdkit molecule.

    Returns
    -------
    None : :class:`NoneType`

    """
    name, ext = os.path.splitext(mol._file)
    pdb_file = name + ".pdb"

    logger.info("\nUpdating molecule with MOPAC optimized "
                "one - {}.\n".format(mol.name))

    new_mol = rdkit.MolFromPDBFile(pdb_file,
                                   sanitize=False,
                                   removeHs=False)
    # Updating the macro_mol.mol infos with the new mol
    mol.mol = new_mol
