"""
Defines optimization functions which use MOPAC.
"""

import os
import subprocess as sp
import psutil
import time
import loggin
import rdkit.Chem.AllChem as rdkit
from uuid import uuid4

logger = logging.getLogger(__name__)

def mopac_opt(macro_mol, mopac_path, settings={}):
    """
    Optimizes the molecule using MOPAC.

    This function runs an optimization. It is possible to provide different
    options, which correspond to the input keywords from MOPAC:
    http://openmopac.net/Manual/index.html

    * need to create a tool for generating mopac inputs

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule who's structure must be optimized.

    mopac_path : str
        The full path of the ``MOPAC`` suite within the user's
        machine. For example, in a default MacOS installation the
        folder will probably be something like
        ``/opt/mopac/MOPAC2016.exe``.

    settings: dict (default = {})
        A dictionary which maps the names of the optimization parameters to
        their values. Valid values are:

            'hamiltonian' : string(default = 'PM7')
                A series of different methods can be selected:
                PM7, PM6, AM1, CIS (CISD, CISDT), MNDO, RM1, etc..

                PM7 is the latest version of the reparametrization of the NDDO
                theory, where all the atomic and diatomic parameters were
                re-optimized - update compared to PM6.
                http://openmopac.net/PM7_accuracy/PM7_accuracy.html

            'method' : string (default = 'OPT')
                The default calculation consists in a geometry optimization.
                You can run single point calculations (SCF) or transition
                search algorithms (TS). Refer to the MOPAC website for specific
                keywords.

            'gradient' : float (default = 0.01)
                The gradient at which the geometry optimization reaches the
                convergence criteria (kcal/mol/Angstrom). For small system high
                precision work, 0.01 is recommended, as these results are
                easily good enough for all high precision work.

            'eps' : float (detault = 80.1)
                Sets the dielectric constant for the solvent. Presence of this
                keyword will cause the COSMO (Conductor-like Screening Model)
                method to be used to approximate the effect of a solvent model
                surrounding the molecule. Solvents with low dielectric constant
                are not likely to work well with this model.
                0 means that the dielectric constant is not included in the
                calculation.
                80.1 can be used to model a water environment at room
                temperature.

            'charge' : list of floats (default = 0)
                When the system being studied is an ion, the charge, n, on the
                ion must be supplied as an integer. For cations n can be 1, 2,
                3, etc.; for anions -1, -2, -3, etc.


            'fileout' : string (default = 'PDBOUT')
                This generates the pdb file with the optimized structure.

            'timeout' : float (default = 172800)
                The amount in seconds the optimization is allowed to run before
                being terminated. The default value is 2 days =
                172,800 seconds.


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
            'hamiltonian': 'PM7',
            'method': 'OPT',
            'gradient': 0.01,
            'eps': 80.1,
            'charge': 0,
            'fileout': 'PDBOUT',
            'timeout': 172800,
            }
    vals.update(settings)

    macro_mol._file = getattr(macro_mol, '_file',
                              '{}.mol'.format(uuid4().int))

    # First write a .mol file of the molecule.
    macro_mol.write(macro_mol._file)
    # MOPAC requires a ``.mop`` file as input. This creates a ``.mop``
    # file holding the molecule.
    _create_mop(macro_mol, vals)
    # Run the optimization
    _run_mopac(macro_mol, mopac_path, settings)
    # Update the rdkit mol info with the ``.pdb`` file generated from
    # the MOPAC run
    _convert_mopout_to_mol(macro_mol)


def _run_mopac(macro_mol, mopac_path, settings, timeout=7200,
                                        logger=logger):

    name, ext = os.path.splitext(macro_mol._file)
    mop_file = name + '.mop'

    print("", time.ctime(time.time()),
          'Running MOPAC - {}.'.format(macro_mol.name), sep='\n')

    # To run MOPAC a command is issued to the console via
    # ``subprocess.Popen``. The command is the full path of the
    # ``mopac`` program.
    file_root, ext = os.path.splitext(mop_file)
    opt_cmd = [mopac_path, file_root]
    opt_proc = psutil.Popen(opt_cmd, stdout=sp.PIPE,
                            stderr=sp.STDOUT,
                            universal_newlines=True)

    try:
        if timeout:
            proc_out, _ = opt_proc.communicate(timeout=timeout)
        else:
            proc_out, _ = opt_proc.communicate()
    except sp.TimeoutExpired:
        logger.warning('\nMinimization took too long and was terminated '
                       'by force - {}\n').format(macro_mol.name))
        _kill_mopac(macro_mol)

    return


def _kill_mopac(macro_mol):
    """
    To kill a MOPAC run for a specific structure it is enough to generate
    a non empty file with the molecule's name with the `.end` extension.
    """
    name, ext = os.path.splitext(macro_mol._file)
    end_file = name + '.end'

    with open(end_file, 'w') as end:
        end.write('SHUT')


def _mop_line(settings):
    """
    Formats the settings dictionary with the correct keywords for MOPAC into
    a string to be added to the MOPAC input.

    Parameters
    ----------
    settings : dict
        Dictionary defined in the mopac_opt function, where all the run details
        are defined.

    Returns
    -------
    mopac_run_str : str
        String containing all the MOPAC keywords correctly formatted for the
        input file.
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


def _create_mop(macro_mol, settings, logger=logger):
    """
    Creates the ``.mop`` file holding the molecule to be optimized.
    The name of the input file will contain info about its charge:
    charge = 0: name_neu
    charge = -1: name_an1
    charge = +1: name_cat 1

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule which is to be optimized. Its molecular
        structure file is converted to a ``.mop`` file. The original
        file is also kept.

    mopac_run_str : str
        This string specifies the MOPAC keywords to be used in the input for
        the calculation.

    logger : FakeLogger or logging.Logger, optional
        Used for logging.

    Modifies
    --------
    This function creates a new ``.mop`` file from the structure file
    in `macro_mol._file`. This new file is placed in the same
    folder as the original file and has the same name with the _charge info.

    Returns
    -------
    str
        The full path of the newly created ``.mop`` file.
    """
    name, ext = os.path.splitext(macro_mol._file)
    mop_file = name + '.mop'
    mol = macro_mol.mol

    logger.log('\nCreating .mop file - {}.'.format(macro_mol.name))

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
            atom_info = "{}   {}   +1  {}   +1  {}   +1 \n".format(atom_symbol,
                                                                   x, y, z)
            mop.write(atom_info)

    return mop_file


def _convert_mopout_to_mol(macro_mol, logger=logger):
    """
    Updates the molecule information (coords) if the opt is successful.
    Takes the ``.pdb`` file of the neutral file generated from the MOPAC run
    and initiates a new rdkit molecule with that coordinates.
    The macro_mol instance is then updated with the new molecule.

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule being optimized. The ``.pdb`` file holding
        its optimized structure is converted to a rdkit molecule.

    logger : FakeLogger or logging.Logger, optional
        Used for logging.

    Modifies
    --------
    This function updates the macro_mol instance.

    Returns
    -------
    None : NoneType

    """
    name, ext = os.path.splitext(macro_mol._file)
    pdb_file = name + ".pdb"

    logger.log("\nUpdating molecule with MOPAC optimized "
               "one - {}.\n".format(macro_mol.name))

    new_mol = rdkit.MolFromPDBFile(pdb_file, sanitize=False,
                                   removeHs=False)
    # Updating the macro_mol.mol infos with the new mol
    macro_mol.mol = new_mol
