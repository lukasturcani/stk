"""
Defines MOPAC optimizers.

"""

import os
import subprocess as sp
import psutil
import logging
from uuid import uuid4

from .optimizers import Optimizer

logger = logging.getLogger(__name__)


class MOPAC(Optimizer):
    """
    Uses MOPAC to optimize molcules.

    Attributes
    ----------
    mopac_path : :class:`str`
        The full path to the MOPAC suite on the user's
        machine. For example, in a default MacOS installation the
        folder will probably be something like
        ``/opt/mopac/MOPAC2016.exe``.

    hamiltonian : :class:`str`, optional
        A series of different methods can be selected:
        ``'PM7'``, ``'PM6'``, ``'AM1'``, ``'CIS'``
        ``'(CISD, CISDT)'``, ``'MNDO'``, ``'RM1'``, etc.

        ``'PM7'`` is the latest version of the reparametrization of
        NDDO theory, where all the atomic and diatomic
        parameters were re-optimized / updated from PM6 [#]_.

    method : :class:`str`, optional
        The default calculation consists of a geometry
        optimization. You can run single point calculations
        (SCF) or transition search algorithms (TS). Refer to
        the MOPAC website for specific keywords.

    minimum_gradient : :class:`float`, optional
        The gradient at which the geometry optimization reaches
        the convergence criteria (``kcal`` / ``mol`` /
        ``Angstrom``). For small system and high precision
        work, ``0.01`` is recommended.

    eps : :class:`float`, optional
        Sets the dielectric constant for the solvent. Presence
        of this keyword will cause the COSMO (Conductor-like
        Screening Model) method to be used to approximate the
        effect of a solvent model surrounding the molecule.
        Solvents with a low dielectric constant are not likely
        to work well with this model. ``0`` means that the
        dielectric constant is not included in the calculation.
        ``80.1`` can be used to model a water environment at
        room temperature.

    charge : :class:`float`, optional
        The charge of the system.

    fileout : :class:`str`, optional
        Determines the output file type.

    timeout : :class:`float`, optional
        The amount in seconds the optimization is allowed to
        run before being terminated. The default value is ``2``
        days or ``172,800`` seconds. ``None`` means there
        is no timeout.

    References
    ----------
    .. [#] http://openmopac.net/PM7_accuracy/PM7_accuracy.html

    """

    def __init__(self,
                 mopac_path,
                 hamiltonian='PM7',
                 method='OPT',
                 minimum_gradient=0.01,
                 eps=80.1,
                 charge=0,
                 fileout='PDBOUT',
                 timeout=172800,
                 use_cache=False):
        """
        Initializes a :class:`MOPAC` instance.

        Parameters
        ----------
        mopac_path : :class:`str`
            The full path to the MOPAC suite on the user's
            machine. For example, in a default MacOS installation the
            folder will probably be something like
            ``/opt/mopac/MOPAC2016.exe``.

        hamiltonian : :class:`str`, optional
            A series of different methods can be selected:
            ``'PM7'``, ``'PM6'``, ``'AM1'``, ``'CIS'``
            ``'(CISD, CISDT)'``, ``'MNDO'``, ``'RM1'``, etc.

            ``'PM7'`` is the latest version of the reparametrization of
            NDDO theory, where all the atomic and diatomic
            parameters were re-optimized / updated from PM6 [#]_.

        method : :class:`str`, optional
            The default calculation consists of a geometry
            optimization. You can run single point calculations
            (SCF) or transition search algorithms (TS). Refer to
            the MOPAC website for specific keywords.

        minimum_gradient : :class:`float`, optional
            The gradient at which the geometry optimization reaches
            the convergence criteria (``kcal`` / ``mol`` /
            ``Angstrom``). For small system and high precision
            work, ``0.01`` is recommended.

        eps : :class:`float`, optional
            Sets the dielectric constant for the solvent. Presence
            of this keyword will cause the COSMO (Conductor-like
            Screening Model) method to be used to approximate the
            effect of a solvent model surrounding the molecule.
            Solvents with a low dielectric constant are not likely
            to work well with this model. ``0`` means that the
            dielectric constant is not included in the calculation.
            ``80.1`` can be used to model a water environment at
            room temperature.

        charge : :class:`float`, optional
            The charge of the system.

        fileout : :class:`str`, optional
            Determines the output file type.

        timeout : :class:`float`, optional
            The amount in seconds the optimization is allowed to
            run before being terminated. The default value is ``2``
            days or ``172,800`` seconds. ``None`` means there
            is no timeout.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule and conformer.

        References
        ----------
        .. [#] http://openmopac.net/PM7_accuracy/PM7_accuracy.html

        """

        self.mopac_path
        self.hamiltonian = hamiltonian
        self.method = method
        self.minimum_gradient = minimum_gradient
        self.eps = eps
        self.charge = charge
        self.fileout = fileout
        self.timeout = timeout
        super().__init__(use_cache=use_cache)

    def optimize(self, mol, conformer=-1):
        """
        Optimizes a molecule.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        None : :class:`NoneType`

        """

        mol._file = f'{uuid4().int}.mol'

        # First write a .mol file of the molecule.
        mol.write(mol._file, conformer=conformer)
        # MOPAC requires a ``.mop`` file as input. This creates a
        # ``.mop`` file holding the molecule.
        self.generate_mop(mol, conformer)
        # Run the optimization.
        self.run_mopac(mol)
        # Update the mol with the generated .pdb file.
        ofile = mol._file.replace('.mol', '.pdb')
        mol.update_from_pdb(ofile, conformer=conformer)

    def run_mopac(self, mol, mopac_path, settings, timeout=7200):
        name, ext = os.path.splitext(mol._file)
        mop_file = f'{name}.mop'

        logger.info(f'Running MOPAC on "{mol.name}".')

        # To run MOPAC a command is issued to the console via
        # ``subprocess.Popen``. The command is the full path of the
        # ``mopac`` program.
        file_root, ext = os.path.splitext(mop_file)
        opt_cmd = [mopac_path, file_root]
        opt_proc = psutil.Popen(opt_cmd,
                                stdout=sp.PIPE,
                                stderr=sp.STDOUT,
                                universal_newlines=True)

        try:
            proc_out, _ = opt_proc.communicate(timeout=timeout)
        except sp.TimeoutExpired:
            logger.warning(
                f'Minimization took too long and was terminated '
                f'by force - {mol.name}'
            )
            self.kill_mopac(mol)

        return

    def kill_mopac(self, mol):
        """
        Kills an in-progress MOPAC run.

        To kill a MOPAC run, a file with the molecule's name and a
        ``.end`` extension is written.

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

    def mop_header(self):
        """
        Returns the header of the ``.mop`` file.

        Returns
        -------
        :class:`str`
            String containing all the MOPAC keywords correctly
            formatted for the input file.

        """

        if self.method == 'OPT':
            grad = f'GNORM={self.minimum_gradient}'
        else:
            grad = ''

        return (
            f'{self.hamiltonian} {self.method} {self.method} {grad} '
            f'EPS={self.eps} CHARGE={self.charge} {self.fileout} LET '
        )

    def generate_mop(self, mol, conformer):
        """
        Creates the ``.mop`` file for the optimization.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule which is to be optimized. Its molecular
            structure file is converted to a ``.mop`` file.

        conformer : :class:`int`
            The conformer of the molecule to be used.

        Returns
        -------
        :class:`str`
            The full path of the newly created ``.mop`` file.

        """

        logger.info(f'Creating ".mop" file - {mol.name}.')
        name, ext = os.path.splitext(mol._file)
        mop_file = f'{name}.mop'

        # Generate the mop file containing the MOPAC run info.
        with open(mop_file, 'w') as mop:
            # Line for the run info.
            mop.write(f'{self.mop_header()}\n')
            # Line with the name of the molecule.
            mop.write(f'{name}\n\n')

            # Write the structural info.
            atom_block = ''
            for i in range(mol.mol.GetNumAtoms()):
                elem = mol.atom_symbol(i)
                x, y, z = mol.atom_coords(i, conformer)
                atom_block += (
                    f'{elem}   {x}   +1  {y}   +1  {z}   +1 \n'
                )

            mop.write(atom_block)

        return mop_file
