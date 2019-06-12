"""
Defines energy calculators.

See :mod:`.energy`.

"""

import rdkit.Chem.AllChem as rdkit
import logging
from functools import wraps
import subprocess as sp
import uuid
import os
import shutil
from ...utilities import valid_xtb_solvent, XTBExtrators


logger = logging.getLogger(__name__)


class EnergyError(Exception):
    """
    Indicates a failed energy calculation.

    """

    ...


def _add_cache_use(energy):
    """
    Makes :meth:`~EnergyCalculator.energy` use the cache.

    Decorates `energy` so that before running it checks if the
    :class:`.Molecule` and conformer have already had their energy
    calculated. If so, and :attr:`~EnergyCalculator.use_cache` is
    ``True``, then the energy value in the
    :attr:`~EnergyCalculator.cache` is returned.

    Parameters
    ----------
    energy : :class:`function`
        A function which is to have cache use added to it.

    Returns
    -------
    :class:`function`
        The decorated function.

    """

    @wraps(energy)
    def inner(self, mol, conformer=-1):
        key = (mol.key, conformer)
        if self.use_cache and key in self.cache:
            logger.info(
                f'Using cached energy value with '
                f'"{mol.name}" conformer {conformer}.'
            )
            return self.cache[key]
        else:
            e = energy(self, mol, conformer)
            if self.use_cache:
                self.cache[key] = e
            return e

    return inner


class EnergyCalculator:
    """
    Calculate the energy of molecules.

    Attributes
    ----------
    cache : :class:`dict`
        A :class:`dict` hodling calculated energy values of the form

        .. code-block:: python

            cache = {
                (mol1, conformer1): 12.2,
                (mol1, conformer2): 2.0,
                (mol2, conformer2): 321.12
            }

        which holds every :class:`Molecule` and conformer whose energy
        was calculated by the :class:`EnergyCalculator`. Here ``mol1``
        and ``mol2`` are :class:`.Molecule` objects and ``conformer1``
        and ``conformer2`` are of type :class:`int`, and are the ids
        of the conformers passed to :meth:`energy`.

    use_cache : :class:`bool`
        If ``True`` :meth:`energy` will not run twice on the same
        molecule and conformer, but will instead return the previously
        calculated value.

    """

    def __init__(self, use_cache=False):
        """
        Initializes a :class:`EnergyCalculator` instance.

        Parameters
        ----------
        use_cache : :class:`bool`, optional
            If ``True`` :meth:`energy` will not run twice on the same
            molecule and conformer, but will instead return the
            previously calculated value.

        """

        self.cache = {}
        self.use_cache = use_cache

    def __init_subclass__(cls, **kwargs):
        cls.energy = _add_cache_use(cls.energy)
        return super().__init_subclass__(**kwargs)

    def energy(self, mol, conformer=-1):
        """
        Calculates the energy of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        conformer : :class:`int`, optinal
            The conformer of `mol` to use.

        Returns
        -------
        :class:`float`
            The energy.

        """

        raise NotImplementedError()


class FormationEnergy(EnergyCalculator):
    """
    Calculates the formation energy of a molecule.

    Attributes
    ----------
    energy_calculator : :class:`EnergyCalculator`
        The :class:`EnergyCalculator` used to calculate the energy of
        all reactant and product molecules.

    reactants : :class:`list` of :class:`.Molecule`
        The reactants. If there are multiples of the same reactant
        then it must appear multiple times in this :class:`list`.

    products : :class:`list` of :class:`.Molecule`
        The molecules which are produced as a result of the formation
        reaction. This :class:`list` must omit the :class:`.Molecule`
        passed to :meth:`energy`. If there are multiples of the same
        product, it must appear multiple times in this :class:`list`.

    reactant_conformers : :class:`list` of :class:`int`
        The conformer ids of the molecules in :attr:`reactants` to use.

    product_conformers : :class:`list` of :class:`int`
        The conformer ids of the molecules in :attr:`products` to use.

    Examples
    --------
    .. code-block:: python

        # Make the molecules needed to calculate formation energy.
        water = StructUnit.smiles_init('[H]O[H]')
        bb1 = StructUnit2(...)
        bb2 = StructUnit3(...)
        cage = Cage([bb1, bb2], FourPlusSix())

        # Make the energy calculator used to calculate energies.
        uff_energy = UFFEnergy()

        # Make the formation energy calculator.
        formation = FormationEnergy(
            energy_calculator=uff_energy,
            reactants=[bb1]*6 + [bb2]*4,
            products=[water]*cage.topology.bonds_made
        )

        # Calculate the formation energy.
        formation_energy = formation.energy(cage)

    """

    def __init__(self,
                 energy_calculator,
                 reactants,
                 products,
                 reactant_conformers=None,
                 product_conformers=None,
                 use_cache=False):
        """
        Initializes a :class:`FormationEnergy` instance.

        Parameters
        ----------
        energy_calculator : :class:`EnergyCalculator`
            The :class:`EnergyCalculator` used to calculate the energy
            of all reactant and product molecules.

        reactants : :class:`list` of :class:`.Molecule`
            The reactants. If there are multiples of the same reactant
            then it must appear multiple times in this :class:`list`.

        products : :class:`list` of :class:`.Molecule`
            The molecules which are produced as a result of the
            formation reaction. This :class:`list` must omit the
            :class:`.Molecule` passed to :meth:`energy`. If there are
            multiples of the same product, it must appear multiple
            times in this :class:`list`.

        reactant_conformers : :class:`list` of :class:`int`, optional
            The conformer ids of the molecules in :attr:`reactants` to
            use.

        product_conformers : :class:`list` of :class:`int`, optional
            The conformer ids of the molecules in :attr:`products` to
            use.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`energy` will not run twice on the same
            molecule and conformer, but will instead return the
            previously calculated value.

        """

        if reactant_conformers is None:
            reactant_conformers = [-1 for i in range(len(reactants))]

        if product_conformers is None:
            product_conformers = [-1 for i in range(len(products))]

        self.energy_calculator = energy_calculator
        self.reactants = reactants
        self.products = products
        self.reactant_conformers = reactant_conformers
        self.product_conformers = product_conformers
        super().__init__(use_cache=use_cache)

    def energy(self, mol, conformer=-1):
        """
        Calculates the formation energy of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose formation energy is to be
            calculated.

        conformer : :class:`int`, optinal
            The conformer of `mol` to use.

        Returns
        -------
        :class:`float`
            The formation energy.

        """

        products = zip(self.products, self.product_conformers)
        product_energy = sum(
            self.energy_calculator.energy(product, conformer)
            for product, conformer in products
        )
        product_energy += self.energy_calculator.energy(mol, conformer)

        reactants = zip(self.reactants, self.reactant_conformers)
        reactant_energy = sum(
            self.energy_calculator.energy(reactant, conformer)
            for reactant, conformer in reactants
        )

        return product_energy - reactant_energy


class MMFFEnergy(EnergyCalculator):
    """
    Uses the MMFF force field to calculate energies.

    Examples
    --------
    .. code-block:: python

        # Create a molecules whose energy we want to know.
        mol1 = StructUnit.smiles_init('CCCNCCCN')
        mol2 = Polymer(...)
        mol3 = Cage(...)

        # Create the energy calculator.
        mmff = MMFFEnergy()

        # Calculate the energies.
        energy1 = mmff.energy(mol1)
        energy2 = mmff.energy(mol2)
        energy3 = mmff.energy(mol3)

    """

    def energy(self, mol, conformer=-1):
        """
        Calculates the energy of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        conformer : :class:`int`, optinal
            The conformer of `mol` to use.

        Returns
        -------
        :class:`float`
            The energy.

        """

        rdkit.GetSSSR(mol.mol)
        mol.mol.UpdatePropertyCache()
        ff = rdkit.MMFFGetMoleculeForceField(
            mol.mol,
            rdkit.MMFFGetMoleculeProperties(mol.mol),
            confId=conformer
        )

        return ff.CalcEnergy()


class UFFEnergy(EnergyCalculator):
    """
    Uses the UFF force field to calculate energies.

    Examples
    --------
    .. code-block:: python

        # Create a molecules whose energy we want to know.
        mol1 = StructUnit.smiles_init('CCCNCCCN')
        mol2 = Polymer(...)
        mol3 = Cage(...)

        # Create the energy calculator.
        uff = UFFEnergy()

        # Calculate the energies.
        energy1 = uff.energy(mol1)
        energy2 = uff.energy(mol2)
        energy3 = uff.energy(mol3)

    """

    def energy(self, mol, conformer=-1):
        """
        Calculates the energy of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        conformer : :class:`int`, optinal
            The conformer of `mol` to use.

        Returns
        -------
        :class:`float`
            The energy.

        """

        mol.mol.UpdatePropertyCache()
        # RingInfo needs to be initialized, else rdkit may raise an
        # error.
        rdkit.GetSSSR(mol.mol)
        ff = rdkit.UFFGetMoleculeForceField(mol.mol, confId=conformer)
        return ff.CalcEnergy()


class XTBEnergy(EnergyCalculator):
    """
    Uses GFN-xTB to calculate energies and other properties.

    By default, :meth:`energy` will extract other properties of the
    :class:`.Molecule` and conformer passed to :meth:`energy`, which
    will be saved as attributes of :class:`.XTBEnergy`.

    Notes
    -----
    When running :meth:`energy`, this calculator changes the
    present working directory with :func:`os.chdir`. The original
    working directory will be restored even if an error is raised so
    unless multi-threading is being used this implementation detail
    should not matter.

    If multi-threading is being used an error could occur if two
    different threads need to know about the current working directory
    as this :class:`.XTBEnergy` can change it from under them.

    Note that this does not have any impact on multi-processing,
    which should always be safe.

    Attributes
    ----------
    xtb_path : :class:`str`
        The path to the xTB executable.

    gfn_version : :class:`str`
        Parameterization of GFN to use in xTB.
        For details see
        https://xtb-docs.readthedocs.io/en/latest/basics.html

    output_dir : :class:`str`
        The name of the directory into which files generated during
        the optimization are written, if ``None`` then
        :func:`uuid.uuid4` is used.

    num_cores : :class:`int`
        The number of cores xTB should use.

    electronic_temperature : :class:`int`
        Electronic temperature to use (in K).

    solvent : :class:`str`
        Solvent to use in GBSA implicit solvation method.
        For details see
        https://xtb-docs.readthedocs.io/en/latest/gbsa.html

    solvent_grid : :class:`str`
        Grid level to use in SASA calculations for GBSA implicit
        solvent.
        Can be one of ``'normal'``, ``'tight'``, ``'verytight'``
        or ``'extreme'``
        For details see
        https://xtb-docs.readthedocs.io/en/latest/gbsa.html

    charge : :class:`int`
        Formal molecular charge.

    unpaired_electrons : :class:`int`
        Number of unpaired electrons.

    use_cache : :class:`bool`
        If ``True`` :meth:`energy` will not run twice on the same
        molecule and conformer.

    unlimited_memory : :class: `bool`
        If ``True`` :meth:`energy` will be run without constraints on
        the stacksize. If memory issues are encountered, this should be
        ``True``, however this may raise issues on clusters.

    total_energies : :class:`dict`
        :class:`dict` of the total energy of :class:`.Molecule` and
        conformer passed to :meth:`energy`.
        Key is ``(mol, conformer)``.

    homo_lumo_gaps : :class:`dict`
        :class:`dict` of the HOMO-LUMO gap of :class:`.Molecule` and
        conformer passed to :meth:`energy`.
        Key is ``(mol, conformer)``.

    fermi_levels : :class:`dict`
        :class:`dict` of the Fermi level of :class:`.Molecule` and
        conformer passed to :meth:`energy`.
        Key is ``(mol, conformer)``.

    homo_lumo_orbitals : :class:`dict`
        :class:`dict` of the HOMO-LUMO orbital properties of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        Key is ``(mol, conformer)``.

    qonly_dipole_moms : :class:`dict`
        :class:`dict` of the `q only` dipole moment of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        Key is ``(mol, conformer)``.

    full_dipole_moms : :class:`dict`
        :class:`dict` of the `full` dipole moment of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        Key is ``(mol, conformer)``.

    qonly_quadrupole_moms : :class:`dict`
        :class:`dict` of the `q only` quadrupole moment of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        Key is ``(mol, conformer)``.

    qdip_quadrupole_moms : :class:`dict`
        :class:`dict` of the `q+dip` quadrupole moment of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        Key is ``(mol, conformer)``.

    full_quadrupole_moms : :class:`dict`
        :class:`dict` of the `full` quadrupole moment of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        Key is ``(mol, conformer)``.

    Examples
    --------
    .. code-block:: python

        mol = StructUnit.smiles_init('NCCNCCN', ['amine'])
        gfnxtb = XTBEnergy('/opt/gfnxtb/xtb')
        gfnxtb.energy(mol)

    Note that for :class:`.MacroMolecule` objects assembled by ``stk``
    :class:`XTBEnergy` should usually be used after optimization with
    some other method. This is because xTB only uses xyz coordinates
    as input and so will not recognize the long bonds created during
    assembly. An optimizer which can minimize these bonds should be
    used before :class:`XTBEnergy`.

    .. code-block:: python

        bb1 = StructUnit2.smiles_init('NCCNCCN', ['amine'])
        bb2 = StructUnit2.smiles_init('O=CCCC=O', ['aldehyde'])
        polymer = Polymer([bb1, bb2], Linear("AB", [0, 0], 3))

        uff = UFF()
        uff.optimize(polymer)
        xtb = XTBEnergy(xtb_path='/opt/gfnxtb/xtb',
                            unlimited_memory=True)
        xtb.energy(polymer)

    Energies and other properties of optimized structures can be
    extracted using :class:`XTBEnergy`. However, be aware that
    :class:`XTBEnergy` will not check the quality of an optimized
    structure.

    .. code-block:: python

        # runs calculation and returns energy
        p_totalenergy = xtb.energy(polymer, conformer)

        # extracts properties from energy calculator for given
        #  conformer
        p_homo_lumo_gap = xtb.homo_lumo_gaps[(polymer, conformer)]
        p_fermi_levels = xtb.fermi_levels[(polymer, conformer)]
        p_homo_lumo_orbitals = xtb.homo_lumo_orbitals[
            (polymer, conformer)
        ]
        p_full_dipole_moms = xtb.full_dipole_moms[
            (polymer, conformer)
        ]
        p_full_quadrupole_moms = xtb.full_quadrupole_moms[
            (polymer, conformer)
        ]

        # the total energy can be extracted at any point from the
        # calculator
        polymer_totalenergy = xtb.total_energies[(polymer, conformer)]

    See Also
    --------
    #. https://xtb-docs.readthedocs.io/en/latest/setup.html

    """
    def __init__(self,
                 xtb_path,
                 gfn_version='2',
                 output_dir=None,
                 num_cores=1,
                 electronic_temperature=300,
                 solvent=None,
                 solvent_grid='normal',
                 charge=0,
                 unpaired_electrons=0,
                 use_cache=False,
                 unlimited_memory=False):
        """
        Initializes a :class:`XTBEnergy` instance.

        Parameters
        ----------
        xtb_path : :class:`str`
            The path to the xTB executable.

        gfn_version : :class:`str`, optional
            Parameterization of GFN to use in xTB.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/basics.html

        output_dir : :class:`str`, optional
            The name of the directory into which files generated during
            the optimization are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        num_cores : :class:`int`, optional
            The number of cores xTB should use.

        electronic_temperature : :class:`int`, optional
            Electronic temperature to use (in K).

        solvent : :class:`str`, optional
            Solvent to use in GBSA implicit solvation method.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html

        solvent_grid : :class:`str`, optional
            Grid level to use in SASA calculations for GBSA implicit
            solvent.
            Can be one of ``'normal'``, ``'tight'``, ``'verytight'``
            or ``'extreme'``
            For details see
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html

        charge : :class:`int`, optional
            Formal molecular charge.

        unpaired_electrons : :class:`int`, optional
            Number of unpaired electrons.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`energy` will not run twice on the same
            molecule and conformer.

        unlimited_memory : :class: `bool`, optional
            If ``True`` :meth:`energy` will be run without constraints
            on the stacksize. If memory issues are encountered, this
            should be ``True``, however this may raise issues on
            clusters.

        """
        self.xtb_path = xtb_path
        self.gfn_version = gfn_version
        self.output_dir = output_dir
        self.num_cores = str(num_cores)
        self.electronic_temperature = str(electronic_temperature)
        self.solvent = solvent
        if self.solvent is not None:
            self.solvent = solvent.lower()
            valid_xtb_solvent(
                gfn_version=self.gfn_version,
                solvent=self.solvent
            )
        self.solvent_grid = solvent_grid
        self.charge = str(charge)
        self.unpaired_electrons = str(unpaired_electrons)
        self.unlimited_memory = unlimited_memory

        self.total_energies = {}
        self.homo_lumo_gaps = {}
        self.fermi_levels = {}
        self.homo_lumo_orbitals = {}
        self.qonly_dipole_moms = {}
        self.full_dipole_moms = {}
        self.qonly_quadrupole_moms = {}
        self.qdip_quadrupole_moms = {}
        self.full_quadrupole_moms = {}
        super().__init__(use_cache=use_cache)

    def _get_properties(self, mol, conformer, output_file):
        """
        Extracts properties from GFN-xTB energy calculation.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        conformer : :class:`int`, optinal
            The conformer of `mol` to use.

        output_file : :class: `str`
            Name of output file with xTB results.

        Returns
        -------
        None : :class:`NoneType`

        """
        xtbext = XTBExtrators(output_file=output_file)

        # Get properties from output file.
        key = (mol, conformer)
        self.total_energies[key] = xtbext.total_energy()
        self.homo_lumo_gaps[key] = xtbext.homo_lumo_gap()
        self.fermi_levels[key] = xtbext.fermi_level()
        self.homo_lumo_orbitals[key] = xtbext.homo_lumo_occ()
        self.qonly_dipole_moms[key] = xtbext.qonly_dipole_mom()
        self.full_dipole_moms[key] = xtbext.full_dipole_mom()
        self.qonly_quadrupole_moms[key] = xtbext.qonly_quadrupole_mom()
        self.qdip_quadrupole_moms[key] = xtbext.qdip_quadrupole_mom()
        self.full_quadrupole_moms[key] = xtbext.full_quadrupole_mom()

    def _write_and_run_command(self, mol, conformer):
        """
        Writes and runs the command for xTB.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        conformer : :class:`int`, optinal
            The conformer of `mol` to use.

        Returns
        -------
        :class:`str`
            Name of output file with xTB results.

        """
        xyz = 'input_structure.xyz'
        out_file = 'energy.output'
        mol.write(xyz, conformer=conformer)
        # Modify the memory limit.
        if self.unlimited_memory:
            cmd = ['ulimit -s unlimited ;']
            # Allow multiple shell commands to be run in one
            # subprocess.
            shell = True
        else:
            cmd = []
            shell = False
        cmd.append(self.xtb_path)
        cmd.append(xyz)
        # Set the GFN Parameterization.
        cmd.append('--gfn')
        cmd.append(self.gfn_version)
        # Set the number of cores.
        cmd.append('--parallel')
        cmd.append(self.num_cores)
        # Add eletronic temp term to cmd.
        cmd.append('--etemp')
        cmd.append(self.electronic_temperature)
        # Write the solvent section of cmd.
        if self.solvent is not None:
            cmd.append('--gbsa')
            cmd.append(self.solvent)
            if self.solvent_grid != 'normal':
                cmd.append(self.solvent_grid)
        # Write the charge section of cmd.
        cmd.append('--chrg')
        cmd.append(self.charge)
        # Write the unpaired_electrons section of cmd.
        cmd.append('--uhf')
        cmd.append(self.unpaired_electrons)

        cmd = ' '.join(cmd)
        f = open(out_file, 'w')
        # Uses the shell if unlimited_memory = True and waits until the
        # subproces is complete. This is required to be able to run the
        # unlimited_memory_cmd and GFN calculation in one command,
        # which is then closed, which minimizes the risk of
        # unrestricting the memory limits.
        sp.call(
            cmd,
            stdin=sp.PIPE,
            stdout=f,
            stderr=sp.PIPE,
            shell=shell
        )
        f.close()
        return out_file

    def energy(self, mol, conformer=-1):
        """
        Calculates the energy `mol` using xTB.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        conformer : :class:`int`, optinal
            The conformer of `mol` to use.

        Returns
        -------
        :class:`float`
            The energy.

        """

        if conformer == -1:
            conformer = mol.mol.GetConformer(conformer).GetId()

        if self.output_dir is None:
            output_dir = str(uuid.uuid4().int)
        else:
            output_dir = self.output_dir
        output_dir = os.path.abspath(output_dir)

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        os.mkdir(output_dir)
        init_dir = os.getcwd()
        try:
            os.chdir(output_dir)
            out_file = self._write_and_run_command(
                mol=mol,
                conformer=conformer
            )
            self._get_properties(
                mol=mol,
                conformer=conformer,
                output_file=out_file
            )
        finally:
            os.chdir(init_dir)
        return self.total_energies[(mol, conformer)]


class XTBFreeEnergy(EnergyCalculator):
    """
    Uses GFN-xTB to calculate free energies and other properties.

    By default, :meth:`energy` will extract other properties of the
    :class:`.Molecule` and conformer passed to :meth:`energy`, which
    will be saved as attributes of :class:`.XTBEnergy`.

    Notes
    -----
    When running :meth:`energy`, this calculator changes the
    present working directory with :func:`os.chdir`. The original
    working directory will be restored even if an error is raised so
    unless multi-threading is being used this implementation detail
    should not matter.

    If multi-threading is being used an error could occur if two
    different threads need to know about the current working directory
    as this :class:`.XTBFreeEnergy` can change it from under them.

    Note that this does not have any impact on multi-processing,
    which should always be safe.

    Attributes
    ----------
    xtb_path : :class:`str`
        The path to the xTB executable.

    gfn_version : :class:`str`
        Parameterization of GFN to use in xTB.
        For details see
        https://xtb-docs.readthedocs.io/en/latest/basics.html

    output_dir : :class:`str`
        The name of the directory into which files generated during
        the optimization are written, if ``None`` then
        :func:`uuid.uuid4` is used.

    num_cores : :class:`int`
        The number of cores xTB should use.

    electronic_temperature : :class:`int`
        Electronic temperature to use (in K).

    solvent : :class:`str`
        Solvent to use in GBSA implicit solvation method.
        For details see
        https://xtb-docs.readthedocs.io/en/latest/gbsa.html

    solvent_grid : :class:`str`
        Grid level to use in SASA calculations for GBSA implicit
        solvent.
        Can be one of ``'normal'``, ``'tight'``, ``'verytight'``
        or ``'extreme'``
        For details see
        https://xtb-docs.readthedocs.io/en/latest/gbsa.html

    charge : :class:`int`
        Formal molecular charge.

    unpaired_electrons : :class:`int`
        Number of unpaired electrons.

    use_cache : :class:`bool`
        If ``True`` :meth:`energy` will not run twice on the same
        molecule and conformer.

    unlimited_memory : :class: `bool`
        If ``True`` :meth:`energy` will be run without constraints on
        the stacksize. If memory issues are encountered, this should be
        ``True``, however this may raise issues on clusters.

    total_energies : :class:`dict`
        :class:`dict` of the total energy of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        Key is ``(mol, conformer)``.

    homo_lumo_gaps : :class:`dict`
        :class:`dict` of the HOMO-LUMO gap of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        Key is ``(mol, conformer)``.

    fermi_levels : :class:`dict`
        :class:`dict` of the Fermi level of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        Key is ``(mol, conformer)``.

    homo_lumo_orbitals : :class:`dict`
        :class:`dict` of the HOMO-LUMO orbital properties of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        Key is ``(mol, conformer)``.

    qonly_dipole_moms : :class:`dict`
        :class:`dict` of the `q only` dipole moment of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        Key is ``(mol, conformer)``.

    full_dipole_moms : :class:`dict`
        :class:`dict` of the `full` dipole moment of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        Key is ``(mol, conformer)``.

    qonly_quadrupole_moms : :class:`dict`
        :class:`dict` of the `q only` quadrupole moment of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        Key is ``(mol, conformer)``.

    qdip_quadrupole_moms : :class:`dict`
        :class:`dict` of the `q+dip` quadrupole moment of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        Key is ``(mol, conformer)``.

    full_quadrupole_moms : :class:`dict`
        :class:`dict` of the `full` quadrupole moment of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        Key is ``(mol, conformer)``.

    total_free_energies : :class:`dict`
        :class:`dict` of the total free energy of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        Key is ``(mol, conformer)``.

    frequencies : :class:`dict`
        :class:`dict` of the vibrational frequencies of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        Key is ``(mol, conformer)``.


    Examples
    --------
    .. code-block:: python

        mol = StructUnit.smiles_init('NCCNCCN', ['amine'])
        xtb = XTBFreeEnergy('/opt/gfnxtb/xtb')
        xtb.energy(mol)

    Note that for :class:`.MacroMolecule` objects assembled by ``stk``
    :class:`XTBFreeEnergy` should usually be used after optimization
    with some other method. This is because xTB only uses xyz
    coordinates as input and so will not recognize the long bonds
    created during assembly. An optimizer which can minimize these
    bonds should be used before :class:`XTBFreeEnergy`.

    .. code-block:: python

        bb1 = StructUnit2.smiles_init('NCCNCCN', ['amine'])
        bb2 = StructUnit2.smiles_init('O=CCCC=O', ['aldehyde'])
        polymer = Polymer([bb1, bb2], Linear("AB", [0, 0], 3))

        uff = UFF()
        uff.optimize(polymer)
        xtb = XTBFreeEnergy(
            xtb_path='/opt/gfnxtb/xtb',
            unlimited_memory=True
        )
        xtb.energy(polymer)

    Energies and other properties of optimized structures can be
    extracted using :class:`XTBFreeEnergy`. For example vibrational
    frequencies, the HOMO-LUMO gap and thermodynamic properties
    (such as the total free energy) can be calculated after an
    optimization with very tight contraints with an implicit solvent
    (THF). Very tight criteria are required to ensure that no negative
    vibrational frequencies are present, however :class:`XTBFreeEnergy`
    will not check for the presence of negative frequencies or poorly
    optimized structures. It is recommended that :class:`XTBFreeEnergy`
    is only used on well optimized structures.

    .. code-block:: python

        xtb = OptimizerSequence(
            stk.UFF(),
            stk.XTB(
                xtb_path='/opt/gfnxtb/xtb',
                unlimited_memory=True,
                opt_level='verytight',
                solvent='THF'
            )
        )
        xtb.optimize(polymer)

        xtb_fe = XTBFreeEnergy(
            xtb_path='/opt/gfnxtb/xtb',
            unlimited_memory=True,
            solvent='THF'
        )
        # runs calculation and returns energy
        poly_totalenergy = xtb_fe.energy(polymer, conformer)

        # extracts properties from energy calculator for given
        # conformer
        poly_free_energy = xtb_fe.total_free_energies[
            (polymer, conformer)
        ]
        poly_freq = xtb_fe.frequencies[(polymer, conformer)]
        poly_homo_lumo_gap = xtb_fe.homo_lumo_gaps[
            (polymer, conformer)
        ]

    See Also
    --------
    #. https://xtb-docs.readthedocs.io/en/latest/setup.html

    """

    def __init__(self,
                 xtb_path,
                 gfn_version='2',
                 output_dir=None,
                 num_cores=1,
                 electronic_temperature=300,
                 solvent=None,
                 solvent_grid='normal',
                 charge=0,
                 unpaired_electrons=0,
                 use_cache=False,
                 unlimited_memory=False):
        """
        Initializes a :class:`XTBFreeEnergy` instance.

        Parameters
        ----------
        xtb_path : :class:`str`
            The path to the xTB executable.

        gfn_version : :class:`str`, optional
            Parameterization of GFN to use in xTB.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/basics.html

        output_dir : :class:`str`, optional
            The name of the directory into which files generated during
            the optimization are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        num_cores : :class:`int`, optional
            The number of cores xTB should use.

        electronic_temperature : :class:`int`, optional
            Electronic temperature to use (in K).

        solvent : :class:`str`, optional
            Solvent to use in GBSA implicit solvation method.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html

        solvent_grid : :class:`str`, optional
            Grid level to use in SASA calculations for GBSA implicit
            solvent.
            Can be one of ``'normal'``, ``'tight'``, ``'verytight'``
            or ``'extreme'``
            For details see
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html

        charge : :class:`int`, optional
            Formal molecular charge.

        unpaired_electrons : :class:`int`, optional
            Number of unpaired electrons.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`energy` will not run twice on the same
            molecule and conformer.

        unlimited_memory : :class: `bool`, optional
            If ``True`` :meth:`energy` will be run without constraints
            on the stacksize. If memory issues are encountered, this
            should be ``True``, however this may raise issues on
            clusters.

        """
        self.xtb_path = xtb_path
        self.gfn_version = gfn_version
        self.output_dir = output_dir
        self.num_cores = str(num_cores)
        self.electronic_temperature = str(electronic_temperature)
        self.solvent = solvent
        if self.solvent is not None:
            self.solvent = solvent.lower()
            valid_xtb_solvent(
                gfn_version=self.gfn_version,
                solvent=self.solvent
            )
        self.solvent_grid = solvent_grid
        self.charge = str(charge)
        self.unpaired_electrons = str(unpaired_electrons)
        self.unlimited_memory = unlimited_memory

        self.total_energies = {}
        self.homo_lumo_gaps = {}
        self.fermi_levels = {}
        self.homo_lumo_orbitals = {}
        self.qonly_dipole_moms = {}
        self.full_dipole_moms = {}
        self.qonly_quadrupole_moms = {}
        self.qdip_quadrupole_moms = {}
        self.full_quadrupole_moms = {}
        self.total_free_energies = {}
        self.frequencies = {}
        super().__init__(use_cache=use_cache)

    def _get_properties(self, mol, conformer, output_file):
        """
        Extracts properties from GFN-xTB energy calculation.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        conformer : :class:`int`, optinal
            The conformer of `mol` to use.

        output_file : :class: `str`
            Name of output file with xTB results.

        Returns
        -------
        None : :class:`NoneType`

        """
        xtbext = XTBExtrators(output_file=output_file)

        # get properties from output string
        key = (mol, conformer)
        self.total_energies[key] = xtbext.total_energy()
        self.homo_lumo_gaps[key] = xtbext.homo_lumo_gap()
        self.fermi_levels[key] = xtbext.fermi_level()
        self.homo_lumo_orbitals[key] = xtbext.homo_lumo_occ()
        self.qonly_dipole_moms[key] = xtbext.qonly_dipole_mom()
        self.full_dipole_moms[key] = xtbext.full_dipole_mom()
        self.qonly_quadrupole_moms[key] = xtbext.qonly_quadrupole_mom()
        self.qdip_quadrupole_moms[key] = xtbext.qdip_quadrupole_mom()
        self.full_quadrupole_moms[key] = xtbext.full_quadrupole_mom()
        self.total_free_energies[key] = xtbext.total_free_energy()
        self.frequencies[key] = xtbext.frequencies()

    def _write_and_run_command(self, mol, conformer):
        """
        Writes and runs the command for xTB.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        conformer : :class:`int`, optinal
            The conformer of `mol` to use.

        Returns
        -------
        :class:`str`
            Name of output file with xTB results.

        """
        xyz = 'input_structure.xyz'
        out_file = 'free_energy.output'
        mol.write(xyz, conformer=conformer)
        # Modify the memory limit.
        if self.unlimited_memory:
            cmd = ['ulimit -s unlimited ;']
            # Allow multiple shell commands to be run in one
            # subprocess.
            shell = True
        else:
            cmd = []
            shell = False
        cmd.append(self.xtb_path)
        cmd.append(xyz)
        # Set the GFN Parameterization.
        cmd.append('--gfn')
        cmd.append(self.gfn_version)
        # Turn on the hessian calculation.
        cmd.append('--hess')
        # Set the number of cores.
        cmd.append('--parallel')
        cmd.append(self.num_cores)
        # Add eletronic temp term to cmd.
        cmd.append('--etemp')
        cmd.append(self.electronic_temperature)
        # Write the solvent section of cmd.
        if self.solvent is not None:
            cmd.append('--gbsa')
            cmd.append(self.solvent)
            if self.solvent_grid != 'normal':
                cmd.append(self.solvent_grid)
        # Write the charge section of cmd.
        cmd.append('--chrg')
        cmd.append(self.charge)
        # Write the unpaired_electrons section of cmd.
        cmd.append('--uhf')
        cmd.append(self.unpaired_electrons)

        cmd = ' '.join(cmd)
        f = open(out_file, 'w')
        # Uses the shell if unlimited_memory = True and waits until the
        # subproces is complete. This is required to be able to run the
        # unlimited_memory_cmd and GFN calculation in one command,
        # which is then closed, which minimizes the risk of
        #  unrestricting the memory limits.
        sp.call(
            cmd,
            stdin=sp.PIPE,
            stdout=f,
            stderr=sp.PIPE,
            shell=shell
        )
        f.close()
        return out_file

    def energy(self, mol, conformer=-1):
        """
        Calculates the energy `mol` using xTB.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        conformer : :class:`int`, optinal
            The conformer of `mol` to use.

        Returns
        -------
        :class:`float`
            The energy.

        """

        if conformer == -1:
            conformer = mol.mol.GetConformer(conformer).GetId()

        if self.output_dir is None:
            output_dir = str(uuid.uuid4().int)
        else:
            output_dir = self.output_dir
        output_dir = os.path.abspath(output_dir)

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        os.mkdir(output_dir)
        init_dir = os.getcwd()
        try:
            os.chdir(output_dir)
            out_file = self._write_and_run_command(
                mol=mol,
                conformer=conformer
            )
            self._get_properties(
                mol=mol,
                conformer=conformer,
                output_file=out_file
            )
        finally:
            os.chdir(init_dir)
        return self.total_energies[(mol, conformer)]
