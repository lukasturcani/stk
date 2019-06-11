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
from ...utilities import valid_XTB_solvent, XTBExts


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
            molecule and conformer, but will instead return the previously
            calculated value.

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
    Uses GFN-xTB to calculate energies and other properties of molecules.

    Notes
    -----
    When running :meth:`energy`, this calculator changes the
    present working directory with :func:`os.chdir`. The original
    working directory will be restored even if an error is raised so
    unless multi-threading is being used this implementation detail
    should not matter.

    If multi-threading is being used an error could occur if two
    different threads need to know about the current working directory
    as this :class:`.EnergyCalculator` can change it from under them.

    Note that this does not have any impact on multi-processing,
    which should always be safe.

    Documentation for xTB available:
    https://xtb-docs.readthedocs.io/en/latest/setup.html

    Attributes
    ----------
    xtb_path : :class:`str`
        The path to the xTB executable.

    gfn_version : :class:`str`
        Parameterization of GFN to use in xTB.
        For details:
            https://xtb-docs.readthedocs.io/en/latest/basics.html

    output_dir : :class:`str`, optional
        The name of the directory into which files generated during
        the optimization are written, if ``None`` then
        :func:`uuid.uuid4` is used.

    num_cores : :class:`int`
        The number of cores for xTB to use. Requires appropriate setup
        of xTB by user.

    etemp : :class:`int`, optional
        Electronic temperature to use (in K). Defaults to 300K.

    solvent : :class:`str`, optional
        Solvent to use in GBSA implicit solvation method.
        For details:
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html

    solvent_grid : :class:`str`, optional
        Grid level to use in SASA calculations for GBSA implicit solvent.
        Options:
            normal, tight, verytight, extreme
        For details:
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html

    charge : :class:`str`, optional
        Formal molecular charge. `-` should be used to indicate sign.

    multiplicity : :class:`str`, optional
        Number of unpaired electrons.

    use_cache : :class:`bool`, optional
        If ``True`` :meth:`energy` will not run twice on the same
        molecule and conformer.

    mem_ulimit : :class: `bool`, optional
        If ``True`` :meth:`energy` will be run without constraints on
        the stacksize. If memory issues are encountered, this should be
        ``True``, however this may raise issues on clusters.

    Examples
    --------
    .. code-block:: python

        mol = StructUnit.smiles_init('NCCNCCN', ['amine'])
        gfnxtb = XTBEnergy('/opt/gfnxtb/xtb')
        gfnxtb.energy(mol)

    Note that for :class:`.MacroMolecule` objects assembled by ``stk``
    :class:`XTBEnergy` should usually be used after optimization with some
    other method. This is because xTB only uses xyz coordinates as input
    and so will not recognize the long bonds created during assembly.
    An optimizer which can minimize these bonds should be used before
    :class:`XTBEnergy`.

    .. code-block:: python

        bb1 = StructUnit2.smiles_init('NCCNCCN', ['amine'])
        bb2 = StructUnit2.smiles_init('O=CCCC=O', ['aldehyde'])
        polymer = Polymer([bb1, bb2], Linear("AB", [0, 0], 3))

        uff = UFF()
        uff.optimize(polymer)
        xtb = XTBEnergy(xtb_path='/opt/gfnxtb/xtb',
                            mem_ulimit=True)
        xtb.energy(polymer)

    Energies and other properties of optimized structures can be
    extracted using :class:`XTBEnergy`. However, be aware that
    :class:`XTBEnergy` will not check the quality of an optimized structure.

    .. code-block:: python

        xtb = OptimizerSequence(
            UFF(),
            XTB(xtb_path='/opt/gfnxtb/xtb',
                mem_ulimit=True,
                opt_level='normal',
                solvent='THF')
        )
        xtb.optimize(polymer)

        xtbenergy = XTBEnergy(xtb_path='/opt/gfnxtb/xtb',
                              mem_ulimit=True,
                              solvent='THF')
        # runs calculation and returns energy
        polymer_totalenergy = xtbenergy.energy(polymer, conformer)

        # extracts properties from energy calculator for given conformer
        polymer_homo_lumo_gap = xtbenergy.homo_lumo_gaps[(polymer, conformer)]
        polymer_fermi_levels = xtbenergy.fermi_levels[(polymer, conformer)]
        polymer_homo_lumo_orbitals = xtbenergy.homo_lumo_orbitals[(polymer, conformer)]
        polymer_Qonly_dipole_moments = xtbenergy.Qonly_dipole_moments[(polymer, conformer)]
        polymer_full_dipole_moments = xtbenergy.full_dipole_moments[(polymer, conformer)]
        polymer_Qonly_quadrupole_moments = xtbenergy.Qonly_quadrupole_moments[(polymer, conformer)]
        polymer_QDip_quadrupole_moments = xtbenergy.QDip_quadrupole_moments[(polymer, conformer)]
        polymer_full_quadrupole_moments = xtbenergy.full_quadrupole_moments[(polymer, conformer)]

        # the total energy can be extracted at any point from the calculator
        polymer_totalenergy = xtbenergy.total_energies[(polymer, conformer)]

    """
    def __init__(self,
                 xtb_path,
                 gfn_version='2',
                 output_dir=None,
                 num_cores=1,
                 etemp=300,
                 solvent=None,
                 solvent_grid='normal',
                 charge=None,
                 multiplicity=None,
                 use_cache=False,
                 mem_ulimit=False):
        """
        Initializes a :class:`XTBEnergy` instance.

        Parameters
        ----------
        xtb_path : :class:`str`
            The path to the xTB executable.

        gfn_version : :class:`str`
            Parameterization of GFN to use in xTB.
            For details:
                https://xtb-docs.readthedocs.io/en/latest/basics.html

        output_dir : :class:`str`, optional
            The name of the directory into which files generated during
            the optimization are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        num_cores : :class:`int`
            The number of cores for xTB to use. Requires appropriate setup
            of xTB by user.

        etemp : :class:`int`, optional
            Electronic temperature to use (in K). Defaults to 300K.

        solvent : :class:`str`, optional
            Solvent to use in GBSA implicit solvation method.
            For details:
                https://xtb-docs.readthedocs.io/en/latest/gbsa.html

        solvent_grid : :class:`str`, optional
            Grid level to use in SASA calculations for GBSA implicit solvent.
            Options:
                normal, tight, verytight, extreme
            For details:
                https://xtb-docs.readthedocs.io/en/latest/gbsa.html

        charge : :class:`str`, optional
            Formal molecular charge. `-` should be used to indicate sign.

        multiplicity : :class:`str`, optional
            Number of unpaired electrons.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`energy` will not run twice on the same
            molecule and conformer.

        mem_ulimit : :class: `bool`, optional
            If ``True`` :meth:`energy` will be run without constraints on
            the stacksize. If memory issues are encountered, this should be
            ``True``, however this may raise issues on clusters.

        """
        self.xtb_path = xtb_path
        self.gfn_version = gfn_version
        self.output_dir = output_dir
        self.num_cores = str(num_cores)
        self.etemp = str(etemp)
        self.solvent = solvent
        if self.solvent is not None:
            self.solvent = solvent.lower()
            valid_XTB_solvent(gfn_version=self.gfn_version,
                                 solvent=self.solvent)
        self.solvent_grid = solvent_grid
        self.charge = charge
        self.multiplicity = multiplicity
        self.mem_ulimit = mem_ulimit
        # properties
        self.total_energies = {}
        self.homo_lumo_gaps = {}
        self.fermi_levels = {}
        self.homo_lumo_orbitals = {}
        self.Qonly_dipole_moments = {}
        self.full_dipole_moments = {}
        self.Qonly_quadrupole_moments = {}
        self.QDip_quadrupole_moments = {}
        self.full_quadrupole_moments = {}
        super().__init__(use_cache=use_cache)

    def __get_properties(self, mol, conformer, output_file):
        """
        Extracts desired properties from GFN-xTB single point energy calculation.

        """
        XTBExt = XTBExts(output_file=output_file)

        # get properties from output string
        self.total_energies[(mol, conformer)] = XTBExt.ext_total_energy()
        self.homo_lumo_gaps[(mol, conformer)] = XTBExt.ext_homo_lumo_gap()
        self.fermi_levels[(mol, conformer)] = XTBExt.ext_fermi_level()
        self.homo_lumo_orbitals[(mol, conformer)] = XTBExt.ext_homo_lumo_occ()
        self.Qonly_dipole_moments[(mol, conformer)] = XTBExt.ext_Qonly_dipole_mom()
        self.full_dipole_moments[(mol, conformer)] = XTBExt.ext_full_dipole_mom()
        self.Qonly_quadrupole_moments[(mol, conformer)] = XTBExt.ext_Qonly_quadrupole_mom()
        self.QDip_quadrupole_moments[(mol, conformer)] = XTBExt.ext_QDip_quadrupole_mom()
        self.full_quadrupole_moments[(mol, conformer)] = XTBExt.ext_full_quadrupole_mom()

    def __write_and_run_command(self, mol, conformer):
        """
        Writes and runs the command for GFN.

        Returns
        -------
        out_file : :class:`str`
            Name of output file with GFN-xTB results.
        """
        xyz = 'input_structure.xyz'
        out_file = 'energy.output'
        mol.write(xyz, conformer=conformer)
        # modify memory limit
        if self.mem_ulimit:
            cmd = ['ulimit -s unlimited ;']
            # allow multiple shell commands to be run in one subprocess
            shell = True
        else:
            cmd = []
            shell = False
        cmd.append(self.xtb_path)
        cmd.append(xyz)
        # set GFN Parameterization
        if self.gfn_version != '2':
            cmd.append('--gfn')
            cmd.append(self.gfn_version)
        # set number of cores
        cmd.append('--parallel')
        cmd.append(self.num_cores)
        # add eletronic temp term
        if self.etemp != '300':
            cmd.append('--etemp')
            cmd.append(self.etemp)
        # write solvent section of cmd
        if self.solvent is not None:
            cmd.append('--gbsa')
            cmd.append(self.solvent)
            if self.solvent_grid != 'normal':
                cmd.append(self.solvent_grid)
        # write charge section of cmd
        if self.charge is not None:
            cmd.append('--chrg')
            cmd.append(self.charge)
        # write multiplicity section of cmd
        if self.multiplicity is not None:
            cmd.append('--uhf')
            cmd.append(self.multiplicity)
        cmd = ' '.join(cmd)
        f = open(out_file, 'w')
        # uses the shell if mem_ulimit = True and waits until
        # subproces is complete. This is required to run the mem_ulimit_cmd
        # and GFN calculation in one command, which is then closed, which
        # minimizes the risk of unrestricting the memory limits.
        sp.call(cmd, stdin=sp.PIPE, stdout=f, stderr=sp.PIPE,
                 shell=shell)
        f.close()
        return out_file

    def energy(self, mol, conformer=-1):
        """
        Calculates the energy of molecule `mol` using xTB.

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
            out_file = self.__write_and_run_command(mol=mol,
                                                    conformer=conformer)
            self.__get_properties(mol=mol,
                                  conformer=conformer,
                                  output_file=out_file)
        finally:
            os.chdir(init_dir)
        return self.total_energies[(mol, conformer)]


class XTBFreeEnergy(EnergyCalculator):
    """
    Uses GFN-xTB to calculate free energies, vibrational frequencies and
    other properties of molecules.

    Notes
    -----
    When running :meth:`energy`, this calculator changes the
    present working directory with :func:`os.chdir`. The original
    working directory will be restored even if an error is raised so
    unless multi-threading is being used this implementation detail
    should not matter.

    If multi-threading is being used an error could occur if two
    different threads need to know about the current working directory
    as this :class:`.EnergyCalculator` can change it from under them.

    Note that this does not have any impact on multi-processing,
    which should always be safe.

    Documentation for xTB available:
    https://xtb-docs.readthedocs.io/en/latest/setup.html

    Attributes
    ----------
    xtb_path : :class:`str`
        The path to the xTB executable.

    gfn_version : :class:`str`
        Parameterization of GFN to use in xTB.
        For details:
            https://xtb-docs.readthedocs.io/en/latest/basics.html

    output_dir : :class:`str`, optional
        The name of the directory into which files generated during
        the optimization are written, if ``None`` then
        :func:`uuid.uuid4` is used.

    num_cores : :class:`int`
        The number of cores for xTB to use. Requires appropriate setup
        of xTB by user.

    etemp : :class:`int`, optional
        Electronic temperature to use (in K). Defaults to 300K.

    solvent : :class:`str`, optional
        Solvent to use in GBSA implicit solvation method.
        For details:
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html

    solvent_grid : :class:`str`, optional
        Grid level to use in SASA calculations for GBSA implicit solvent.
        Options:
            normal, tight, verytight, extreme
        For details:
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html

    charge : :class:`str`, optional
        Formal molecular charge. `-` should be used to indicate sign.

    multiplicity : :class:`str`, optional
        Number of unpaired electrons.

    use_cache : :class:`bool`, optional
        If ``True`` :meth:`energy` will not run twice on the same
        molecule and conformer.

    mem_ulimit : :class: `bool`, optional
        If ``True`` :meth:`energy` will be run without constraints on
        the stacksize. If memory issues are encountered, this should be
        ``True``, however this may raise issues on clusters.

    Examples
    --------
    .. code-block:: python

        mol = StructUnit.smiles_init('NCCNCCN', ['amine'])
        xtb = XTBFreeEnergy('/opt/gfnxtb/xtb')
        xtb.energy(mol)

    Note that for :class:`.MacroMolecule` objects assembled by ``stk``
    :class:`XTBFreeEnergy` should usually be used after optimization with some
    other method. This is because xTB only uses xyz coordinates as input
    and so will not recognize the long bonds created during assembly.
    An optimizer which can minimize these bonds should be used before
    :class:`XTBFreeEnergy`.

    .. code-block:: python

        bb1 = StructUnit2.smiles_init('NCCNCCN', ['amine'])
        bb2 = StructUnit2.smiles_init('O=CCCC=O', ['aldehyde'])
        polymer = Polymer([bb1, bb2], Linear("AB", [0, 0], 3))

        uff = UFF()
        uff.optimize(polymer)
        xtb = XTBFreeEnergy(xtb_path='/opt/gfnxtb/xtb',
                            mem_ulimit=True)
        xtb.energy(polymer)

    Energies and other properties of optimized structures can be
    extracted using :class:`XTBFreeEnergy`. For example vibrational frequencies,
    the HOMO-LUMO gap and thermodynamic properties (such as the total free
    energy) can be calculated after an optimization with very tight contraints
    with an implicit solvent (THF). Very tight criteria are required to ensure
    that no negative vibrational frequencies are present, however
    :class:`XTBFreeEnergy` will not check for the presence of negative
    frequencies or poorly optimized structures. It is recommended that
    :class:`XTBFreeEnergy` is only used on well optimized structures.

    .. code-block:: python

        xtb = OptimizerSequence(
            UFF(),
            XTB(xtb_path='/opt/gfnxtb/xtb',
                mem_ulimit=True,
                opt_level='verytight',
                solvent='THF')
        )
        xtb.optimize(polymer)

        xtbfreeenergy = XTBFreeEnergy(xtb_path='/opt/gfnxtb/xtb',
                                      mem_ulimit=True,
                                      solvent='THF')
        # runs calculation and returns energy
        polymer_totalenergy = xtbfreeenergy.energy(polymer, conformer)

        # extracts properties from energy calculator for given conformer
        polymer_free_energy = xtbfreeenergy.total_free_energies[(polymer, conformer)]
        polymer_freq = xtbfreeenergy.frequencies[(polymer, conformer)]
        polymer_homo_lumo_gap = xtbfreeenergy.homo_lumo_gaps[(polymer, conformer)]
    """

    def __init__(self,
                 xtb_path,
                 gfn_version='2',
                 output_dir=None,
                 num_cores=1,
                 etemp=300,
                 solvent=None,
                 solvent_grid='normal',
                 charge=None,
                 multiplicity=None,
                 use_cache=False,
                 mem_ulimit=False):
        """
        Initializes a :class:`XTBFreeEnergy` instance.

        Parameters
        ----------
        xtb_path : :class:`str`
            The path to the xTB executable.

        gfn_version : :class:`str`
            Parameterization of GFN to use in xTB.
            For details:
                https://xtb-docs.readthedocs.io/en/latest/basics.html

        output_dir : :class:`str`, optional
            The name of the directory into which files generated during
            the optimization are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        num_cores : :class:`int`
            The number of cores for xTB to use. Requires appropriate setup
            of xTB by user.

        etemp : :class:`int`, optional
            Electronic temperature to use (in K). Defaults to 300K.

        solvent : :class:`str`, optional
            Solvent to use in GBSA implicit solvation method.
            For details:
                https://xtb-docs.readthedocs.io/en/latest/gbsa.html

        solvent_grid : :class:`str`, optional
            Grid level to use in SASA calculations for GBSA implicit solvent.
            Options:
                normal, tight, verytight, extreme
            For details:
                https://xtb-docs.readthedocs.io/en/latest/gbsa.html

        charge : :class:`str`, optional
            Formal molecular charge. `-` should be used to indicate sign.

        multiplicity : :class:`str`, optional
            Number of unpaired electrons.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`energy` will not run twice on the same
            molecule and conformer.

        mem_ulimit : :class: `bool`, optional
            If ``True`` :meth:`energy` will be run without constraints on
            the stacksize. If memory issues are encountered, this should be
            ``True``, however this may raise issues on clusters.

        """
        self.xtb_path = xtb_path
        self.gfn_version = gfn_version
        self.output_dir = output_dir
        self.num_cores = str(num_cores)
        self.etemp = str(etemp)
        self.solvent = solvent
        if self.solvent is not None:
            self.solvent = solvent.lower()
            valid_XTB_solvent(gfn_version=self.gfn_version,
                                 solvent=self.solvent)
        self.solvent_grid = solvent_grid
        self.charge = charge
        self.multiplicity = multiplicity
        self.mem_ulimit = mem_ulimit
        # properties
        self.total_energies = {}
        self.homo_lumo_gaps = {}
        self.fermi_levels = {}
        self.homo_lumo_orbitals = {}
        self.Qonly_dipole_moments = {}
        self.full_dipole_moments = {}
        self.Qonly_quadrupole_moments = {}
        self.QDip_quadrupole_moments = {}
        self.full_quadrupole_moments = {}
        self.total_free_energies = {}
        self.frequencies = {}
        super().__init__(use_cache=use_cache)

    def __get_properties(self, mol, conformer, output_file):
        """
        Extracts desired properties from xTB single point energy calculation.

        """
        XTBExt = XTBExts(output_file=output_file)

        # get properties from output string
        self.total_energies[(mol, conformer)] = XTBExt.ext_total_energy()
        self.homo_lumo_gaps[(mol, conformer)] = XTBExt.ext_homo_lumo_gap()
        self.fermi_levels[(mol, conformer)] = XTBExt.ext_fermi_level()
        self.homo_lumo_orbitals[(mol, conformer)] = XTBExt.ext_homo_lumo_occ()
        self.Qonly_dipole_moments[(mol, conformer)] = XTBExt.ext_Qonly_dipole_mom()
        self.full_dipole_moments[(mol, conformer)] = XTBExt.ext_full_dipole_mom()
        self.Qonly_quadrupole_moments[(mol, conformer)] = XTBExt.ext_Qonly_quadrupole_mom()
        self.QDip_quadrupole_moments[(mol, conformer)] = XTBExt.ext_QDip_quadrupole_mom()
        self.full_quadrupole_moments[(mol, conformer)] = XTBExt.ext_full_quadrupole_mom()
        self.total_free_energies[(mol, conformer)] = XTBExt.ext_total_free_energy()
        self.frequencies[(mol, conformer)] = XTBExt.ext_frequencies()

    def __write_and_run_command(self, mol, conformer):
        """
        Writes and runs the command for GFN.

        Returns
        -------
        out_file : :class:`str`
            Name of output file with GFN-xTB results.

        """
        xyz = 'input_structure.xyz'
        out_file = 'free_energy.output'
        mol.write(xyz, conformer=conformer)
        # modify memory limit
        if self.mem_ulimit:
            cmd = ['ulimit -s unlimited ;']
            # allow multiple shell commands to be run in one subprocess
            shell = True
        else:
            cmd = []
            shell = False
        cmd.append(self.xtb_path)
        cmd.append(xyz)
        # set GFN Parameterization
        if self.gfn_version != '2':
            cmd.append('--gfn')
            cmd.append(self.gfn_version)
        # turn on hessian calculation if free energy requested
        cmd.append('--hess')
        # set number of cores
        cmd.append('--parallel')
        cmd.append(self.num_cores)
        # add eletronic temp term
        if self.etemp != '300':
            cmd.append('--etemp')
            cmd.append(self.etemp)
        # write solvent section of cmd
        if self.solvent is not None:
            cmd.append('--gbsa')
            cmd.append(self.solvent)
            if self.solvent_grid != 'normal':
                cmd.append(self.solvent_grid)
        # write charge section of cmd
        if self.charge is not None:
            cmd.append('--chrg')
            cmd.append(self.charge)
        cmd = ' '.join(cmd)
        f = open(out_file, 'w')
        # uses the shell if mem_ulimit = True and waits until
        # subproces is complete. This is required to run the mem_ulimit_cmd
        # and GFN calculation in one command, which is then closed, which
        # minimizes the risk of unrestricting the memory limits.
        sp.call(cmd, stdin=sp.PIPE, stdout=f, stderr=sp.PIPE,
                 shell=shell)
        f.close()
        return out_file

    def energy(self, mol, conformer=-1):
        """
        Calculates the energy of molecule `mol` using xTB.

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
            out_file = self.__write_and_run_command(mol=mol,
                                                    conformer=conformer)
            self.__get_properties(mol=mol,
                                  conformer=conformer,
                                  output_file=out_file)
        finally:
            os.chdir(init_dir)
        return self.total_energies[(mol, conformer)]
