"""
Defines energy calculators.

Energy calculators are objects which calculate the energy of molecules.
Each :class:`EnergyCalculator` is initialized with some settings and
calculates the energy of a molecule with
:meth:`~EnergyCalculator.energy`.

.. code-block:: python

    # Energy calculators work with any Molecule objects, such as
    # StructUnit, Polymer, Cage, etc.
    mol1 = StructUnit(...)
    mol2 = Cage(...)
    mol3 = Polymer(...)

    # Create the energy calculator.
    mmff = MMFFEnergy()

    # Calculate energies of various molecules.
    mol1_energy = mmff.energy(mol1)
    # We can optionally specify a conformer.
    mol2_energy = mmff.energy(mol2, conformer=1)
    mol3_energy = mmff.energy(mol3)

By default, calling :meth:`~EnergyCalculator.energy` twice on the
will calculate the energy a second time. However, we can use the
:attr:`~EnergyCalculator.use_cache` option to prevent recalculations
when the same molecule and conformer are given to the same energy
calculator a second time.

.. code-block:: python

    caching_mmff = MMFFEnergy(use_cache=True)
    # Calculate the energy the first time.
    energy1 = caching_mmff.energy(mol1)
    # The second time, the enegy is returned directly from memory, a
    # second calculation is not run.
    energy2 = caching_mmff.enegy(mol2)

.. _`adding energy calculators`:

Extending stk: Making new energy calculators.
---------------------------------------------

New energy calculators can be made by simply making a class which
inherits the :class:`EnergyCalculator` class. In addition to this,
the new class must define a :meth:`~EnergyCalculator.energy` method.
The method must take 2 arguments, a mandatory `mol` argument and an
optional `conformer` argument. The method will then calculate and
return the energy. There are no requirements regarding how it should
go about calculating the energy.

"""

import rdkit.Chem.AllChem as rdkit
import logging
import re
from functools import wraps


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

        products = zip(self.products,  self.product_conformers)
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


class GFNXTB(EnergyCalculator):
    """
    Uses GFN-xTB to calculate energies.

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

    def __init__(self, gfnxtb_path, output_dir=None, free=False, num_cores=1,
                 energy_type='total', solvent=None, charge=None
                 ):
        """
        Initializes a :class:`GFNXTB` instance.

        Parameters
        ----------
        gfnxtb_path : :class:`str`
            The path to the GFN-xTB or GFN2-xTB executable.

        output_dir : :class:`str`, optional
            The name of the directory into which files generated during
            the optimization are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        free : :class:`bool`, optional
            If ``True`` :meth:`optimize` will perform a numerical hessian
            calculation on the optimized structure to give Free energy also.

        opt_level : :class:`int`, optional
            Optimization level to use.
            -2 =  , -1 = , 0 = , 1 = , 2 = .

        num_cores : :class:`int`
            The number of cores for GFN-xTB to use.

        energy_type : :class:`str`, optional
            Type of energy to extract.
            `total` : total energy of system
            `free` : total free energy of system. Reqiures free=True

        solvent : :class:`str`, optional
            XXXX

        charge : :class:`str`, optional
            XXXX
        """
        self.gfnxtb_path = gfnxtb_path
        self.output_dir = output_dir
        self.free = free
        self.num_cores = str(num_cores)
        self.energy_type = energy_type
        self.solvent = solvent
        if self.solvent is not None:
            self.valid_solvent()
        self.charge = charge

    def extract_energy(self, output_file):
        """
        Extract desired energy from GFN2-xTB output file.

        Format depends on version. Works with version 190418.


        Obtained results (in a.u.):
            - free energies (FE)
            - absolute energy (TE)
            - SCC energy (SCE)
            - G (GT)
            - H (HT)

        """
        print('CHECK OTHER ENERGY TYPES')
        if self.free is False:
            if self.energy_type in ['total']:
                raise(f'{self.energy_type} requires hessian calculation - free=True')

        for line in open(output_file, 'r'):
            # regex:
            nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
            # free energy in a.u.
            if '| TOTAL FREE ENERGY' in line and self.energy_type == 'total':
                FE_au = nums.search(line.rstrip()).group(0)
                energy_au = FE_au
                break
            if '| TOTAL ENERGY' in line and self.energy_type == 'free':
                TE_au = nums.search(line.rstrip()).group(0)
                energy_au = TE_au
                break
            if '| TOTAL ENTHALPY' in line and self.energy_type == 'enthalpy':
                HT_au = nums.search(line.rstrip()).group(0)
                energy_au = HT_au
                break
        return energy_au

    def energy(self, mol, conformer=-1):
        ' get from optimizer code once implemented'
