"""
Defines energy calculators.

Energy calculators are objects which calculate the energy of molecules.
Each :class:`.EnergyCalculator` is initialized with some settings and
calculates the energy of a molecule with
:meth:`~.EnergyCalculator.get_energy`.

.. code-block:: python

    import stk

    # Energy calculators work with any Molecule objects, such as
    # BuildingBlock or ConstructedMolecule.
    mol1 = stk.BuildingBlock('[Br]CC[Br]', ['bromine'])

    chain = stk.polymer.Linear('A', [0], 12)
    mol2 = stk.ConstructedMolecule([mol1], chain)

    # Create the energy calculator.
    mmff = MMFFEnergy()

    # Calculate energies of various molecules.
    mol1_energy = mmff.get_energy(mol1)
    mol2_energy = mmff.get_energy(mol2)

By default, calling :meth:`~.EnergyCalculator.get_energy` twice on the
same molecule will calculate the energy a second time. However, we can
use the :attr:`~.EnergyCalculator.use_cache` option to prevent
recalculations when the same molecule is given to the same energy
calculator a second time

.. code-block:: python

    caching_mmff = MMFFEnergy(use_cache=True)
    # Calculate the energy the first time.
    energy1 = caching_mmff.get_energy(mol1)
    # The second time, the energy is returned directly from memory, a
    # second calculation is not run.
    energy2 = caching_mmff.get_energy(mol1)

Available energy calculators.
-----------------------------

#. :class:`.MMFFEnergy`
#. :class:`.UFFEnergy`
#. :class:`.MacroModelEnergy`
#. :class:`.MOPACEnergy`
#. :class:`.FormationEnergy`
#. :class:`.XTBEnergy`
#. :class:`.XTBFreeEnergy`


.. _`adding energy calculators`:

Extending stk: Making new energy calculators.
---------------------------------------------

New energy calculators can be made by simply making a class which
inherits the :class:`.EnergyCalculator` class. In addition to this,
the new class must define a :meth:`~.EnergyCalculator.get_energy`
method. The method must take 1 parameter, `mol`. The method will then
calculate and return the energy. There are no requirements regarding
how it should go about calculating the energy. New energy calculators
can be added into the :mod:`.energy_calculators` submodule or into a
new submodule.


"""

from .energy_calculators import *
from .macromodel import *
from .mopac import *
