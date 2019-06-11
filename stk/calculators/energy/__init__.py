"""
Defines energy calculators.

Energy calculators are objects which calculate the energy of molecules.
Each :class:`.EnergyCalculator` is initialized with some settings and
calculates the energy of a molecule with
:meth:`~.EnergyCalculator.energy`.

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

By default, calling :meth:`~.EnergyCalculator.energy` twice on the
same molecule will calculate the energy a second time. However, we can
use the :attr:`~.EnergyCalculator.use_cache` option to prevent
recalculations when the same molecule and conformer are given to the
same energy calculator a second time.

.. code-block:: python

    caching_mmff = MMFFEnergy(use_cache=True)
    # Calculate the energy the first time.
    energy1 = caching_mmff.energy(mol1)
    # The second time, the energy is returned directly from memory, a
    # second calculation is not run.
    energy2 = caching_mmff.energy(mol1)

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
the new class must define a :meth:`~.EnergyCalculator.energy` method.
The method must take 2 arguments, a mandatory `mol` argument and an
optional `conformer` argument. The method will then calculate and
return the energy. There are no requirements regarding how it should
go about calculating the energy. New energy calculators can be added
into the :mod:`.energy_calculators` submodule or into a new submodule.


"""

from .energy_calculators import *
from .macromodel import *
from .mopac import *
