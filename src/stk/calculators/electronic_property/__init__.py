"""
Defines electronic property calculators.

Electronic property calculators calculate the electronic properties of
molecules, namely the dipole moment, electron affinity and ionization
potential. Each :class:`.ElectronicPropertyCalculator` is initialized
with some settings and calculates the electronic properties of a
molecule with :meth:`~.ElectronicPropertyCalculator.dipole_moment`,
:meth:`~.ElectronicPropertyCalculator.electron_affinity` and
:meth:`~.ElectronicPropertyCalculator.ionization_potential`.

.. code-block:: python

    # Electronic property calculators work with any Molecule objects ,
    # such as BuildingBlock, Polymer, Cage etc.
    mol1 = BuildingBlock(...)
    mol2 = Cage(...)
    mol3 = Polymer(...)

    # Create the electronic property calculator.
    gfnxtb = GFNXTBElectronicProperties(gfnxtb_path='/opt/xtb/bin/xtb')
    mol1_ea = gfnxtb.electron_affinity(mol1)
    # We can optionally specify a conformer.
    mol2_ip = gfnxtb.ionization_potential(mol2, conformer=1)
    mol3_dm = gfnxtb.dipole_moment(mol3)

By default, calling any of the :class:`.ElectronicPropertyCalculator`
methods twice on the same molecule will calculate the property a
second time. However, we can use the
:attr:`~.ElectronicPropertyCalculator.use_cache` option to prevent
recalculations when the same molecule and conformer are given to the
same electronic property calculator a second time.

.. code-block:: python

    caching_gfnxtb = GFNXTBElectronicProperties(
        gfnxtb_path='/opt/xtb/bin/xtb',
        use_cache=True
    )
    # Calculate the EA the first time.
    ea1 = caching_gfnxtb.electron_affinity(mol1)
    # The second time, the EA is returned directly from memory, a
    # second calculation in not run.
    ea2 = caching_gfnxtb.electron_affinity(mol1)


Available electronic prperty calculators.
-----------------------------------------



.. _`adding electronic property calculators`:

Extending stk: Making new electronic property calculators.
----------------------------------------------------------

New electronic property calculators can be made by simply making a
class which inherits the :class:`.ElectronicPropertyCalculator` class.
In addition to this, the new class must define
:meth:`~.ElectronicPropertyCalculator.dopole_moment`,
:meth:`~.ElectronicPropertyCalculator.electron_affinity` and
:meth:`~.ElectronicPropertyCalculator.ionization_potential` methods.
Each method must take 2 arguments, a mandatory `mol` argument and an
optional `conformer` argument. The method will then calculate and
return the relevant property. There are no requirements regarding how
it should go about calculating the property. New energy calculators can
be added into the :mod:`.electronic_property_calculators` submodule or
into a new submodule.

"""

from .electronic_property_calculators import *
from .mopac import *
