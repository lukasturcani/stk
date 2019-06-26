"""
Defines optimizers.

Optimizers are objects used to optimize molecules. Each optimizer is
initialized with some settings and used to optimize a molecule
with :meth:`~.Optimizer.optimize`.

.. code-block:: python

    import rdkit.Chem.AllChem as rdkit
    mol = StructUnit2.smiles_init('NCCCN', ['amine'])
    mmff = MMFF()
    mmff.optimize(mol)

    # Optionally, a conformer can be provided.
    mmff.optimize(mol, conformer=2)

    # Optimizers also work with MacroMolecule objects.
    polymer = Polymer([mol], Linear('A', [0], n=3))
    etkdg = ETKDG()
    etkdg.optimize(polymer)

Sometimes it is desirable to chain multiple optimizations, one after
another. For example, before running an optimization, it may be
desirable to embed a molecule first, to generate an initial structure.
:class:`.OptimizerSequence` may be used for this.

.. code-block:: python

    # Create a new optimizer which chains the previously defined
    # mmff and etkdg optimizers.
    optimizer_sequence = OptimizerSequence(etkdg, mmff)

    # Run each optimizer in sequence.
    optimizer_sequence.optimize(polymer)

By default, running :meth:`.Optimizer.optimize` twice on the same
molecule will perform an optimization a second time on a molecule. If
we want to skip optimizations on molecules which have already been
optimized we can use the :attr:`~.Optimizer.use_cache` flag.

.. code-block:: python

    caching_etkdg = ETKDG(use_cache=True)
    # First optimize call runs an optimization.
    caching_etkdg.optimize(polymer)
    # Second call does nothing.
    caching_etkdg.optimize(polymer)

Caching is done on a per :class:`.Optimizer` basis. Just because the
molecule has been cached by one :class:`.Optimizer` instance does not
mean that a different :class:`.Optimizer` instance will no longer
optimize the molecule.

Available optimizers.
---------------------

#. :class:`.NullOptimizer`
#. :class:`.MMFF`
#. :class:`.UFF`
#. :class:`.ETKDG`
#. :class:`.XTB`
#. :class:`.MacroModelForceField`
#. :class:`.MacroModelMD`
#. :class:`.MOPAC`
#. :class:`.OptimizerSequence`
#. :class:`.CageOptimizerSequence`
#. :class:`.TryCatchOptimizer`
#. :class:`.RaisingOptimizer`

.. _`adding optimizers`:

Extending stk: Making new optimizers.
-------------------------------------

New optimizers can be made by simply making a class which inherits the
:class:`.Optimizer` class. In addition to this, the new class must
define a :meth:`~.Optimizer.optimize` method. The method must take 2
arguments a mandatory `mol` argument and an optional `conformer`
argument. :meth:`~.Optimizer.optimize` will take a molecule and change
its structure in whatever way it likes. Beyond this there are no
requirements. New optimizers can be added into the :mod:`.optimizers`
submodule or into a new submodule.

"""

from .optimizers import *
from .macromodel import *
from .mopac import *
