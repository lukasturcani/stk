Extending Calculators
=============================

This guide takes the form a cookbook. Common problems are presented
followed by a suggested solution.

Adding timing information to calculator calls.
==============================================

Problem
.......

You want to see how long different optimization settings take on
different molecules.

Solution
........

Create a new optimizer class. The optimizers of this class will take
take the optimizer you want to time as an initialization parameter.

.. code-block:: python

    import stk
    import time


    class TimedOptimizer(stk.Optimizer):
        def __init__(self, optimizer, use_cache=False):
            self.optimizer = optimizer
            super().__init__(use_cache=use_cache)

        def optimize(self, mol, conformer=-1):
            start = time.time()
            r = self.optimizer.optimize(mol, conformer)
            opt_name = self.optimizer.__class__.__name__
            print(
                f'{opt_name} takes {time.time()-start} '
                f'seconds on {mol.name}.'
            )
            return r


    class TimedEnergyCalculator(stk.EnergyCalculator):
        def __init__(self, energy_calculator, use_cache=False):
            self.energy_calculator = energy_calculator
            super().__init__(use_cache=use_cache)

        def energy(self, mol, conformer=-1):
            start = time.time()
            energy = self.energy_calculator.energy(mol, conformer)
            calc_name = self.energy_calculator.__class__.__name__
            print(
                f'{calc_name} takes {time.time()-start} '
                f'seconds on {mol.name}.'
            )
            return energy

    # Create some timed calculators.
    mmff = stk.MMFF()
    timed_mmff = TimedOptimizer(mmff)

    uff_energy = stk.UFFEnergy()
    timed_uff_energy = TimedEnergyCalculator(uff_energy)

    mol = stk.StructUnit(...)
    # These calls will print how long the function execution took.
    timed_mmff.optimize(mol)
    timed_uff_energy.energy(mol)


The initial instinct may be to use a decorator to decorate
:meth:`optimize` and :meth:`energy`. However, this approach is not
not compatible with multiprocessing as the resulting calculators
cannot be pickled.

Specifying a conformer with :meth:`.Population.optimize`.
=========================================================

Problem
.......

When using :meth:`.Population.optimize`, you cannot specify which
conformer of the molecules in the population is used. The value
in the `conformer` argument of :meth:`~.Optimizer.optimize` defaults
to ``-1``.

Solution
--------

Create a new optimizer class which allows you to set which conformer
to use on initialization.

.. code-block:: python

    import stk


    class ConformerOptimizer(stk.Optimizer):
        def __init__(self, optimizer, conformer, use_cache=False):
            self.optimizer = optimizer
            self.conformer = conformer
            super().__init__(use_cache=use_cache)

        def optimize(self, mol, conformer=-1):
            return self.optimizer.optimize(mol, self.conformer)


    mmff = MMFF()
    # The conformer optimized by this optimizer will be 10.
    mmff_conf_10 = ConformerOptimizer(mmff, 10)

    # Make a population of molecules.
    pop = stk.Population(...)

    # Optimize conformer 10 of all the molecules.
    pop.optimize(mmff_conf_10)


Once again, a decorator cannot be used because it is incompatible with
pickle and multiprocessing. Here a new optimizer class is defined.
Notice that the signature of :meth:`optimize` is unchanged, however
the value passed to the conformer argument will be ignored. Instead
the conformer set during initialization of :class:`ConformerOptimizer`
will be used.
