.. _cookbook:

========
Cookbook
========

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
            self._optimizer = optimizer
            super().__init__(use_cache=use_cache)

        def optimize(self, mol):
            start = time.time()
            r = self._optimizer.optimize(mol)
            opt_name = self.optimizer.__class__.__name__
            print(
                f'{opt_name} takes {time.time()-start} '
                f'seconds on {mol}.'
            )
            return r


    class TimedEnergyCalculator(stk.EnergyCalculator):
        def __init__(self, energy_calculator, use_cache=False):
            self._energy_calculator = energy_calculator
            super().__init__(use_cache=use_cache)

        def energy(self, mol):
            start = time.time()
            energy = self._energy_calculator.energy(mol, conformer)
            calc_name = self.energy_calculator.__class__.__name__
            print(
                f'{calc_name} takes {time.time()-start} '
                f'seconds on {mol}.'
            )
            return energy

    # Create some timed calculators.
    mmff = stk.MMFF()
    timed_mmff = TimedOptimizer(mmff)

    uff_energy = stk.UFFEnergy()
    timed_uff_energy = TimedEnergyCalculator(uff_energy)

    mol = stk.BuildingBlock('NCCCN')
    # These calls will print how long the function execution took.
    timed_mmff.optimize(mol)
    timed_uff_energy.energy(mol)
