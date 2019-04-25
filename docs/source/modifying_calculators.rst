Modifying ``stk`` Calculators
=============================

This guide takes the form a cookbook to give examples of how ``stk``
calculators can be extended in user code to add functionality. The
basic principle is the use of decorators to extend functionality.

Adding timing information to calculator calls.
==============================================

Problem
.......

You want to see how long different optimization settings take on
different molecules.


Solution
........

Create a decorator which wraps a function and times it.

.. code-block:: python

    import stk
    import time
    from functools import wraps

    def timed_fn(fn):
        """
        Adds timing info to functions.

        """

        @wraps(fn)
        def inner(self, mol, conformer=-1):
            start = time.time()
            r = fn(*args, **kwargs)
            print(
                f'{fn.__name__} takes {time.time()-start} '
                f'seconds on {mol.name}.'
            )
            return r

        return inner

    # Create some calculators.
    mmff = stk.MMFF()
    energy_calculator = stk.UFFEnergy()

    # Add timing info to the calculators.
    mmff.optimize = timed_fn(mmff.optimize)
    energy_calculator.energy = timed_fn(energy_calculator.energy)

    mol = stk.StructUnit(...)
    # These calls will print how long the function execution took.
    mmff.optimize(mol)
    energy_calculator.energy(mol)


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

Change the default argument with a decorator.

.. code-block:: python

    import stk
    from functools import wraps

    def set_default_conformer(optimize, conformer):
        """
        Sets default conformer of `optimize` to `conformer`.

        """

        @wraps(fn)
        def inner(self, mol):
            return optimize(self, mol, conformer)

        return inner

    mmff = MMFF()
    # The conformer optimized by this optimizer will be 10.
    mmff.optimize = set_default_conformer(mmff.optimize, 10)

    # Make a population of molecules.
    pop = stk.Population(...)

    # Optimize conformer 10 of all the molecules.
    pop.optimize(mmff)
