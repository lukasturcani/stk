Calculators
===========

``stk`` makes extensive use of calculator objects. Calculators are a
very general pattern for performing operations on molecules. All
calculators are objects with a special method, used to perform the
calculation. For example, an optimizer has a special method called
:meth:`~.Optimizer.optimize`, which optimizes a molecule.

Other types of calculators will have different special methods used to
perform calculations on molecules. For example, an
:class:`.EnergyCalculator` will define a
:meth:`~.EnergyCalculator.get_energy` method. This method calculates
the energy of a molecule.

All evolutionary algorithm operations are also implemented through
calculators, take for example a :class:`.Mutator` such as
:class:`.RandomTopologyGraph`

.. code-block:: python

    import stk

    topology_graphs = [
        stk.cage.EightPlusTwelve(),
        stk.cage.FourPlusSix(),
        stk.cage.TwentyPlusThirty()

    random_top = RandomTopologyGraph(topology_graphs)

    cage = stk.ConstructedMolecule([bb1, bb2], stk.cage.SixPlusNine())

    # Perform a mutation.
    mutant = random_top.mutate(cage)

A :class:`.Mutator` is a calculator which implements a
:meth:`~.Mutator.mutate` method. This method takes a molecule and
returns a mutant, *i.e.* a modified version, of that molecule.
In this example, the ``random_top`` is a :class:`.Mutator`, which
replaces the topology of a molecule with a new one, in order to
generate the mutant. In this case, the new topologies would be one of
:class:`.EightPlusTwelve`, :class:`.FourPlusSix` or
:class:`.TwentyPlusThirty`.

Calculators often support additional options to modify their behaviour.
For example, calculators of type :class:`.Optimizer` or
:class:`.EnergyCalculator` support caching. This means that if the
same molecule and conformer is supplied to the calculator, it will not
run the optimization or energy calculation again, it will return the
previously calculated value.

Combining Calculators
.....................

Calculators can be combined to create complex behaviour on the fly.
For example, we may wish to make a 2 step optimization process. First,
we perform an optimization using the MMFF force field and the run
a MacroModel molecular dynamics conformer search with
:class:`.MacroModelMD`. The obvious way to do
this is to run the two optimizers in sequence

.. code-block:: python

    mmff = stk.MMFF()
    macromodel = stk.MacroModelMD('/opt/schrodinger2017-4')
    mmff.optimize(mol)
    macromodel.optimize(mol)

However, there is a better way! We can use an optimizer called
:class:`.OptimizerSequence`. The :meth:`~.OptimizerSequence.optimize`
method of this optimizer calls the :meth:`~.Optimizer.optimize` methods
of the optimizers it was initialized with

.. code-block:: python

    opt_sequence = OptimizerSequence(mmff, macromodel)
    # Optimize with mmff and then with macromodel.
    opt_sequence.optimize(mol)

This pattern is quite common and powerful. For example, we can take
three different :class:`.Mutator` objects. Each of these defines
a different :meth:`~.Mutator.mutate` method. We want to apply one of
these mutations at random. We can simply use :class:`.RandomMutation`

.. code-block:: python

    random_bb = stk.RandomBuildingBlock(...)
    similar_bb = stk.SimilarBulidingBlock(...)
    random_topology = stk.RandomTopology(...)

    random_mutation = stk.RandomMutation(
        random_bb,
        similar_bb,
        random_topology
    )

    # Use one of the mutate() methods of random_bb, similar_bb and
    # random_topology at random.
    mutant1 = random_mutation.mutate(mol)
    # The next call use a different mutation to the call above.
    mutant2 = random_mutation.mutate(mol)

The :meth:`.RandomMutation.mutate` method randomly selects a
:class:`.Mutator` it was initialized with to carry out the mutation
on its behalf.

Making New Calculators
......................

New calculators can be added very simply and they can be defined in
user code.

A simple example of adding a new calculator in user code.

.. code-block:: python

    import stk

    class NewEnergyCalculator(stk.EnergyCalculator):
        def get_energy(self, mol):
            return 15

    energy_calc = NewEnergyCalculator()
    # Calculate the energy with the new calculator.
    energy_calc.get_energy(mol)


Note that you can also modify the behaviour of existing calculators in
your code, see the :doc:`cookbook`.

Saving Calculator Results
.........................

If you want to save the results which your calculator found you should
set ``use_cache=True`` and then you can simply use :mod:`pickle`
to dump and load the calculator object.
