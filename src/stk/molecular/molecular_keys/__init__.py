"""
Molecular Keys
==============

#. :class:`.BuildingBlockKey`
#. :class:`.ConstructedMoleculeKey`
#. :class:`.MoleculeKey`

Molecular keys are used by :mod:`stk` to provide some identity metric
for different molecules. Each class molecular key provides a different
metric for identity. There are 3 abstract base classes from which a
molecular key can be subclassed. However, each molecular key will
subclass at most 1 of them.

An very important thing to note is that the abstract base classes of
molecular keys are not related by inheritance. For example,
a :class:`.BuildingBlockKey` IS NOT A :class:`.MoleculeKey`, despite
the fact that a :class:`.BuildingBlock` IS A :class:`.Molecule`.
This means you cannot use :class:`.BuildingBlockKey` where a
:class:`.MoleculeKey` is required.

More information about each type of molecular key can be seen in the
documentation of each of the abstract base classes.

Why Are the Abstract Base Classes Not Related by Inheritance?
-------------------------------------------------------------

This section details the reasoning behind the lack of inheritance
relationship between the molecular keys, which may appear
counter-intuitive at first. This is purely a discussion of the
implementation, and users of :mod:`stk` can safely ignore this section,
as long as they follow the guideline on molecular key usage provided
above.

The reason the abstract base classes are not related is because
they define completely different methods. When you look at the
abstract base classes initially, you may think they all define the
same method, namely :meth:`get_key` method. However, this is
extremely misleading. Looking at the type signatures of each of these
methods, it is clear that they are not interchangeable.

For example, :class:`.MoleculeKey.get_key` takes a
:class:`.Molecule` object and returns a key. However,
:class:`.BuildingBlockKey` takes a :class:`.BuildingBlock` object
and returns a key. This alone makes them completely unrelated
methods. The easiest way to explain why this is the case is through
an example.

For the sake of the example, lets define a really simple
:class:`.BuildingBlockKey` subclass. Note that because
were are defining a :class:`.BuildingBlockKey`, we can use the
entire interface of :class:`.BuildingBlock`. This simple subclass
simply looks at the number of functional groups the building block
has, and uses that as a key

.. code-block:: python

    import stk

    class MyBuildingBlockKey(stk.BuildingBlockKey):
        def get_name(self):
            return 'MyBuildingBlockKey'

        def get_key(self, building_block):
            return building_block.get_num_functional_groups()

Already it should be obvious why a :class:`.BuildingBlockKey` IS NOT
A :class:`.MoleculeKey`. If we give provide :class:`MyBuildingBlockKey`
to a function which expects a :class:`.MoleculeKey`, at some point
a :class:`.Molecule` will get passed to the :meth:`get_key` method
of :class:`MyBuildingBlockKey`. Once this happens,
:meth:`.~BuildingBlock.get_num_functional_groups` will be called.
However, a :class:`.Molecule` instance is not guaranteed to have this
method. If a plain :class:`.Molecule` instance was used or a
:class:`.ConstructedMolecule` was used, a runtime error will occur,
because such as method does not exist.

Therefore, you cannot give a :class:`.BuildingBlockKey` where a
:class:`.MoleculeKey` is expected, because eventually you will pass a
:class:`.Molecule`, which is not also a :class:`.BuildingBlock` to the
:class:`.BuildingBlockKey` and that will cause a runtime error if the
:class:`.BuildingBlockKey` makes use of methods found on a
:class:`.BuildingBlock` and not a plain :class:`.Molecule`.

"""

from .molecule_keys import *
from .building_block_keys import *
from .constructed_molecule_keys import *
