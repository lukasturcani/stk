"""
Molecular Key Makers
====================

#. :class:`.ConstructedMoleculeKeyMaker`
#. :class:`.MoleculeKeyMaker`

Molecular key makers are used by :mod:`stk` to provide some identity
metric for different molecules. Each molecular key maker provides a
different metric for identity. There are 2 abstract base classes from
which a molecular key maker can be subclassed.

A very important thing to note, is that the abstract base classes of
molecular key makers have their inheritance relationship reversed. For
example, a :class:`.ConstructedMoleculeKeyMaker` IS NOT A
:class:`.MoleculeKeyMaker`, despite the fact that a
:class:`.ConstructedMolecule` IS A
:class:`.Molecule`. This means you cannot use
:class:`.ConstructedMoleculeKeyMaker` where a
:class:`.MoleculeKeyMaker` is required. However a
:class:`.MoleculeKeyMaker` IS A :class:`.ConstructedMoleculeKeyMaker`.
This means that you can use a :class:`.MoleculeKeyMaker` where a
:class:`.ConstructedMoleculeKeyMaker` is required.

More information about each type of molecular key maker can be seen in
the documentation of each of the abstract base classes.

Why is the Inheritance Relationship Reversed?
---------------------------------------------

This section details the reasoning behind the reversed
relationship between the molecular key makers, which may appear
counter-intuitive at first. This is purely a discussion of the
implementation, and users of :mod:`stk` can safely ignore this section,
as long as they follow the guideline on molecular key maker usage
provided above.

The reason the abstract base classes have their relationship reversed
is they define completely different methods. When you look at the
abstract base classes initially, you may think they all define the
same method, namely :meth:`get_key`. However, this is
extremely misleading. Looking at the type signatures of each of these
methods, it is clear that they are not interchangeable.

For example, :class:`.MoleculeKeyMaker.get_key` takes a
:class:`.Molecule` object and returns a key. However,
:class:`.ConstructedMoleculeKeyMaker` takes a
:class:`.ConstructedMolecule` object and returns a key. This alone
makes them completely unrelated methods. The easiest way to explain
why this is the case is through an example.

For the sake of the example, lets define a really simple
:class:`.ConstructedMoleculeKeyMaker` subclass. Note that because
were are defining a :class:`.ConstructedMoleculeKeyMaker`, we can use
the entire interface of :class:`.ConstructedMolecule`. This subclass
simply looks at the number of building blocks the constructed molecule
has, and uses that as a key

.. code-block:: python

    import stk

    class MyKeyMaker(stk.ConstructedMoleculeKeyMaker):
        def get_key_name(self):
            return 'MyConstructedMoleculeKeyMaker'

        def get_key(self, constructed_molecule):
            building_blocks = tuple(
                constructed_molecule.get_building_blocks()
            )
            return len(building_blocks)

Already, it should be obvious why a
:class:`.ConstructedMoleculeKeyMaker` IS NOT A
:class:`.MoleculeKeyMaker`. If we give :class:`MyKeyMaker` to a
function, which expects a :class:`.MoleculeKeyMaker`, at some point,
a :class:`.Molecule` will get passed to the :meth:`get_key` method
of :class:`MyKeyMaker`. Once this happens,
:meth:`~.ConstructedMolecule.get_building_blocks` will be called.
However, a :class:`.Molecule` instance is not guaranteed to have this
method. If a plain :class:`.Molecule` instance was used, or a
was used, a runtime error will occur, because such as method does not
exist on plain :class:`.Molecule` objects.

Therefore, you cannot give a :class:`.ConstructedMoleculeKeyMaker`
where a :class:`.MoleculeKeyMaker` is expected, because eventually you
will pass a :class:`.Molecule`, which is not also a
:class:`.ConstructedMolecule` to the
:class:`.ConstructedMoleculeKeyMaker` and that will cause a runtime
error if the :class:`.ConstructedMoleculeKeyMaker` makes use of methods
found on a :class:`.ConstructedMolecule` and not a plain
:class:`.Molecule`.

However, the opposite is not true. Because any
:class:`.MoleculeKeyMaker` will accept any :class:`.Molecule`,
including a :class:`.ConstructedMolecule`, you can pass a
:class:`.MoleculeKeyMaker` where a
:class:`.ConstructedMoleculeKeyMaker` is required, because even if you
do pass a :class:`.ConstructedMolecule` instance to a
:class:`.MoleculeKeyMaker`, that is guaranteed to be safe. As a result,
a :class:`.MoleculeKeyMaker` is a type of
:class:`.ConstructedMoleculeKeyMaker`

"""

from .molecule import *
from .constructed_molecule import *
