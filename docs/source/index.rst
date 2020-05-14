.. stk documentation master file, created by
   sphinx-quickstart on Mon Nov 13 11:15:04 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :hidden:
   :caption: stk
   :maxdepth: 2

   Construction Overview <construction_overview>
   Basic Examples <basic_examples>

.. toctree::
   :hidden:
   :caption: Molecules
   :maxdepth: 2

   Building Block <stk.molecular.molecules.building_block>
   Constructed Molecule <stk.molecular.molecules.constructed_molecule>
   Functional Groups <stk.molecular.functional_groups.functional_groups.functional_group>
   Functional Group Factories <stk.molecular.functional_groups.factories.functional_group_factory>
   Key Makers <stk.molecular.key_makers.molecule>
   Reactions <stk.molecular.reactions.reactions.reaction.reaction>
   Reaction Factories <stk.molecular.reactions.factories.reaction_factory>

.. toctree::
   :hidden:
   :caption: Topology Graphs
   :maxdepth: 2

   Polymer <polymer>
   Organic Cage <stk.molecular.topology_graphs.cage.cage>
   Covalent Organic Framework <stk.molecular.topology_graphs.cof.cof>
   Metal Complex <stk.molecular.topology_graphs.metal_complex.metal_complex>
   Host Guest Complex <stk.molecular.topology_graphs.host_guest.complex>
   Macrocycle <stk.molecular.topology_graphs.macrocycle.macrocycle>
   [n]Rotaxane <stk.molecular.topology_graphs.rotaxane.nrotaxane>
   Adding Topology Graphs <stk.molecular.topology_graphs.topology_graph.topology_graph.topology_graph>

.. toctree::
   :hidden:
   :caption: Molecular Databases
   :maxdepth: 2

    Molecule <stk.databases.molecule>
    Constructed Molecule <stk.databases.constructed_molecule>
    Value <stk.databases.value>

.. toctree::
   :hidden:
   :caption: Evolutionary Algorithm
   :maxdepth: 2

   Overview <stk.ea.evolutionary_algorithm.evolutionary_algorithm>
   Basic Example <basic_ea_example>
   Intermediate Example <intermediate_ea_example>
   Fitness Calculators <stk.ea.fitness_calculators.fitness_calculator>
   Fitness Normalizers <stk.ea.fitness_normalizers.fitness_normalizer>
   Selection <stk.ea.selection.selectors.selector>
   Mutation <mutation>
   Crossover <crossover>
   Plotting <plotters>

.. toctree::
   :hidden:
   :caption: Developer Guide
   :maxdepth: 2

   Software Architecture <software_architecture>


============
Introduction
============

GitHub: https://www.github.com/lukasturcani/stk


Installation
------------

To get :mod:`.stk`, you can install it with pip::

    $ pip install stk

Make sure you also install :mod:`rdkit`, which is a dependency::

    $ conda install -c rdkit rdkit


Overview
--------

:mod:`stk` is a Python library which allows construction and
manipulation of complex molecules, as well as automatic
molecular design, and the creation of molecular, and molecular property,
databases.

Molecular Construction
......................

:mod:`.stk` provides tools for constructing a variety of molecular
structures, including organic cages, covalent organic frameworks,
polymers and macrocycles, among others. While additional molecular
structures are always being added, :mod:`stk` provides tools which
allow users to easily specify new kinds of molecular structures,
no matter how simple or complex, in case something they want to build
is not built-in. :mod:`stk` also makes it easy to specify which isomer
of a given structure should be built.

Here is a preview of some of the more complex structures
:mod:`stk` can construct

.. image:: https://i.imgur.com/PvkfoGs.jpg

Automatic Molecular Design
..........................

:mod:`.stk` provides an evolutionary algorithm, which can be used for
discovery of new molecules, which have properties desired by the user.
This evolutionary algorithm works with any molecule :mod:`stk` can
construct, even those defined by users.

Databases
.........

:mod:`.stk` provides tools for the creation of molecular databases,
and storage of molecular property values. :mod:`stk` comes with
support for the creation of MongoDB databases, which can be used to
store molecules constructed by users, or discovered by the
evolutionary algorithm. In addition, property values calculated for
those molecules can also be stored and retrieved from the database.

Usable Defaults
...............

A goal of :mod:`.stk` is to provide simple interfaces and require
minimal information from users to achieve basic and common tasks,
but also allow extensive customization and extension in
order to fulfill complex requirements.

Extensibility and Customization
...............................

Every part of :mod:`.stk` can be extended and customized in user
code, and every user-made extension is indistinguishable from
natively implemented features. This means users can use
:mod:`.stk` to construct new classes of molecules, add new
kinds of molecular databases, and add or customize evolutionary
algorithm operations, all
without looking at :mod:`.stk` source code. All such extensions will
work with the rest of :mod:`.stk` as though they were part of the
library itself.

Documentation and Examples
..........................

Every use-case and extension or customization of :mod:`.stk` has
documentation and examples which will guide users.
:mod:`.stk` is built around abstract base classes, which means
all user extensions to :mod:`.stk` involve creating a new class and
defining, usually, a single method.


Future Releases
---------------

:mod:`.stk` is under active development. Important features in the
future will include:

Molecular Database Visualization
................................

Currently :mod:`stk` allows users to store molecules in a database.
However, a current goal is to develop an app which can be used to
visually inspect and browse the molecules in the database. This should
help multi-member research groups, where a computational scientist
is responsible for suggesting molecules for synthesis, while an
experimentalist actually has to make them. By providing an app to
inspect the database, all collaborators can examine the currently
deposited molecules in real-time.

Distributed Evolutionary Algorithms
...................................

Evolutionary algorithms are very simple to parallelize. It's just a
matter of calculating each fitness value on a separate CPU core.
However, this idea can be taken further. Instead of calculating the
fitness function on a separate CPU core, calculate it on a separate
computer. :mod:`stk` will support fitness functions, which instead
of calculating the fitness value locally, send the molecule over the
internet to a server, so that the server is responsible for calculating
the fitness value, which it then sends back to the client. You will
be able to keep adding servers as your computational requirements
increase, and the load will be distributed fairly between them.
This will allow continuous horizontal scaling of :mod:`stk`
evolutionary algorithms.

More Molecular Structures
.........................

:mod:`.stk` is always being expanded with new molecular structures.
If there is a specific kind of molecule you would like to be able
to construct with :mod:`.stk`, and you don't feel like implementing
it yourself, simply create an issue on
https://github.com/lukasturcani/stk/issues, and describe what you would
like to build.

What Next?
----------

Something you might like to do first, is look at the
`construction overview`_, which can give you a picture of how
:mod:`stk` goes about constructing molecules. It will also introduce
you to the basic concepts and types found within :mod:`stk`. Next, the
`basic examples`_, will allow
you to get a feel for how to use :mod:`stk`. After that, examples of
molecular construction can be seen by looking at the documentation of
the different topology graphs. In general, you will find examples on
how to use a class, in that classes documentation.
Once you are comfortable with construction, you can start looking at
how to deposit and retrieve the molecules from databases. Finally,
:mod:`stk` provides multiple examples on how to use its evolutionary
algorithm.

.. _`construction overview`: construction_overview.html
.. _`basic examples`: basic_examples.html
