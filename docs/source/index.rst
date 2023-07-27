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
   Video Tutorials <video_tutorials>

.. toctree::
   :hidden:
   :caption: Molecules
   :maxdepth: 2

   Building Block <_autosummary/stk.BuildingBlock>
   Constructed Molecule <_autosummary/stk.ConstructedMolecule>
   Functional Groups <functional_groups>
   Functional Group Factories <functional_group_factories>
   Key Makers <key_makers>
   Reactions <reactions>
   Reaction Factories <reaction_factories>
   Writers <writers>

.. toctree::
   :hidden:
   :caption: Topology Graphs
   :maxdepth: 2

   Polymers <polymer>
   Small Molecules <small>
   Organic & Metal-Organic Cages <cage>
   Covalent Organic Frameworks <cof>
   Metal Complexes <metal_complex>
   Host Guest Complex <_autosummary/stk.host_guest.Complex>
   Macrocycle <_autosummary/stk.macrocycle.Macrocycle>
   [n]Rotaxane <_autosummary/stk.rotaxane.NRotaxane>
   Adding Topology Graphs <_autosummary/stk.TopologyGraph>
   Optimizers <optimizers>

.. toctree::
   :hidden:
   :caption: Molecular Databases
   :maxdepth: 2

    Molecule <molecule_databases>
    Constructed Molecule <constructed_molecule_databases>
    Value <value_databases>

.. toctree::
   :hidden:
   :caption: Evolutionary Algorithm
   :maxdepth: 2

   Overview <_autosummary/stk.EvolutionaryAlgorithm>
   Basic Example <basic_ea_example>
   Intermediate Example <intermediate_ea_example>
   Fitness Calculators <fitness_calculators>
   Fitness Normalizers <fitness_normalizers>
   Selection <selection>
   Mutation <mutation>
   Crossover <crossover>
   Plotting <plotters>

.. toctree::
   :hidden:
   :caption: Modules
   :maxdepth: 1

   Modules <modules.rst>


============
Introduction
============

| GitHub: https://www.github.com/lukasturcani/stk
| Discord: https://discord.gg/zbCUzuxe2B


Installation
------------

To get :mod:`.stk`, you can install it with pip::

  $ pip install stk

Overview
--------

:mod:`stk` is a Python library which allows construction and
manipulation of complex molecules, as well as automatic
molecular design, and the creation of molecular, and molecular property,
databases.

Molecular Construction
......................

:mod:`.stk` provides tools for constructing a variety of molecular
structures, including organic and metal-organic cages,
covalent organic frameworks, polymers and metal complexes, to name a
few. While additional molecular
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

Molecular Database Visualization
................................

.. image:: https://i.imgur.com/8MCBUGZ.png

:mod:`.stk` has a sibling project called `stk-vis`_, which is a
cross-platform application, which lets you connect to a database
created by :mod:`.stk`, and view the molecules deposited into it.
`stk-vis`_ also shows you any molecular properties you deposited,
and lets you browse the building blocks of any constructed molecules.
It's ideal for collaboration between multiple people, because one
person can deposit a molecule into the database, and another person
can immediately see and examine it. `stk-vis`_ is a stand-alone
application and does not require coding or :mod:`.stk` to use.

Features of `stk-vis`_ include

* 3D interactive molecular rendering.
* Image of the 2D molecular graph.
* Tabulation of any molecular properties deposited into the
  database.
* Sorting molecules according to property values to quickly find ones
  with the best and worst properties.
* You can list the building blocks of any constructed molecules.
  If those building blocks are also constructed molecules, you can
  list their building blocks as well, and so on.
* Writing molecules to files.

You read more about `stk-vis`_ here:

  https://github.com/lukasturcani/stk-vis

.. _`stk-vis`: https://github.com/lukasturcani/stk-vis

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

:mod:`.stk` is under active development. You can get alerted
when a new release comes out by going to the `GitHub page`_ and
click on the ``watch`` button in the top right corner. Then select
``Releases only`` from the dropdown menu.

.. _`GitHub page`: https://github.com/lukasturcani/stk

Important features in the future will include:


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
