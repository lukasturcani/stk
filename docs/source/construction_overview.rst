Construction Overview
=====================

Introduction
------------

The key idea behind :mod:`stk` is that the construction of a molecule can
be broken down into two fundamental pieces of information, its
building blocks and its topology graph. The building blocks of a
molecule are molecules, or molecular fragments, which are used for
construction. The smallest possible building block is a single atom
and constructed molecules can become the building blocks of other
constructed molecules. The topology graph is an abstract representation
of a constructed molecule. The nodes of the graph represent the
positions of the building blocks and the edges of the graph represent
which building blocks have bonds formed between them during
construction.

To use :mod:`stk` you only have to choose which building blocks and
topology graph you wish to use and :mod:`stk` will take care of everything
else, take for example the construction of a linear polymer

.. code-block:: python

    import stk

    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            building_blocks=(
                stk.BuildingBlock(
                    smiles='Brc1ccc(Br)cc1',
                    functional_groups=[stk.BromoFactory()],
                ),
                stk.BuildingBlock(
                    smiles='BrC#CBr',
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            repeating_unit='AB',
            num_repeating_units=2,
        )
    )
    # You can write the molecule to a file if you want to view it.
    polymer.write('polymer.mol')

which will produce:

.. image:: https://i.imgur.com/jQ7s8qp.png

Because the topology graph is an idealized representation of the
constructed molecule, the bonds formed during construction often have
unrealistic lengths. This means that constructed molecules will need to
undergo structure optimization. There is no single correct way to go
about this, because the appropriate methodology for structure
optimization will depend on various factors, such as the nature of the
constructed molecule, the desired accuracy, and time constraints.
In addition, there are countless options already available,
be it Python libraries such as :mod:`rdkit` or :mod:`ase`, or
some sort of computational chemistry software. Since
:mod:`stk` cannot hope to provide a good solution to this problem,
it does try to make it easy for you to convert an
:mod:`stk` :class:`.Molecule` into whatever format you need to make
use of other software. This means you can access atoms and
bonds with :meth:`.Molecule.get_atoms` and :meth:`.Molecule.get_bonds`,
you can convert any :mod:`stk` :class:`.Molecule` into an
:mod:`rdkit` molecule with :meth:`.Molecule.to_rdkit_mol` or you
can write it to a file with :meth:`.Molecule.write`.

.. figure:: https://i.imgur.com/UlCnTj9.png
    :align: center

    The general construction workflow of ``stk``.

Topology Graph
--------------

The abstraction provided by the topology graph has a number of
powerful benefits. Firstly, because every vertex is responsible for the
placement of a building block, it is extremely easy to construct
different structural isomers of the constructed molecule. The vertex
can be told to perform different transformations on the building block,
so that its orientation in the constructed molecule changes. For the
end user, selecting the transformation from a set of
predefined ones is easy. Also, since the transformation is restricted
to a single building block on a single vertex, it easy for developers
to define.

The second major benefit of the topology graph is that the vertices and
edges can hold additional state useful for the construction of a
molecule. An example of this is in the construction of different
structural isomers, but another can be seen in the construction of
periodic systems. For example, :mod:`stk` allows you to construct
covalent organic frameworks. With the topology graph this is trivial
to implement, simply label some of the edges as periodic and they
will construct periodic bonds instead of regular ones.

Thirdly, the topology graph allows users to
easily modify the construction of molecules by placing different
building blocks on different vertices.

Finally, the topology graph breaks down the construction of
a molecule into independent steps. Each vertex
represents a single, independent operation on a building block while
each edge represents a single, independent operation on a collection
of building blocks. As a result, each vertex and edge represents a
single operation, which can be executed in parallel. This allows
:mod:`stk` to scale efficiently to large topology graphs and take
advantage of multiple cores even during the construction of a single
molecule.

Building Blocks
---------------

Building blocks in :mod:`stk` are molecules, or molecular fragments,
which are placed on the nodes of the :class:`.TopologyGraph`. After
building blocks are placed the nodes, they are connected to
each other through a :mod:`~.reaction.reaction`. :mod:`stk` supports
multiple reactions and users can define their own. Reactions can add or
remove atoms and bonds between building blocks, which are connected by
edges in the topology graph.

When it comes to reactions, an important question, which must be
addressed, is, which atoms of a :class:`.BuildingBlock` are modified
by a :mod:`~.reaction.reaction`? In :mod:`stk`, the answer to this is
a :mod:`~.functional_groups.functional_group`. When a user of
:mod:`stk` creates a :class:`.BuildingBlock`, they also specify which
functional groups are present in the :class:`.BuildingBlock`. This
lets :mod:`stk` know which atoms the user intends to transform during
construction.

There are many different types of
:mod:`~.functional_groups.functional_group` present
in :mod:`stk`, for example, :class:`.Bromo`, :class:`.Alcohol` or
:class:`.Aldehyde`. When a user creates a :class:`.BuildingBlock`,
they can specify multiple functional groups at at time using
a :mod:`~.functional_group_factory`. A
:mod:`~.functional_group_factory` finds all the functional groups of a
specific type, and adds them to the :class:`.BuildingBlock`. For
example, if we want to create
a :class:`.BuildingBlock`, and you want to react its bromo groups
during construction, you can use a :class:`.BromoFactory`.

.. code-block:: python

    import stk

    building_block = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])

In the example above, ``building_block`` will have two
:class:`.Bromo` functional groups. When ``building_block`` is used
for construction, it is the atoms held by the :class:`.Bromo`
groups, which will be modified. If we have a building block with
aldehyde functional groups, we could have used an
:class:`.AldehydeFactory`.

.. code-block:: python

    building_block2 = stk.BuildingBlock('O=CCC=O', [stk.AldehydeFactory()]

Finally, if we had a mix of functional groups, we could have used
a mix of factories

.. code-block:: python

    building_block2 = stk.BuildingBlock(
        smiles='O=CCCBr',
        functional_groups=[stk.AldehydeFactory(), stk.BromoFactory()],
    )


Based on the specific functional groups found on an
edge of the :class:`.TopologyGraph`, :mod:`stk` will select an
appropriate :mod:`~.reaction.reaction` to join them. You can also force
:mod:`stk` to use a different :mod:`~.reaction.reaction` of your
choosing, which is covered in the `basic examples`_.

.. _`basic examples`: basic_examples.html
