Construction Overview
=====================

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
                stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
                stk.BuildingBlock('BrCNCBr', [stk.BromoFactory()]),
            ),
            repeating_unit='ABBBA',
            num_repeating_units=2
        )
    )
    # You can write the molecule to a file if you want to view it.
    polymer.write('polymer.mol')

which will produce

.. image:: https://i.imgur.com/XmKRRun.png


Because the topology graph is an idealized representation of the
constructed molecule, the bonds formed during construction often have
unrealistic lengths. This means that constructed molecules will need to
undergo structure optimization. There is no single correct way to go
about this, because the appropriate methodology for structure
optimization will depend various factors, such as the nature of the
constructed molecule, the desired accuracy and time constraints.
In addition, there are countless ways to go about doing this,
be it Python libraries such as :mod:`rdkit` or :mod:`ase`, or
some sort of computational chemistry software. Since
:mod:`stk` cannot hope to provide a good solution to this problem,
it does try to make it easy for you to convert an
:mod:`stk` :class:`.Molecule` into whatever format you need to make
use of this other software. This means you can access atoms and
bonds with :meth:`.Molecule.get_atoms` and :meth:`.Molecule.get_bonds`,
you can convert any :mod:`stk` :class:`.Molecule` into an
:mod:`rdkit` molecule with :meth:`.Molecule.to_rdkit_mol` or you
can write it to a file with :meth:`.Molecule.write`.

.. figure:: https://i.imgur.com/UlCnTj9.png
    :align: center

    The general construction workflow of ``stk``.

The abstraction provided by the topology graph has a number of
powerful benefits. Firstly, because every vertex is responsible for the
placement of a building block, it is extremely easy to construct
different structural isomers of the constructed molecule. The vertex
can be told to perform different transformations on the building block,
so that its orientation in the constructed molecule changes. For the
end user, selecting the transformation and therefore structural isomer
is relatively easy. Take the example of an organic cage, which can be
constructed with the following code


.. code-block:: python

    # Create the building blocks.
    bb1 = stk.BuildingBlock('O=CC(C=O)C(Cl)C=O', ['aldehyde'])
    bb2 = stk.BuildingBlock('O=CC(C=O)C=O', ['aldehyde'])
    bb3 = stk.BuildingBlock('NCC(Cl)N', ['amine'])
    bb4 = stk.BuildingBlock('NCCN', ['amine'])

    # Create the topology graph.
    tetrahedron = stk.cage.FourPlusSix()

    # Because there are multiple building blocks with the same
    # number of functional groups, they need to be explicitly
    # placed on vertices, as there are multiple valid combinations.
    building_block_vertices = {
        bb1: tetrahedron.vertices[:1],
        bb2: tetrahedron.vertices[1:4],
        bb3: tetrahedron.vertices[4:5],
        bb4: tetrahedron.vertices[5:]
    }

    # Create the molecule.
    cage = stk.ConstructedMolecule(
        building_blocks=[bb1, bb2, bb3, bb4],
        topology_graph=tetrahedron,
        building_block_vertices=building_block_vertices
    )
    # You can write the molecule to a file if you want to view it.
    cage.write('cage.mol')

and looks like this

.. figure:: https://i.imgur.com/MAFrzAl.png


You can see that the green atoms on adjacent building blocks
point toward the different edges. However, by specifying a different
edge to align with, the building block will be rotated

.. code-block:: python

    # Vertex 0 gets aligned to the third edge it's connected to.
    isomer_graph = stk.cage.FourPlusSix(vertex_alignments={0: 2})
    building_block_vertices = {
        bb1: isomer_graph.vertices[:1],
        bb2: isomer_graph.vertices[1:4],
        bb3: isomer_graph.vertices[4:5],
        bb4: isomer_graph.vertices[5:]
    }
    isomer = stk.ConstructedMolecule(
        building_blocks=[bb1, bb2, bb3, bb4],
        topology_graph=tetrahedron,
        building_block_vertices=building_block_vertices
    )
    isomer.write('cage_isomer.mol')

.. figure:: https://i.imgur.com/cg9n69u.png


The same thing can be done to any other building block on the cage to
perform a rotation on it. You can also write a loop, to create all the
structural isomers of a single cage in one swoop

.. code-block:: python

    import itertools as it

    edges = (
        range(len(v.edges)) for v in stk.cage.FourPlusSix.vertex_data
    )
    # Create 5184 structural isomers.
    isomers = []
    for i, aligners in enumerate(it.product(*edges)):
        tetrahedron = stk.cage.FourPlusSix(
            vertex_alignments={
                vertex.id: edge
                for vertex, edge
                in zip(stk.cage.FourPlusSix.vertex_data, aligners)
            }
        )
        isomer = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2, bb3, bb4],
            topology_graph=tetrahedron,
            building_block_vertices={
                bb1: tetrahedron.vertices[:1],
                bb2: tetrahedron.vertices[1:4],
                bb3: tetrahedron.vertices[4:5],
                bb4: tetrahedron.vertices[5:]
            }
        )
        isomers.append(isomer)


The second major benefit of the topology graph is that the vertices and
edges can hold additional state useful for the construction of a
molecule. An example of this is in the construction of different
structural isomers, but another can be seen in the construction of
periodic systems. For example, :mod:`stk` allows you to construct
covalent organic frameworks. With the topology graph this is trivial
to implement, simply label some of the edges a periodic and they
will construct periodic bonds instead of regular ones.

The third benefit of the topology graph is that it allows users to
easily modify the construction of molecules by placing different
building blocks on different vertices. The user can use the
*building_block_vertices* parameter with any topology graph.

The fourth benefit of the topology graph is that the construction of
a molecule is broken down into a independent steps. Each vertex
represents a single, independent operation on a building block while
each edge represents a single, independent operation on a collection
of building blocks. As a result, each vertex and edge represents a
single operation, which can be executed in parallel. This allows
:mod:`stk` to scale efficiently to large topology graphs and take
advantage of multiple cores even during the construction of a single
molecule.
