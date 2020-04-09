==============
Basic Examples
==============


Creating Building Blocks
========================

There are a number of ways to create :class:`.BuildingBlock` molecules,
but the two most common ways are direct creation through SMILES
strings and loading them from molecule structure files.


.. code-block:: python

    import stk

    bb1 = stk.BuildingBlock('NCCCN')
    bb2 = stk.BuildingBlock.init_from_file('path/to/file.mol')


Look at the documentation of :class:`.BuildingBlock` to see other
available initialization methods. They are easy to find, as they all
begin with :mod:`init`.

Specifying Multiple Functional Groups
=====================================

When you create a :class:`.BuildingBlock`, you need also specify
which atoms are modified during construction of a
:class:`.ConstructedMolecule`. This is achieved by providing the
:class:`.BuildingBlock` with
:mod:`~.functional_groups.functional_group` instances. To save you
the pain of creating function groups one by one, you can use a
:mod:`~.functional_group_factory`. If you have a building block
with bromo groups, and you want the bromo groups to be modified
during construction, you would use a :class:`.BromoFactory`

.. code-block:: python

    import stk

    bb = stk.BuildingBlock('BrCCCBr', [stk.BromoFactory()])

The ``bb``, in the example above, would have two :class:`.Bromo`
functional groups. Similarly, if you have a building block with
aldehyde groups

.. code-block:: python

    bb2 = stk.BuildingBlock('O=CCCC=O', [stk.AldehydeFactory()])

In this example, ``bb2`` will have two :class:`.Aldehyde` functional
groups. Finally, if you have both aldehyde and bromo groups on a
molecule, and you want both to be modified during construction,
you would use both of the factories

.. code-block:: python

    bb3 = stk.BuildingBlock(
        smiles='O=CCCBr',
        functional_groups=[stk.AldehydeFactory(), stk.BromoFactory()],
    )

In the example above, ``bb3`` has one :class:`.Bromo` and one
:class:`.Aldehyde` functional group.

Constructing Molecules
======================

To construct molecules, you need to create a new
:class:`.ConstructedMolecule`. The required input consists of
a :class:`.TopologyGraph`, which, in turn,  requires
:class:`.BuildingBlock` instances.

.. code-block:: python

    import stk

    # React the amine functional groups during construction.
    bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
    # React the aldehyde functional groups during construction.
    bb2 = stk.BuildingBlock('O=CC(C=O)C=O', [stk.AldehydeFactory()])
    # Build a polymer.
    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            building_blocks=(bb1, bb2),
            repeating_unit='AB',
            num_repeating_units=12,
        ),
    )

    # Build a longer polymer.
    longer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            building_blocks=(bb1, bb2),
            repeating_unit='AB',
            num_repeating_units=23,
        ),
    )

Each topology graph requires different input parameters,
for example, organic cage topology graphs only require the
:class:`.BuildingBlock` instances.

.. code-block:: python

    cage = stk.ConstructedMolecule(stk.cage.FourPlusSix((bb1, bb2)))


Read the documentation for each kind of :class:`.TopologyGraph`, for
more examples on how to initialize it, and to see what optional
parameters you have available.

Using RDKit to Optimize Molecular Structures
============================================

Molecules used by :mod:`stk` can be structure optimized both before and
after construction. One easy way to do is, is with the
:mod:`rdkit` library. You can optimize any :mod:`stk`
:class:`.Molecule`, such as a :class:`.BuildingBlock`

.. code-block:: python

    import stk
    import rdkit.Chem.AllChem as rdkit

    bb = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])

    # Optimize with the MMFF force field.

    rdkit_bb = bb.to_rdkit_mol()
    rdkit.SanitizeMol(rdkit_bb)
    rdkit.MMFFOptimizeMolecule(rdkit_bb)

    # stk molecules are immutable. with_position_matrix returns a
    # a clone, holding the new position matrix.
    bb = bb.with_position_matrix(
        position_matrix=rdkit_bb.GetConformer().GetPositions(),
    )

or a :class:`.ConstructedMolecule`

.. code-block:: python

    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear((bb, ), 'A', 15),
    )

    # Optimize with the MMFF force field.

    rdkit_polymer = polymer.to_rdkit_mol()
    rdkit.SanitizeMol(rdkit_polymer)
    rdkit.MMFFOptimizeMolecule(rdkit_polymer)

    # stk molecules are immutable. with_position_matrix returns a
    # a clone, holding the new position matrix.
    polymer = polymer.with_position_matrix(
        position_matrix=rdkit_polymer.GetConformer().GetPositions(),
    )

Writing Molecular Files
=======================

The simplest way to save molecules is to write them to a file.
This works with any :class:`.Molecule`, including both the
:class:`.BuildingBlock`

.. code-block:: python

    import stk

    bb = stk.BuildingBlock(
        smiles='ICCBr',
        functional_groups=[stk.BromoFactory()],
    )
    bb.write('bb.mol')


and the :class:`.ConstructedMolecule`

.. code-block:: python

    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear((bb, ), 'A', 10),
    )
    polymer.write('polymer.mol')

You can see what file formats are supported by reading the
documentation for :meth:`~.Molecule.write`.

Placing and Retrieving Molecules From a Database
================================================

Requirements
------------

:mod:`stk` allows you to place molecules into a
:class:`.MoleculeDatabase`. Out-of-the-box, :mod:`stk` comes
with support for a :class:`.MoleculeMongoDb`. In order to use it
locally, you have to install MongoDB on your computer. You will
then also have to install :mod:`pymongo` with::

    $ pip install pymongo

Documentation for making sure your local MongoDB is working properly
can be found here__.

__ https://api.mongodb.com/python/current/

You can also use a remote MongoDB, in which case you do not have to
install it locally, but you will still need to install
:mod:`pymongo`.

Molecules and Building Blocks
-----------------------------

To place molecules into the database, first create the database

.. code-block:: python

    import stk
    import pymongo


    # Connect to a MongoDB. This example connects to a local
    # MongoDB, but you can connect to a remote DB too with
    # MongoClient() - read the documentation for pymongo to see how
    # to do that.
    client = pymongo.MongoClient()
    db = stk.MoleculeMongoDb(client)

You then create and place a molecule into the database,
for example, a :class:`.BuildingBlock`

.. code-block:: python

    bb = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])
    db.put(bb)


To restore it, by default, you would use the InChIKey

.. code-block:: python

    import rdkit.Chem.AllChem as rdkit

    restored = db.get({
        'InChIKey': rdkit.MolToInchiKey(rdkit.MolFromSmiles('BrCCBr')),
    })

However, you can customize it. For example, the documentation of
:class:`.MoleculeMongoDb`, shows how you can use SMILES to retrieve
your molecules.

The ``restored`` molecule is only a :class:`.Molecule` instance,
and not a :class:`.BuildingBlock` instance, which means that it lacks
functional groups. You can restore your functional groups however

.. code-block:: python

    restored_bb = stk.BuildingBlock.init_from_molecule(
        molecule=restored,
        functional_groups=[stk.BromoFactory()],
    )

Constructed Molecules
---------------------

You can use the same database for placing
:class:`.ConstructedMolecule` instances

.. code-block:: python

    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear((bb, ), 'A', 2),
    )
    db.put(polymer)

and restore them in the same way

.. code-block:: python

    restored = db.get({
        'InChIKey': rdkit.MolToInchiKey(rdkit.MolFromSmiles(
            'BrCCCCBr'
        )),
    })

However, once again, ``restored`` will only be a :class:`.Molecule`
instance, and not a :class:`.ConstructedMolecule` instance.

If you want to store and retrieve :class:`.ConstructedMolecule`
instances, you have to create a :class:`.ConstructedMoleculeMongoDb`

.. code-block:: python

    constructed_db = stk.ConstructedMoleculeMongoDb(client)
    constructed_db.put(polymer)
    restored_polymer = constructed_db.get({
        'InChIKey': rdkit.MolToInchiKey(rdkit.MolFromSmiles(
            'BrCCCCBr'
        )),
    })

Unlike ``restored``, ``restored_polymer`` is a
:class:`.ConstructedMolecule` instance.

Placing and Retrieving Molecular Property Values From a Database
================================================================

Requirements
------------

Using a :class:`.ValueMongoDb` has the same requirements as the
previous example.

Storing Values
--------------

Unlike the previous example, you can deposit values for both
a :class:`.BuildingBlock` and a :class:`.ConstructedMolecule` in the
same database. First lets create one

.. code-block:: python

    import stk

    # Connect to a MongoDB. This example connects to a local
    # MongoDB, but you can connect to a remote DB too with
    # MongoClient() - read the documentation for pymongo to see how
    # to do that.
    client = pymongo.MongoClient()

    # You have to choose name for your collection.
    energy_db = stk.ValueMongoDb(client, 'energy')

Here, ``energy_db`` will store energy values. Lets create a function
to calculate the energy of a molecule.

.. code-block:: python

    import rdkit.Chem.AllChem as rdkit

    def get_energy(molecule):
        rdkit_molecule = molecule.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_molecule)
        ff = rdkit.UFFGetMoleculeForceField(rdkit_molecule)
        return ff.CalcEnergy()

Now we can deposit the energy value in the database

.. code-block:: python

    bb = stk.BuildingBlock('BrCCCCBr')
    energy_db.put(bb, get_energy(bb))

And we can retrieve it

.. code-block:: python

    energy = energy_db.get(bb)


If we make the same molecule in some other way, for example we
can make ``BrCCCCBr`` as a constructed molecule

.. code-block:: python

    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            building_blocks=(
                stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
            ),
            repeating_unit='A',
            num_repeating_units=2,
        ),
    )

we can still retrieve the value

.. code-block:: python

    # You get the correct energy out, because polymer and bb are
    # actually the same molecule.
    bb_energy = energy_db.get(polymer)


You can also you a :class:`.ConstructedMolecule` to deposit values
into the database, for example

.. code-block:: python

    atom_count_db = stk.ValueMongoDb(client, 'atom_counts')
    atom_count_db.put(polymer, polymer.get_num_atoms())


These values will also be accessible in a later session

.. code-block:: python

    # Assume this a new Python session.
    import stk
    import pymongo


    client = pymongo.MongoClient()
    energy_db = stk.ValueMongoDb(client, 'energy')
    atom_count_db = stk.ValueMongoDb(client, 'atom_counts')

    bb = stk.BuildingBlock('BrCCCCBr')
    bb_energy = energy_db.get(bb)
    bb_atom_count = energy.get(bb)

Specifying Functional Groups Individually
=========================================

If you want to be more precise about which functional groups get
created, you can provide them directly to the :class:`.BuildingBlock`.
For example, if you have multiple bromo groups on a molecule, but
you only want to use one during construction

.. code-block:: python

    import stk

    bb = stk.BuildingBlock(
        smiles='BrCCCBr',
        functional_groups=[
            stk.Bromo(
                # The number is the atom's id.
                bromine=stk.Br(0),
                atom=stk.C(1),
                # bonders are atoms which have bonds added during
                # construction.
                bonders=(stk.C(1), ),
                # deleters are atoms which are deleted during
                # construction.
                deleters=(stk.Br(0), ),
            ),
        ],
    )

When creating a :class:`.Bromo` functional group, you have to
specify things like which atoms have bonds added during construction,
and which ones are removed during construction. These are specified by
the `bonders` and `deleters` parameters, respectively. You can add
as many functional groups to :class:`.BuildingBlock` as you like
in this way, and you can mix different types of
:mod:`~.functional_groups.functional_group`. You can even mix
a :mod:`~.functional_groups.functional_group` instances with
:mod:`~.functional_group_factory` instances.

Changing Bonder and Deleter Atoms in Functional Group Factories
===============================================================

In the previous example, you saw that during creation of a
:class:`.Bromo` instance, you can specify which atoms have bonds
added during construction, and which atoms are deleted during
construction. You might like to customize this in the functional groups
created by a :mod:`~.functional_group_factory`.

Take, for example, a :class:`.CarboxylicAcid` functional group. There
are two likely ways you would like to modify
this group, ``C(=O)O``, during construction. In the first way, you want
to add a bond to the carbon atom, and delete the ``OH`` group, which is
treated as a leaving group. This is what
:class:`.CarboxylicAcidFactory` will do by default

.. code-block:: python

    import stk

    bb = stk.BuildingBlock(
        smiles='O=C(O)CCC(=O)O',
        functional_groups=[stk.CarboxylicAcidFactory()],
    )

Here, ``bb`` will have two :class:`.CarboxylicAcid` functional groups.
In each, the deleter atoms will be the oxygen and hydrogen atom of
the ``OH`` group, and the bonder atom will be the carbon atom.

Now, the second way you might want to modify a carobxylic acid group,
is to only delete the hydrogen atom of the ``OH`` group during
construction, and add a bond to the oxygen atom of the
``OH`` group. This means the hydrogen atom is the deleter atom and
the oxygen atom is the bonder atom. You can tell the
:class:`.CarboxylicAcidFactory` to create :class:`.CarboxylicAcid`
instances of this kind

.. code-block:: python

    bb2 = stk.BuildingBlock(
        smiles='O=C(O)CCC(=O)O',
        functional_groups=[
            stk.CarboxylicAcidFactory(
                # Atom number 3 corresponds to the OH oxygen atom in a
                # carboxylic acid group. THIS IS NOT THE ATOM'S ID IN
                # THE MOLECULE.
                bonders=(3, ),
                # Atom number 4 corresponds to the hydrogen atom in a
                # carboxylic acid group. THIS IS NOT THE ATOM'S ID IN
                # THE MOLECULE.
                deleters=(4, ),
            ),
        ],
    )

Here, ``bb2`` will also have two :class:`.CarboxylicAcid` functional
groups. In each, the deleter atom will be the hydrogen of the
``OH`` group and the bonder atom will be the oxygen atom of the
``OH`` group.

You might be wondering: "How do I know which number to use for
which atom in the functional group, so that I can specify the correct
atoms to be the bonders or deleters?" The docstring of
:class:`.CarboxylicAcidFactory` will tell you which number corresponds
to which atom in the functional group. The same is true for any
other :mod:`~.functional_group_factory`. Note that the number you
provide to the factory, is not the id of the atom found in the
molecule!!

Extending stk
=============

There are a lot of way to extend :mod:`stk`, for example by adding
new functional groups, topology graphs, mutation operations and so on.
However, because every part of :mod:`stk` is built around an abstract
base class, all you need to do is find the appropriate abstract base
class, and create a new subclass for it. The abstract base class
will provide documentation and examples on how to create a subclass.
You can easily find the abstract base classes by looking at the
sidebar.
