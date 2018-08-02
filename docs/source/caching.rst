Caching
=======

One important feature of ``stk``, which is often the source of unexpected
behaviour, is that it caches every molecule it creates. In addition,
every time the the user attempts to build or load the same molecule, the
cached copy is returned from memory. Here are some examples

.. code-block:: python

    # Example 1
    bb1 = StructUnit('molecule1.mol')
    bb1 # <StructUnit at 0x0012931>
    bb2 = StructUnit('molecule1.mol')
    bb2 # <StructUnit at 0x0012931>
    bb1 is bb2 # True

    # Example 2 - when picking a different functional group a new
    # object will be created.
    bb3 = StructUnit('molecule1.mol', 'amine')
    bb3 # <StructUnit at 0x13131>
    bb3 is bb1 # False

    bb4 = StructUnit('molecule1.mol', 'aldehdye')
    bb4 # <StructUnit at 0x243134>
    bb3 is bb4 # False

    # Example 3 - if the molecule has the same structure, even if
    # the file is different, the cached object will be returned.
    bb5 = StructUnit('same_struct_as_bb1.mol2')
    bb6 = StructUnit('same_struct_as_bb1.mol2', 'amine')
    bb7 = StructUnit('same_struct_as_bb1.mol2', 'aldehdye')

    bb5 is bb1 and bb5 is bb2 # True
    bb6 is bb3 # True
    bb7 is bb4 # True

    # Example 4 - each class has its own cache.
    bb8 = StructUnit2('molecule1.mol')
    bb9 = StructUnit3('molecule1.mol')

    bb8 is bb1 # False
    bb8 is bb9 # False

and with :class:`.MacroMolecule`

.. code-block:: python

    # Example 5 - using the same objects for initialization results
    # in the same object being returned from the cache. Assembly step is
    # skipped!
    mmol1 = Cage([bb1, bb2], FourPlusSix()) # Slow, have to build.
    mmol2 = Cage([bb1, bb2], FourPlusSix()) # Fast, returned from memory.

    mmol1 is mmol2 # True

    # Changing the topology will mean a new object is built.
    mmol3 = Cage([bb1, bb2], FourPlusSix(A_alignments=[1, 0, 2, 1]))
    mmol1 is mmol3 # False


The user has the option to turn the cache off and on through the
``CACHE_SETTINGS`` variable.

.. code-block:: python

    CACHE_SETTINGS['ON'] = False # Cache is off.
    bb10 = StructUnit('molecele1.mol')
    bb10 is bb1 # False

    CACHE_SETTINGS['ON'] = True # Cache is back on.
    bb11 = StructUnit('molecule1.mol')
    bb11 is bb1 # True
    bb11 is bb10 # False
