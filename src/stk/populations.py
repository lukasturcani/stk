"""
Populations
===========

:class:`.Population` objects are container specialized to perform
operations on groups of ``stk`` molecules, often in parallel.

"""

import itertools as it
import os
from os.path import join
import numpy as np
import json
import psutil
from functools import wraps
import logging
import pathos

from .utilities import dedupe, dice_similarity
from .molecular import ConstructedMolecule, Molecule


logger = logging.getLogger(__name__)


class Population:
    """
    A container for  :class:`.Molecule` objects.

    :class:`Population` instances can be nested.

    In addition to holding :class:`.Molecule` objects, the
    :class:`Population` class can be used to create large numbers of
    these instances through the class methods beginning with "init".

    :class:`.Molecule` instances held by a :class:`Population` can have
    their structures optimized in parallel through the
    :meth:`optimize` method.

    It supports all expected and necessary container operations such as
    iteration, indexing and membership checks (via the ``is in``
    operator).

    Attributes
    ----------
    direct_members : :class:`list` of :class:`.Molecule`
        Held here are direct members of the :class:`Population`.
        In other words, these are the molecules not held by any
        subpopulations. As a result, not all members of a
        :class:`Population` are stored in this attribute.

    subpopulations : :class:`list` of :class:`Population`
        A :class:`list` holding the subpopulations.

    Examples
    --------
    A :class:`Population` can be iterated through just like a
    :class:`list`

    .. code-block:: python

        import stk

        # Create a population.
        pop = stk.Population(
            stk.BuildingBlock(...),
            stk.ConstructedMolecule(...),
            stk.BuildingBlock(...),
            stk.BuildingBlock(...),

            stk.Population(
                stk.BuildingBlock(...),
                stk.ConstructedMolecule(...)
            )

            stk.ConstructedMolecule(...),
            stk.BuildingBlock(...)

        )

        for member in pop:
            do_stuff(member)

    When iterating through a :class:`Population` you will also iterate
    through nested members, that is members which are held by
    subpopulations. If you only wish to iterate through direct
    members, you can

    .. code-block:: python

        for member in pop.direct_members:
            do_stuff(member)

    You can also get access to members by using indices.
    Indices have access to all members in the population

    .. code-block:: python

        first_member = pop[0]
        second_member = pop[1]

    Indices will first access direct members of the population and then
    access members in the subpopulations. Indices access nested members
    depth-first

    .. code-block:: python

        pop2 = stk.Population(bb1, bb2, stk.Population(bb3, bb4))
        # Get bb1.
        pop2[0]
        # Get bb2.
        pop2[1]
        # Get bb3.
        pop2[2]
        # Get bb4.
        pop2[3]

    You can get a subpopulation by taking a slice

    .. code-block:: python

        # new_pop is a new Population instance and has no nesting.
        new_pop = pop[2:4]

    You can take the length of a population to get the total number
    of members

    .. code-block:: python

        len(pop)

    Adding populations creates a new population with both of the added
    populations as subpopulations

    .. code-block:: python

        # added has no direct members and two subpopulations, pop and
        # pop2.
        added = pop + pop2

    Subtracting populations creates a new, flat population.

    .. code-block:: python

        # subbed has all objects in pop except those also found in
        # pop2.
        subbed = pop - pop2

    You can check if an object is already present in the population.

    .. code-block:: python

        bb1 = stk.BuildingBlock(...)
        bb2 = stk.BuildingBlock(...)
        pop3 = stk.Population(bb1)

        # Returns True.
        bb1 in pop3
        # Returns False.
        bb2 in pop3
        # Returns True.
        bb2 not in pop3

    If you want to run multiple :meth:`optimize` calls in a row, use
    the "with" statement. This keeps a single process pool open, and
    means you do not create a new one for each :meth:`optimize` call.
    It also automatically closes the pool for you when the block
    exits

    .. code-block:: python

        population = stk.Population(...)
        # Keep a process pool open through the "with" statement.
        with population.open_process_pool(8):
            # All optimize calls within this block will use the
            # same process pool.
            population.optimize(stk.UFF())
            population.add_members(...)
            population.optimize(stk.UFF())
        # Process pool is automatically cleaned up when the block
        # exits.

    """

    def __init__(self, *args):
        """
        Initialize a :class:`Population`.

        Parameters
        ----------
        *args : :class:`.Molecule`, :class:`Population`
            A population is initialized with the :class:`.Molecule` and
            :class:`Population` instances it should hold.

        Examples
        --------
        .. code-block:: python

            bb1 = stk.BuildingBlock('CCC')
            bb2 = stk.BuildingBlock('NCCNCNC')
            bb3 = stk.BuildingBlock('[Br]CCC[Br]')
            pop1 = stk.Population(bb1, bb2, bb3)

            bb4 = stk.BuildingBlock('NNCCCN')
            # pop2 has pop1 as a subpopulation and bb4 as a direct
            # member.
            pop2 = stk.Population(pop1, bb4)

        """

        self.direct_members = []
        self.subpopulations = []
        self._process_pool = None

        for arg in args:
            if isinstance(arg, Population):
                self.subpopulations.append(arg)
            else:
                self.direct_members.append(arg)

    @classmethod
    def init_all(
        cls,
        building_blocks,
        topology_graphs,
        num_processes=None,
        duplicates=False,
        use_cache=False
    ):
        """
        Make all possible molecules from groups of building blocks.

        Parameters
        ----------
        building_blocks : :class:`list` of :class:`.Molecule`
            A :class:`list` holding nested building blocks, for
            example

            .. code-block:: python

                bbs1 = [
                    stk.BuildingBlock(...),
                    stk.BuildingBlock(...),
                    ...
                ]
                bbs2 = [
                    stk.ConstructedMolecule(...),
                    stk.BuildingBlock(...),
                    ...,
                ]
                bbs3 = [
                    stk.BuildingBlock(...),
                    stk.BuildingBlock(...),
                    ...
                ]
                building_blocks = [bbs1, bbs2, bbs3]

            To construct a new :class:`.ConstructedMolecule`, a
            :class:`.Molecule` is picked from each of the sublists
            in `building_blocks`. The picked :class:`.Molecule`
            instances are then supplied to
            :class:`.ConstructedMolecule`

            .. code-block:: python

                # mol is a new ConstructedMolecule. bb1 is selected
                # from bbs1, bb2 is selected from bbs2 and bb3 is
                # selected from bbs3.
                mol = stk.ConstructedMolecule(
                    building_blocks=[bb1, bb2, bb3],
                    topology_graph=topology_pick
                )

            The order a :class:`.Molecule` instance is given to
            the :class:`.ConstructedMolecule` is determined by the
            sublist of `building_blocks` it was picked from. Note that
            the number of sublists in `building_blocks` is not fixed.
            It merely has to be compatible with the
            `topology_graphs`.

        topology_graphs : :class:`list` of :class:`.TopologyGraph`
            The topology graphs of `.ConstructedMolecule` being made.

        num_processes : :class:`int`, optional
            The number of parallel processes to create when
            constructing the molecules. If ``None``, creates a process
            for each core on the computer.

        duplicates : :class:`bool`, optional
            If ``False``, duplicate structures are removed from
            the population.

        use_cache : :class:`bool`, optional
            Toggles use of the molecular cache.

        Returns
        -------
        :class:`Population`
            A :class:`.Population` holding `.ConstructedMolecule`
            instances.

        Examples
        --------
        Construct all possible cage molecules from some precursors

        .. code-block:: python

            import stk

            amines = [
                stk.BuildingBlock('NCCCN', ['amine']),
                stk.BuildingBlock('NCCCCCN', ['amine']),
                stk.BuildingBlock('NCCOCCN', ['amine']),
            ]
            aldehydes = [
                stk.BuildingBlock('O=CCC(C=O)CC=O', ['aldehyde']),
                stk.BuildingBlock('O=CCC(C=O)CC=O', ['aldehyde']),
                stk.BuildingBlock('O=C(C=O)COCC=O', ['aldehyde']),
            ]
            # A total of 9 cages will be created.
            cages = stk.Population.init_all(
                building_blocks=[amines, aldehydes],
                topology_graphs=[stk.cage.FourPlusSix()]
            )

        Use the constructed cages and a new bunch of building blocks to
        create all possible cage complexes.

        .. code-block:: python

            encapsulants = [
                stk.BuildingBlock('[Br][Br]'),
                stk.BuildingBlock('[F][F]'),
            ]

            # Every combination of cage and encapsulant.
            complexes = stk.Population.init_all(
                building_blocks=[cages, encapsulants],
                topology_graphs=[stk.host_guest_complex.Complex()]
            )

        """

        bbs, topologies = [], []
        mols = it.product(*building_blocks, topology_graphs)
        for *mol_bbs, topology in mols:
            bbs.append(mol_bbs),
            topologies.append(topology)

        with pathos.pools.ProcessPool(num_processes) as pool:
            mols = pool.map(ConstructedMolecule, bbs, topologies)

        # Update the cache.
        if use_cache:
            for i, mol in enumerate(mols):
                # If the molecule did not exist already, add it to the
                # cache.
                if (
                    not ConstructedMolecule.has_cached_mol(
                        identity_key=mol.get_identity_key()
                    )
                ):
                    mol.update_cache()
                # If the molecule did exist already, use the cached
                # version.
                else:
                    mols[i] = ConstructedMolecule.get_cached_mol(
                        identity_key=mol.get_identity_key()
                    )

        p = cls(*mols)
        if not duplicates:
            p.remove_duplicates()
        return p

    def clone(self):
        """
        Return a clone.

        The clone will share the :class:`.Molecule` objects, copies of
        :class:`.Molecule` objects will not be made.

        Returns
        -------
        :class:`Population`
            The clone.

        Examples
        --------
        .. code-block:: python

            import stk

            # Make an intial population.
            pop = stk.Population(stk.BuildingBlock('NCCN'))
            # Make a clone.
            clone = pop.clone()

        """

        clone = self.__class__(*self.direct_members)
        clone.subpopulations = [
            pop.clone() for pop in self.subpopulations
        ]
        return clone

    @classmethod
    def init_diverse(
        cls,
        building_blocks,
        topology_graphs,
        size,
        random_seed=None,
        use_cache=False
    ):
        """
        Construct a chemically diverse :class:`.Population`.

        All constructed molecules are held in :attr:`direct_members`.

        In order to construct a :class:`.ConstructedMolecule`, a random
        :class:`.Molecule` is selected from each sublist in
        `building_blocks`. Once the first construction is complete,
        the next :class:`.Molecule` selected from each sublist is the
        one with the most different Morgan fingerprint to the prior
        one. The third construction uses randomly selected
        :class:`.Molecule` objects again and so on. This is done until
        `size` :class:`.ConstructedMolecule` instances have been
        constructed.

        Parameters
        ----------
        building_blocks : :class:`list` of :class:`.Molecule`
            A :class:`list` holding nested building blocks, for
            example

            .. code-block:: python

                bbs1 = [
                    stk.BuildingBlock(...),
                    stk.BuildingBlock(...),
                    ...
                ]
                bbs2 = [
                    stk.ConstructedMolecule(...),
                    stk.BuildingBlock(...),
                    ...
                ]
                bbs3 = [
                    stk.BuildingBlock(...),
                    stk.BuildingBlock(...),
                    ...
                ]
                building_blocks = [bbs1, bbs2, bbs3]

            To construct a new :class:`.ConstructedMolecule`, a
            :class:`.Molecule` is picked from each of the sublists
            in `building_blocks`. The picked :class:`.Molecule`
            instances are then supplied to the
            :class:`.ConstructedMolecule`

            .. code-block:: python

                # mol is a new ConstructedMolecule. bb1 is selected
                # from bbs1, bb2 is selected from bbs2 and bb3 is
                # selected from bbs3.
                mol = stk.ConstructedMolecule(
                    building_blocks=[bb1, bb2, bb3],
                    topology_graph=topology_pick
                )

            The order a :class:`.Molecule` instance is given to
            :class:`.ConstructedMolecule` is determined by the
            sublist of `building_blocks` it was picked from. Note that
            the number of sublists in `building_blocks` is not fixed.
            It merely has to be compatible with the `topology_graphs`.

        topology_graphs : :class:`iterable` of :class:`.TopologyGraph`
            An iterable holding topology grpahs which should be
            randomly selected for the construction of a
            :class:`.ConstructedMolecule`.

        size : :class:`int`
            The desired size of the :class:`.Population`.

        random_seed : :class:`int`, optional
            Seed for the random number generator to get replicable
            results.

        use_cache : :class:`bool`, optional
            Toggles use of the molecular cache.

        Returns
        -------
        :class:`.Population`
            A population filled with the constructed molecules.

        Examples
        --------
        Construct a diverse :class:`.Population` of cage
        molecules from some precursors

        .. code-block:: python

            import stk

            amines = [
                stk.BuildingBlock('NCCCN', ['amine']),
                stk.BuildingBlock('NCCCCCN', ['amine']),
                stk.BuildingBlock('NCCOCCN', ['amine']),
            ]
            aldehydes = [
                stk.BuildingBlock('O=CCC(C=O)CC=O', ['aldehyde']),
                stk.BuildingBlock('O=CCC(C=O)CC=O', ['aldehyde']),
                stk.BuildingBlock('O=C(C=O)COCC=O', ['aldehyde']),
            ]
            # A total of 4 cages will be created.
            cages = stk.Population.init_diverse(
                building_blocks=[amines, aldehydes],
                topology_graphs=[stk.cage.FourPlusSix()],
                size=4
            )

        Use the constructed cages and a new bunch of building blocks to
        create some diverse cage complexes.

        .. code-block:: python

            encapsulants = [
                stk.BuildingBlock('[Br][Br]'),
                stk.BuildingBlock('[F][F]'),
            ]

            # 4 combinations of cage and encapsulant.
            complexes = stk.Population.init_diverse(
                building_blocks=[cages, encapsulants],
                topology_graphs=[stk.host_guest_complex.Complex()],
                size=4
            )

        """

        pop = cls()

        # Shuffle the sublists.
        generator = np.random.RandomState(random_seed)
        for db in building_blocks:
            generator.shuffle(db)

        # Go through every possible constructed molecule.
        for *bbs, top in it.product(*building_blocks, topology_graphs):

            # Generate the random constructed molecule.
            mol = ConstructedMolecule(
                building_blocks=bbs,
                topology_graph=top,
                use_cache=use_cache
            )
            if mol not in pop:
                pop.direct_members.append(mol)

            if len(pop) == size:
                break

            # Get the most different BuildingBlock to the previously
            # selected one, per sublist.
            diff_bbs = [
                min(db, key=lambda mol: dice_similarity(bb, mol))
                for bb, db in zip(bbs, building_blocks)

            ]
            mol = ConstructedMolecule(
                building_blocks=diff_bbs,
                topology_graph=top,
                use_cache=use_cache
            )
            if mol not in pop:
                pop.direct_members.append(mol)

            if len(pop) == size:
                break

        assert len(pop) == size
        return pop

    @classmethod
    def init_random(
        cls,
        building_blocks,
        topology_graphs,
        size,
        random_seed=None,
        use_cache=False
    ):
        """
        Construct molecules for a random :class:`.Population`.

        All molecules are held in :attr:`direct_members`.

        From the supplied building blocks a random :class:`.Molecule`
        is selected from each sublist to form a
        :class:`.ConstructedMolecule`. This is done until `size`
        :class:`.ConstructedMolecule` objects have been constructed.

        Parameters
        ----------
        building_blocks : :class:`list` of :class:`.Molecule`
            A :class:`list` holding nested building blocks, for
            example

            .. code-block:: python

                bbs1 = [
                    stk.BuildingBlock(...),
                    sk.BuildingBlock(...),
                    ...
                ]
                bbs2 = [
                    stk.ConstructedMolecule(...),
                    stk.BuildingBlock(...),
                    ...
                ]
                bbs3 = [
                    stk.BuildingBlock(...),
                    stk.BuildingBlock(...),
                    ...
                ]
                building_blocks = [bbs1, bbs2, bbs3]

            To construct a new :class:`.ConstructedMolecule`, a
            :class:`.Molecule` is picked from each of the sublists
            in `building_blocks`. The picked :class:`.Molecule`
            instances are then supplied to
            :class:`.ConstructedMolecule`

            .. code-block:: python

                # mol is a new ConstructedMolecule. bb1 is selected
                # from bbs1, bb2 is selected from bbs2 and bb3 is
                # selected from bbs3.
                mol = stk.ConstructedMolecule(
                    building_blocks=[bb1, bb2, bb3],
                    topology_graph=topology_pick
                )

            The order a :class:`.Molecule` instance is given to
            the :class:`.ConstructedMolecule` is determined by the
            sublist of `building_blocks` it was picked from. Note that
            the number of sublists in `building_blocks` is not fixed.
            It merely has to be compatible with the `topology_graphs`.

        topology_graphs : :class:`iterable` of :class:`.TopologyGraph`
            An :class:`iterable` holding topology graphs which should
            be randomly selected during initialization of
            :class:`.ConstructedMolecule`.

        size : :class:`int`
            The size of the population to be initialized.

        random_seed : :class:`int`, optional
            Seed for the random number generator to get replicable
            results.

        use_cache : :class:`bool`, optional
            Toggles use of the molecular cache.

        Returns
        -------
        :class:`.Population`
            A population filled with random
            :class:`.ConstructedMolecule` instances.

        Examples
        --------
        Construct 5 random cage molecules from some
        precursors

        .. code-block:: python

            import stk

            amines = [
                stk.BuildingBlock('NCCCN', ['amine']),
                stk.BuildingBlock('NCCCCCN', ['amine']),
                stk.BuildingBlock('NCCOCCN', ['amine']),
            ]
            aldehydes = [
                stk.BuildingBlock('O=CCC(C=O)CC=O', ['aldehyde']),
                stk.BuildingBlock('O=CCC(C=O)CC=O', ['aldehyde']),
                stk.BuildingBlock('O=C(C=O)COCC=O', ['aldehyde']),
            ]
            # A total of 5 cages will be created.
            cages = stk.Population.init_random(
                building_blocks=[amines, aldehydes],
                topology_graphs=[stk.cage.FourPlusSix()],
                size=5
            )

        Use the constructed cages and a new bunch of building blocks to
        create some random cage complexes.

        .. code-block:: python

            encapsulants = [
                stk.BuildingBlock('[Br][Br]'),
                stk.BuildingBlock('[F][F]'),
            ]

            # Random combinations of cage and encapsulant.
            complexes = stk.Population.init_random(
                building_blocks=[cages, encapsulants],
                topology_graphs=[stk.host_guest_complex.Complex()],
                size=5
            )

        """

        pop = cls()
        generator = np.random.RandomState(random_seed)

        # Shuffle the sublists.
        for db in building_blocks:
            generator.shuffle(db)

        # Go through every possible constructed molecule.
        for *bbs, top in it.product(*building_blocks, topology_graphs):
            # Generate the random constructed molecule.
            mol = ConstructedMolecule(
                building_blocks=bbs,
                topology_graph=top,
                use_cache=use_cache
            )
            if mol not in pop:
                pop.direct_members.append(mol)

            if len(pop) == size:
                break

        assert len(pop) == size
        return pop

    @classmethod
    def init_from_list(cls, pop_list, use_cache=False):
        """
        Initialize a population from a :class:`list` representation.

        Parameters
        ----------
        pop_list : :class:`list`
            A :class:`list` which represents a :class:`Population`.
            Like the ones created by :meth:`to_list`. For example in,

            .. code-block:: python

                pop_list = [{...}, [{...}], [{...}, {...}], {...}]

            ``pop_list`` represents the :class:`Population`, sublists
            represent its subpopulations and the :class:`dict`
            ``{...}`` represents the members.

        use_cache : :class:`bool`, optional
            Toggles use of the molecular cache.

        Returns
        -------
        :class:`Population`
            The population represented by `pop_list`.

        """

        pop = cls()
        for item in pop_list:
            if isinstance(item, list):
                sp = cls.init_from_list(item, use_cache=use_cache)
                pop.subpopulations.append(sp)
            else:
                pop.direct_members.append(
                    Molecule.init_from_dict(item, use_cache=use_cache)
                )
        return pop

    def add_members(self, molecules, duplicate_key=None):
        """
        Add :class:`.Molecule` instances to the :class:`.Population`.

        The added :class:`.Molecule` instances are added as direct
        members of the population, they are not placed into any
        subpopulations.

        Parameters
        ----------
        molecules : :class:`iterable` of :class:`.Molecule`
            The molecules to be added as direct members.

        duplicate_key : :class:`callable`, optional
            If not ``None``, ``duplicate_key(mol)`` is evalued on each
            molecule in `members`. If a molecule with the same
            `duplicate_key` is already present in the population, the
            molecule is not added.

        Returns
        -------
        None : :class:`NoneType`

        """

        if duplicate_key is None:
            self.direct_members.extend(molecules)
        else:
            keys = {duplicate_key(mol) for mol in self}
            for mol in molecules:
                key = duplicate_key(mol)
                if key not in keys:
                    keys.add(key)
                    self.direct_members.append(mol)

    def add_subpopulation(self, population):
        """
        Add a clone of `population` to :attr:`subpopulations`.

        Only a clone of the `population` container is made. The
        molecules it holds are not copies.

        Parameters
        ----------
        population : :class:`Population`
            The population to be added as a subpopulation.

        Returns
        -------
        None : :class:`NoneType`

        """

        self.subpopulations.append(population.clone())

    def _get_all_members(self):
        """
        Yield both direct and nested members of the population.

        Yields
        ------
        :class:`.Molecule`
            The next :class:`.Molecule` instance held within the
            :class:`.Population`.

        """

        # Go through `direct_members` attribute and yield ``Molecule``
        # instances held within one by one.
        for ind in self.direct_members:
            yield ind

        # Go thorugh `populations` attribute and for each
        # ``Population`` instance within, yield ``Molecule``
        # instances from its `all_members` generator.
        for pop in self.subpopulations:
            yield from pop._get_all_members()

    def set_mol_ids(self, n, overwrite=False):
        """
        Give each member of the population an id starting from `n`.

        This method adds an :attr:`id` attribute to each
        :class:`.Molecule` instance held by the population.

        Parameters
        ----------
        n : :class:`int`
            A number. Members of this :class:`Population` are given a
            unique number as an id, starting from `n` and
            incremented by one between members.

        overwrite : :class:`bool`, optional
            If ``True``, existing ids are replaced.

        Returns
        -------
        :class:`int`
            The value of the last id assigned, plus 1.

        """

        for mem in self:
            if not hasattr(mem, 'id') or overwrite:
                mem.id = n
                n += 1
        return n

    def dump(
        self,
        path,
        include_attrs=None,
        ignore_missing_attrs=False
    ):
        """
        Dump the :class:`.Population` to a file.

        Parameters
        ----------
        path : :class:`str`
            The full path of the file to which the :class:`Population`
            should be dumped.

        include_attrs : :class:`list` of :class:`str`, optional
            The names of attributes of the molecules to be added to
            the JSON. Each attribute is saved as a string using
            :func:`repr`.

        ignore_missing_attrs : :class:`bool`, optional
            If ``False`` and an attribute in `include_attrs` is not
            held by a :class:`.Molecule`, an error will be raised.

        Returns
        -------
        None : :class:`NoneType`

        """

        with open(path, 'w') as f:
            content = self.to_list(include_attrs, ignore_missing_attrs)
            json.dump(content, f, indent=4)

    @classmethod
    def load(cls, path, use_cache=False):
        """
        Initialize a :class:`Population` from one dumped to a file.

        Parameters
        ----------
        path : :class:`str`
            The full path of the file holding the dumped population.

        use_cache : :class:`bool`, optional
            Toggles use of the moleular cache.

        Returns
        -------
        :class:`Population`
            The population stored in the dump file.

        """

        with open(path, 'r') as f:
            pop_list = json.load(f)

        return cls.init_from_list(pop_list, use_cache)

    def open_process_pool(self, num_processes=None):
        """
        Open a process pool.

        Parameters
        ----------
        num_processes : :class:`int`, optional
            The number of processes in the pool. If ``None``, then
            creates a process for each core on the computer.

        Returns
        -------
        :class:`.Population`
            The population.

        Raises
        ------
        :class:`RuntimeError`
            If a process pool is already open.

        """

        if self._process_pool is not None:
            raise RuntimeError('A process pool is already open.')

        if num_processes is None:
            num_processes = psutil.cpu_count()

        if num_processes != 1:
            self._process_pool = pathos.pools.ProcessPool(
                nodes=num_processes
            )
        return self

    def close_process_pool(self):
        """
        Close an open process pool.

        Returns
        -------
        :class:`.Population`
            The population.

        """

        if self._process_pool is not None:
            self._process_pool.close()
            self._process_pool.clear()
            self._process_pool = None
        return self

    def _optimize_parallel(self, optimizer, num_processes):
        opt_fn = _Guard(optimizer, optimizer.optimize)

        # Only send molecules which need to have a calculation
        # performed to the process pool - this should improve
        # performance.
        if optimizer.is_caching():
            # Use a list here because to_evaluate is iterated through
            # twice.
            to_evaluate = [
                mol for mol in self if not optimizer.is_in_cache(mol)
            ]
        else:
            to_evaluate = self

        # Use an existing process pool, if it exists.
        opened_pool = False
        if self._process_pool is None:
            opened_pool = True
            self.open_process_pool(num_processes)

        # Run the optimization.
        evaluated = self._process_pool.map(opt_fn, to_evaluate)

        if opened_pool:
            self.close_process_pool()

        # Update the structures in the population.
        for input_mol, result in zip(to_evaluate, evaluated):
            if isinstance(result, Exception):
                raise result

            output_mol, _ = result
            input_mol.set_position_matrix(
                position_matrix=output_mol.get_position_matrix()
            )
            if optimizer.is_caching():
                optimizer.add_to_cache(input_mol)

    def _optimize_serial(self, optimizer):
        for member in self:
            optimizer.optimize(member)

    def optimize(self, optimizer, num_processes=None):
        """
        Optimize the structures of molecules in the population.

        The molecules are optimized serially or in parallel depending
        if `num_processes` is ``1`` or more. The serial version may be
        faster in cases where all molecules have already been
        optimized and the `optimizer` will skip them.
        In this case creating a parallel process pool creates
        unnecessary overhead.

        Parameters
        ----------
        optimizer : :class:`.Optimizer`
            The optimizer used to carry out the optimizations.

        num_processes : :class:`int`, optional
            The number of parallel processes to create. Optimization
            will run serially if ``1``. If ``None``, creates a
            process for each core on the computer. This parameter will
            be ignored if the population has an open process pool.

        Returns
        -------
        None : :class:`NoneType`

        """

        if num_processes is None:
            num_processes = psutil.cpu_count()

        if self._process_pool is None and num_processes == 1:
            self._optimize_serial(optimizer)
        else:
            self._optimize_parallel(optimizer, num_processes)

    def remove_duplicates(self, across_subpopulations=True, key=id):
        """
        Remove duplicates from the population.

        The question of which molecule is preserved when duplicates are
        removed is difficult to answer. The iteration through a
        population is depth-first, so a rule such as "the molecule in
        the topmost population is preserved" is not the case here.
        Rather, the first molecule found is preserved.

        However, this question is only relevant if duplicates in
        different subpopulations are being removed. In this case it is
        assumed that it is more important to have a single instance
        than to worry about which subpopulation it is in.

        If the duplicates are being removed from within subpopulations,
        each subpopulation will end up with a single instance of all
        molecules held before. There is no "choice".

        Parameters
        ----------
        across_subpopulations : :class:`bool`, optional
            When ``False`` duplicates are only removed from within a
            given subpopulation. If ``True``, all duplicates are
            removed, regardless of which subpopulation they are in.

        key : :class:`callable`, optional
            Two molecules are considered the same if the values
            returned by ``key(molecule)`` are the same.

        Returns
        -------
        None : :class:`NoneType`

        """

        return self._remove_duplicates(
            across_subpopulations=across_subpopulations,
            key=key,
            seen=set()
        )

    def _remove_duplicates(self, across_subpopulations, key, seen):
        if not across_subpopulations:
            seen = set()

        self.direct_members = list(
            dedupe(self.direct_members, key, seen)
        )
        for subpop in self.subpopulations:
            subpop._remove_duplicates(across_subpopulations, key, seen)

    def remove_members(self, key):
        """
        Remove all members where ``key(member)`` is ``True``.

        Parameters
        ----------
        key : :class:`callable`
            A callable which takes 1 argument. Each member of the
            population is passed as the argument to `key` in turn. If
            the result is ``True`` then the member is removed from the
            population.

        Returns
        -------
        None : :class:`NoneType`

        """

        self.direct_members = [
            ind for ind in self.direct_members if not key(ind)
        ]
        for subpop in self.subpopulations:
            subpop.remove_members(key)

    def to_list(self, include_attrs=None, ignore_missing_attrs=False):
        """
        Convert the population to a :class:`list` representation.

        Parameters
        ----------
        include_attrs : :class:`list` of :class:`str`, optional
            The names of attributes to be added to the molecular
            representations. Each attribute is saved as a string using
            :func:`repr`.

        ignore_missing_attrs : :class:`bool`, optional
            If ``False`` and an attribute in `include_attrs` is not
            held by a :class:`.Molecule`, an error will be raised.

        Returns
        -------
        :class:`list`
            A :class:`list` representation of the :class:`Population`.

        """

        if include_attrs is None:
            include_attrs = []

        pop = [
            m.to_dict(include_attrs, ignore_missing_attrs)
            for m in self.direct_members
        ]
        for sp in self.subpopulations:
            pop.append(sp.to_list(include_attrs, ignore_missing_attrs))
        return pop

    def write(self, path):
        """
        Write the ``.mol`` files of members to a directory.

        Parameters
        ----------
        path : :class:`str`
            The full path of the directory into which the ``.mol`` file
            is written.

        Returns
        -------
        None : :class:`NoneType`

        """

        # If the directory does not exist, create it.
        if not os.path.exists(path):
            os.mkdir(path)

        for i, member in enumerate(self):
            member.write(join(path, f'{i}.mol'))

    def __iter__(self):
        """
        Iterate through members of the :class:`Population`.

        When :class:`Population` instances are iterated through, they
        yield both direct and nested members.

        Yields
        ------
        :class:`.Molecule`
            The next molecule in the population.

        """

        yield from self._get_all_members()

    def __getitem__(self, key):
        """
        Get a member or members by index.

        Molecules held by the :class:`Population` instance can be
        accessed by index. Slices are also supported and they return
        a new :class:`Population` instance holding the members with the
        requested indices. Using slices will return a flat
        :class:`Population` instance, meaing no nesting is preserved.

        Indexing accesses both direct and nested members of the
        :class:`Population`.

        Parameters
        ----------
        key : :class:`int`, :class:`slice`
            An :class:`int` or :class:`slice` can be used depending on
            if a single members needs to be returned or a collection of
            them.

        Returns
        -------
        :class:`.Molecule`
            If the supplied `key` is an :class:`int`. Returns the
            :class:`.Molecule` instance at the corresponding index.

        :class:`Population`
            If the supplied `key` is a :class:`slice`. The returned
            :class:`Population` instance holds members at the given
            indices.

        Raises
        ------
        :class:`IndexError`
            If the supplied `key` is out of range.

        :class:`TypeError`
            If the supplied `key` is not an :class:`int` or
            :class:`slice`.

        """

        if isinstance(key, int):
            for i, member in enumerate(self):
                if i == key:
                    return member
            raise IndexError('Population index out of range.')

        if isinstance(key, slice):
            mols = it.islice(self, key.start, key.stop, key.step)
            return self.__class__(*mols)

        raise TypeError(
            'Index must be an integer or slice, not '
            f'{key.__class__.__name__}.'
        )

    def __len__(self):
        """
        Return the total number of members in the population.

        Returns
        -------
        :class:`int`
            The number of members held by the population, including
            those held within its subpopulations.

        """

        size = len(self.direct_members)
        stack = list(self.subpopulations)
        while stack:
            subpop = stack.pop()
            stack.extend(subpop.subpopulations)
            size += len(subpop.direct_members)
        return size

    def __sub__(self, other):
        """
        Remove members of `other` from the population.

        Subtracting one population from another,

        .. code-block:: python

            pop3 = pop1 - pop2

        returns a new population, ``pop3``. The returned population
        contains all molecules in ``pop1`` except those also found in
        ``pop2``. This refers to all molecules, including those held
        within any subpopulations. The returned population is flat.
        This means any nesting in ``pop1`` is not preserved.

        Parameters
        ----------
        other : :class:`Population`
            A collection of :class:`.Molecule` instances to be removed
            from the population, if held by it.

        Returns
        -------
        :class:`Population`
            A flat population of :class:`.Molecule` instances which
            are in the population but not in `other`.

        """

        new_pop = self.__class__()
        new_pop.add_members(
            molecules=(mol for mol in self if mol not in other)
        )
        return new_pop

    def __add__(self, other):
        """
        Add two populations.

        Parameters
        ----------
        other : :class:`Population`
            A population to be joined with.

        Returns
        -------
        :class:`Population`
            A new :class:`Population` instance which holds two
            subpopulations and no direct members. The two
            subpopulations are the two :class:`Population` instances on
            which the ``+`` operator was applied.

        """

        return self.__class__(self, other)

    def __contains__(self, item):
        return any(mol is item for mol in self)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close_process_pool()

    def __str__(self):
        name = f'Population {id(self)}\n'
        dash = '-'
        members = 'Members'
        subpopulations = 'Subpopulations'

        if self.direct_members:
            direct_members = '\n'.join(
                f'\t{mol}' for mol in self.direct_members
            )
        else:
            direct_members = '\tNone'

        if self.subpopulations:
            subpops = ', '.join(
                str(id(sp)) for sp in self.subpopulations
            )
        else:
            subpops = 'None'

        subpop_strs = '\n'.join(str(sp) for sp in self.subpopulations)

        return (
            f'{name}\n{dash*len(name)}\n'
            f'\t{members}\n\t{dash*len(members)}\n'
            f'{direct_members}\n\n'
            f'\t{subpopulations}\n\t{dash*len(subpopulations)}\n'
            f'{subpops}\n'
            f'{subpop_strs}'
        )

    def __repr__(self):
        return str(self)


class EAPopulation(Population):
    """
    A population which also stores fitness values of molecules.

    Attributes
    ----------
    direct_members : :class:`list` of :class:`.Molecule`
        Held here are direct members of the :class:`Population`.
        In other words, these are the molecules not held by any
        subpopulations. As a result, not all members of a
        :class:`Population` are stored in this attribute.

    subpopulations : :class:`list` of :class:`Population`
        A :class:`list` holding the subpopulations.

    """

    def get_fitness_values(self):
        """
        Return the fitness values of molecules.

        Returns
        -------
        :class:`dict`
            Maps a :class:`.Molecule` to its fitness value.

        """

        return dict(self._fitness_values)

    def set_fitness_values_from_dict(self, fitness_values):
        """
        Set the fitness values of molecules.

        Parameters
        ----------
        fitness_values : :class:`dict`
            Maps molecules in the population to their fitness
            values.

        Returns
        -------
        :class:`.EAPopulation`
            The population is returned.

        """

        self._fitness_values = dict(fitness_values)
        for pop in self.subpopulations:
            pop.set_fitness_values_from_dict(fitness_values)

    def set_fitness_values_from_calculators(
        self,
        fitness_calculator,
        fitness_normalizer=None,
        num_processes=None,
    ):
        """
        Set the fitness values of molecules.

        Parameters
        ----------
        fitness_calculator : :class:`.FitnessCalculator`
            Used to calculate the initial fitness values.

        fitness_normalizer : :class:`.FitnessNormalizer`, optional
            Used to normalize the fitness values.

        num_processes : :class:`int`, optional
            The number of parallel processes to create. Calculations
            will run serially if ``1``. If ``None``, creates a
            process for each core on the computer. This parameter will
            be ignored if the population has an open process pool.

        Returns
        -------
        :class:`.EAPopulation`
            The population is returned.

        """

        if num_processes is None:
            num_processes = psutil.cpu_count()

        if self._process_pool is None and num_processes == 1:
            self._set_fitness_values_serial(fitness_calculator)
        else:
            self._set_fitness_values_parallel(
                fitness_calculator=fitness_calculator,
                num_processes=num_processes,
            )

        if fitness_normalizer is not None:
            self._fitness_values = fitness_normalizer.normalize(self)

        for pop in self.subpopulations:
            pop.set_fitness_values_from_dict(self._fitness_values)

    def _set_fitness_values_serial(self, fitness_calculator):
        self._fitness_values = {
            mol: fitness_calculator.get_fitness(mol)
            for mol in self
        }

    def _set_fitness_values_parallel(
        self,
        fitness_calculator,
        num_processes,
    ):

        fitness_fn = _Guard(
            calculator=fitness_calculator,
            fn=fitness_calculator.get_fitness,
        )

        self._fitness_values = {}
        # Use a list here because the to_evaluate is iterated through
        # twice.
        to_evaluate = list(
            self._handle_cached_mols(fitness_calculator)
        )
        # Use an existing process pool, if it exists.
        opened_pool = False
        if self._process_pool is None:
            opened_pool = True
            self.open_process_pool(num_processes)

        # Apply the function, in parallel.
        evaluated = self._process_pool.map(fitness_fn, to_evaluate)

        if opened_pool:
            self.close_process_pool()

        # Collect results.
        for mol, result in zip(to_evaluate, evaluated):

            if isinstance(result, Exception):
                raise result

            _, fitness = result
            self._fitness_values[mol] = fitness

            if fitness_calculator.is_caching():
                fitness_calculator.add_to_cache(mol, fitness)

        return self

    def _handle_cached_mols(self, fitness_calculator):
        # Only send molecules which need to have a calculation
        # performed to the process pool - this should improve
        # performance.
        if fitness_calculator.is_caching():
            for mol in self:
                if fitness_calculator.is_in_cache(mol):
                    self._fitness_values[mol] = (
                        fitness_calculator.get_fitness(mol)
                    )
                else:
                    yield mol
        else:
            yield from self


class _Guard:
    """
    A decorator for parallelized functions.

    This decorator should be applied to all functions which are to
    be used with :mod:`pathos`. It prevents functions from
    raising if they fail, which prevents the :mod:`pathos` pool
    from hanging.

    """

    def __init__(self, calculator, fn):
        self._calc_name = calculator.__class__.__name__
        wraps(fn)(self)

    def __call__(self, mol):
        """
        Decorates and calls the function.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be passed to the function.

        Returns
        -------
        :class:`tuple`
            The input molecule is the first element of the tuple and
            the value returned by the function is the second element.

        """

        fn = self.__wrapped__.__name__
        cls = self._calc_name
        try:
            logger.info(f'Running "{cls}.{fn}()" on "{mol}"')
            return mol, self.__wrapped__(mol)

        except Exception as ex:
            errormsg = (
                f'"{cls}.{fn}()" failed on molecule "{mol}"'
            )
            logger.error(errormsg, exc_info=True)
            return ex
