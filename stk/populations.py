"""
Defines :class:`Population`.

"""

import itertools as it
import os
from os.path import join
import numpy as np
import json
from glob import iglob
import multiprocessing as mp
import psutil
from functools import wraps
import logging
from threading import Thread

from .utilities import dedupe, OPTIONS, daemon_logger, logged_call


logger = logging.getLogger(__name__)


class Population:
    """
    A container for  :class:`.Molecule` objects.

    :class:`Population` instances can be nested.

    In addtion to holding :class:`.Molecule` objects, the
    :class:`Population` class can be used to create large numbers of
    these instances through the :meth:`init\_` class methods.

    :class:`.Molecule` instances held by a :class:`Population` can have
    their structures optimized in parallel through the
    :meth:`optimize` method.

    The only operations directly implemented by this class are those
    relevant to its role as a container. It supports all expected and
    necessary container operations such as iteration, indexing,
    membership checks (via the ``is in`` operator) as would be
    expected. Details of the various implementations and a full list of
    supported operations can be found by examining the included
    methods.

    Attributes
    ----------
    populations : :class:`list` of :class:`Population` instances
        A list of other instances of the :class:`Population` class.
        This allows the implementation of subpopulations or
        evolutionary islands. This attribute is also used for grouping
        molecules within a given population for organizational
        purposes.

    members : :class:`list` of :class:`.Molecule` instances
        Held here are the members of the population which are not held
        within any subpopulations. This means that not all members of a
        population are stored in this attribute. To access all members
        of a population the generator :meth:`all_members` should be
        used.

    """

    def __init__(self, *args):
        """
        Intializes a population.

        Parameters
        ----------
        *args : :class:`.Molecule`, :class:`Population`
            A population is initialized with the :class:`.Molecule` and
            :class:`Population` instances it should hold. These are
            placed into the :attr:`members` or :attr:`populations`
            attributes, respectively.

        """

        self.populations = []
        self.members = []

        for arg in args:
            if isinstance(arg, Population):
                self.populations.append(arg)
            else:
                self.members.append(arg)

    @classmethod
    def init_all(cls,
                 macromol_class,
                 building_blocks,
                 topologies,
                 processes=None,
                 duplicates=False):
        """
        Creates all possible molecules from provided building blocks.

        Parameters
        ----------
        macromol_class : :class:`type`
            The class of the :class:`.MacroMolecule` objects being
            built.

        building_blocks : :class:`list`
            A :class:`list` of the form

            .. code-block:: python

                building_blocks = [[StructUnit2(), StructUnit2(), ...],
                                   [StructUnit3(), StructUnit3(), ...],
                                   [StructUnit2(), StructUnit2(), ...]]

            To assemble a new :class:`.MacroMolecule`, a
            :class:`.StructUnit` is picked from each of the sublists
            in `building_blocks`. The picked :class:`.StructUnit`
            instances are then supplied to the macromolecule:

            .. code-block:: python

                macro_mol = MacroMolecule([pick1, pick2, pick3],
                                          Topology())

            The order of picked :class:`.StructUnit` instances
            corresponds to the order of the sublists.

        topologies : :class:`list` of :class:`.Topology`
            The topologies of macromolecules being made.

        processes : :class:`int`, optional
            The number of parallel processes to create when building
            the molecules.

        duplicates : :class:`bool`, optional
            If ``False``, duplicate structures are removed from
            the population.

        Returns
        -------
        :class:`Population`
            A population holding all possible macromolecules from
            assembled from `databases`.

        """

        args = []
        for *bbs, topology in it.product(*building_blocks, topologies):
            args.append((bbs, topology))

        with mp.Pool(processes) as pool:
            mols = pool.starmap(macromol_class, args)

        # Update the cache.
        for i, mol in enumerate(mols):
            # If the molecule did not exist already add it to the
            # cache.
            if mol.key not in macromol_class.cache:
                macromol_class.cache[mol.key] = mol
            # If the molecule did exist already, use the cached
            # version.
            else:
                mols[i] = macromol_class.cache[mol.key]

        p = cls(*mols)
        if not duplicates:
            p.remove_duplicates()
        return p

    @classmethod
    def init_copy(cls, population):
        """
        Makes a copy of `population`

        Molecules in `population` are not copied.

        Parameters
        ----------
        population : :class:`Population`
            The population to copy.

        Returns
        -------
        :class:`Population`
            A copy of `population`.

        """

        copy = cls(*population.members)
        copy.populations = [
            cls.init_copy(pop) for pop in population.populations
        ]
        return copy

    @classmethod
    def init_diverse(cls,
                     macromol_class,
                     building_blocks,
                     topologies,
                     size):
        """
        Assembles a population of :class:`.MacroMolecule`.

        All molecules are held in the :attr:`members`.

        From the supplied sublists of building blocks, a random
        molecule is selected to initialize a :class:`.MacroMolecule`
        per sublist. The next molecule selected from the same sublist
        is the one most with the most different Morgan fingerprint. The
        next molecule is picked at random again and so on. This is done
        until `size` :class:`.MacroMolecule` instances have been
        formed.

        Parameters
        ----------
        macromol_class : :class:`type`
            The class of :class:`.MacroMolecule` to be assembled.

        building_blocks : :class:`list`
            A :class:`list` of the form

            .. code-block:: python

                building_blocks = [[StructUnit2(), StructUnit2(), ...],
                                   [StructUnit3(), StructUnit3(), ...],
                                   [StructUnit2(), StructUnit2(), ...]]

            To assemble a new :class:`.MacroMolecule`, a
            :class:`.StructUnit` is picked from each of the sublists
            in `building_blocks`. The picked :class:`.StructUnit`
            instances are then supplied to the macromolecule:

            .. code-block:: python

                macro_mol = MacroMolecule([pick1, pick2, pick3],
                                          Topology())

            The order of picked :class:`.StructUnit` instances
            corresponds to the order of the sublists.

        topolgies : :class:`iterable` of :class:`.Topology`
            An iterable holding topologies which should be randomly
            selected during initialization of :class:`.MacroMolecule`.

        size : :class:`int`
            The size of the population to be initialized.

        Returns
        -------
        :class:`.Population`
            A population filled with generated molecules.

        """

        pop = cls()

        # Shuffle the sublists.
        for db in building_blocks:
            np.random.shuffle(db)

        # Go through every possible macromolecule.
        for *bbs, top in it.product(*building_blocks, topologies):

            # Generate the random macromolecule.
            macro_mol = macromol_class(bbs, top)
            if macro_mol not in pop:
                pop.members.append(macro_mol)

            if len(pop) == size:
                break

            # Make an iterators which goes through all rdkit molecules
            # in the sublists.
            mol_iters = [(struct_unit.mol for struct_unit in db) for
                         db in building_blocks]
            # Make a dictionary which maps every rdkit molecule to its
            # StructUnit, for every sublist in building_blocks.
            mol_maps = [{struct_unit.mol: struct_unit for struct_unit in db}
                        for db in building_blocks]

            # Get the most different StructUnit to the previously
            # selected one, per sublist. Take index of 1 because the
            # index of 0 will the molecule itself.
            diff_mols = [bb.similar_molecules(mols)[1][1] for
                         bb, mols in zip(bbs, mol_iters)]
            diff_bbs = [mol_map[mol] for
                        mol_map, mol in zip(mol_maps, diff_mols)]

            macro_mol = macromol_class(diff_bbs, top)
            if macro_mol not in pop:
                pop.members.append(macro_mol)

            if len(pop) == size:
                break

        assert len(pop) == size
        return pop

    @classmethod
    def init_from_files(cls, folder, moltype, glob_pattern='*'):
        """
        Creates a population from files in `folder`.

        Parameters
        ----------
        folder : :class:`str`
            The path to a folder holding molecular structure files
            used to initialize :class:`~.Molecule` objects held by
            the population.

        moltype : :class:`type`
            An initializer for the molecular structure files. For
            example, :class:`.StructUnit` or :class:`.StructUnit2`.
            If `folder` contains ``.json`` dump files of
            :class:`.MacroMolecule` then :meth:`.Molecule.load` could
            also be used.

        glob_pattern : :class:`str`, optional
            A glob used for selecting specific files within `folder`.

        Returns
        -------
        :class:`Population`
            A population made from files in `folder`.

        """

        return cls(*(moltype(x) for x in
                     iglob(join(folder, glob_pattern))))

    @classmethod
    def init_random(cls,
                    macromol_class,
                    building_blocks,
                    topologies,
                    size):
        """
        Assembles a population of :class:`.MacroMolecule`.

        All molecules are held in :attr:`members`.

        From the supplied building blocks a random molecule is selected
        per sublist to for a :class:`.MacroMolecule`. This is done
        until `size` :class:`.MacroMolecule` have been formed.

        Parameters
        ----------
        macromol_class : :class:`type`
            The class of :class:`.MacroMolecule` to be assembled.

        building_blocks : :class:`list`
            A :class:`list` of the form

            .. code-block:: python

                building_blocks = [[StructUnit2(), StructUnit2(), ...],
                                   [StructUnit3(), StructUnit3(), ...],
                                   [StructUnit2(), StructUnit2(), ...]]

            To assemble a new :class:`.MacroMolecule`, a
            :class:`.StructUnit` is picked from each of the sublists
            in `building_blocks`. The picked :class:`.StructUnit`
            instances are then supplied to the macromolecule:

            .. code-block:: python

                macro_mol = MacroMolecule([pick1, pick2, pick3],
                                          Topology())

            The order of picked :class:`.StructUnit` instances
            corresponds to the order of the sublists.

        topolgies : :class:`iterable` of :class:`.Topology`
            An iterable holding topologies which should be randomly
            selected during initialization of :class:`.MacroMolecule`.

        size : :class:`int`
            The size of the population to be initialized.

        Returns
        -------
        :class:`.Population`
            A population filled with random cages.

        """

        pop = cls()

        # Shuffle the sublists.
        for db in building_blocks:
            np.random.shuffle(db)

        # Go through every possible macromolecule.
        for *bbs, top in it.product(*building_blocks, topologies):

            # Generate the random macromolecule.
            macro_mol = macromol_class(bbs, top)
            if macro_mol not in pop:
                pop.members.append(macro_mol)

            if len(pop) == size:
                break

        assert len(pop) == size
        return pop

    def add_members(self, population, duplicates=False):
        """
        Adds :class:`.Molecule` instances into :attr:`members`.

        The :class:`.Molecule`  instances held within `population`, are
        added into :attr:`members`. Any nesting is removed. This is
        because all :class:`.Molecule`  instances are added into
        :attr:`members` directly, regardless of how nested they are
        within `population`.

        Parameters
        ----------
        population : :class:`iterable` of :class:`.Molecule`
            :class:`.Molecule` instances held by this container are
            added into :attr:`members`.

        duplicates : :class:`bool`, optional
            Indicates whether multiple instances of the same molecule
            are allowed to be added. Note, the sameness of a molecule
            is judged by :meth:`.Molecule.same`.

            When ``False`` only molecules which are not already
            held by the population will be added. ``True`` allows more
            than one instance of the same molecule to be added.

        Returns
        -------
        None : :class:`NoneType`

        """

        if duplicates:
            self.members.extend(mol for mol in population)
        else:
            self.members.extend(mol for mol in population if
                                mol not in self)

    def add_subpopulation(self, population):
        """
        Appends a copy of `population` to :attr:`populations`.

        Only a copy of the `population` container is made. The
        molecules it holds are not copies.

        Parameters
        ----------
        population : :class:`Population`
            The population to be added as a subpopulation.

        Returns
        -------
        None : :class:`NoneType`

        """

        pop = self.__class__(*population.members)
        for sp in population.populations:
            pop.add_subpopulation(sp)

        self.populations.append(pop)

    def all_members(self):
        """
        Yields all molecules in the population and its subpopulations.

        Yields
        ------
        :class:`.Molecule`
            The next :class:`.Molecule` instance held within the
            population or its subpopulations.

        """

        # Go through `members` attribute and yield ``Molecule``
        # instances held within one by one.
        for ind in self.members:
            yield ind

        # Go thorugh `populations` attribute and for each
        # ``Population`` instance within, yield ``Molecule``
        # instances from its `all_members` generator.
        for pop in self.populations:
            yield from pop.all_members()

    def assign_names_from(self, n, overwrite=False):
        """
        Give each member of the population a name starting from `n`.

        Notes
        -----
        This method modifies the :attr:`.Molecule.name` attribute of
        :class:`.Molecule` instances held by the population.

        Parameters
        ----------
        n : :class:`int`
            A number. Members of this population are given a unique
            number as a name, starting from this number and incremented
            by one between members.

        overwrite : :class:`bool`, optional
            If ``True``, existing names are replaced.

        Returns
        -------
        :class:`int`
            The value of the last name assigned, plus 1.

        """

        for mem in self:
            if not mem.name or overwrite:
                mem.name = str(n)
                n += 1

        return n

    def dump(self, path, include_attrs=None):
        """
        Dumps the population to a file.

        The population is dumped in the JSON format in the following
        way. The population is represented as a list,

        .. code-block:: python

            [mem1.json(), mem2.json(), [mem3.json(), [mem4.json()]]]

        where each member of the population held directly in the
        `members` attribute is placed an an element in the list. Any
        subpopulations are held as sublists.

        Parameters
        ----------
        path : :class:`str`
            The full path of the file to which the population should
            be dumped.

        include_attrs : :class:`list` of :class:`str`, optional
            The names of attributes of the molecules to be added to
            the JSON. Each attribute is saved as a string using
            :func:`repr`.

        Returns
        -------
        None : :class:`NoneType`

        """

        with open(path, 'w') as f:
            json.dump(self.to_list(include_attrs), f, indent=4)

    @classmethod
    def from_list(cls, pop_list, member_init):
        """
        Initializes a population from a :class:`list` representation.

        Parameters
        ----------
        pop_list : :class:`list`
            A :class:`list` which represents a population. Like the
            ones created by :meth:`to_list`. For example in,

            .. code-block:: python

                pop_list = [{...}, [{...}], [{...}, {...}], {...}]

            ``pop_list`` represents the population, sublists represent
            its subpopulations and the :class:`dict` ``{...}``
            represents the members.

        member_init : :class:`function`
            The initialization function for the population's members.
            It converts the member represenations in `pop_list` into
            desired objects.

        Returns
        -------
        :class:`Population`
            The population represented by `pop_list`.

        """

        pop = cls()
        for item in pop_list:
            if isinstance(item, dict):
                pop.members.append(member_init(item))
            elif isinstance(item, list):
                pop.populations.append(cls.from_list(item, member_init))

            else:
                raise TypeError(('Population list must consist only'
                                 ' of strings and lists.'))
        return pop

    def has_structure(self, mol):
        """
        Returns ``True`` if molecule with `mol` structure is held.

        Parameters
        ----------
        mol : :class:`.Molecule`
            A molecule whose structure is being evaluated for presence
            in the population.

        Returns
        -------
        :class:`bool`
            ``True`` if a molecule with the same structure as `mol`
            is held by the population.

        """

        return any(x.same(mol) for x in self)

    @classmethod
    def load(cls, path, member_init):
        """
        Initializes a :class:`Population` from one dumped to a file.

        Parameters
        ----------
        path : :class:`str`
            The full path of the file holding the dumped population.

        member_init : :class:`function`
            The initialization function for the population's members.
            It converts the member representations in the file to
            :class:`.Molecule` objects. For example
            :meth:`.Molecule.from_dict` when loading ``.json`` files
            generated by :meth:`dump`.

        Returns
        -------
        :class:`Population`
            The population stored in the dump file.

        """

        with open(path, 'r') as f:
            pop_list = json.load(f)

        return cls.from_list(pop_list, member_init)

    def max(self, key):
        """
        Calculates the maximum in the population given a key.

        This method applies ``key(member)`` on every member of the
        population and returns the maximum of returned values.

        For example, if the maximum molecular cavity in the
        population is desired:

        .. code-block:: python

            population.max(lambda macro_mol: macro_mol.cavity_size())

        Parameters
        ----------
        key : :class:`function`
            A function which should take a :class:`.Molecule` instance
            as its argument and returns a numerical value.

        Returns
        -------
        :class:`float`
            The maximum of the values returned by the function `key`
            when it is applied to all members of the population.

        """

        return np.max([key(member) for member in self], axis=0)

    def mean(self, key):
        """
        Calculates the mean in the population given a key.

        This method applies ``key(member)`` on every member of the
        population and returns the mean of returned values.

        For example, if the mean of molecular cavities in the
        population is desired:

        .. code-block:: python

            population.mean(lambda macro_mol: macro_mol.cavity_size())

        Parameters
        ----------
        key : :class:`function`
            A function which should take a :class:`.Molecule` instance
            as its argument and returns a numerical value.

        Returns
        -------
        :class:`float`
            The mean of the values returned by the function `key` when
            it is applied to all members of the population.

        """

        return np.mean([key(member) for member in self], axis=0)

    def min(self, key):
        """
        Calculates the minimum in the population given a key.

        This method applies ``key(member)`` on every member of the
        population and returns the minimum of returned values.

        For example, if the minimum molecular cavity in the
        population is desired:

        .. code-block:: python

            population.min(lambda macro_mol: macro_mol.cavity_size())

        Parameters
        ----------
        key : :class:`function`
            A function which should take a :class:`.Molecule` instance
            as its argument and returns a numerical value.

        Returns
        -------
        :class:`float`
            The minimum of the values returned by the function `key`
            when it is applied to all members of the population.

        """

        return np.min([key(member) for member in self], axis=0)

    def _optimize_parallel(self, optimizer, processes):
        manager = mp.Manager()
        logq = manager.Queue()
        log_thread = Thread(target=daemon_logger, args=(logq, ))
        log_thread.start()

        opt_fn = _Guard(optimizer, optimizer.optimize)

        # Apply the function to every member of the population, in
        # parallel.
        with mp.get_context('spawn').Pool(processes) as pool:
            optimized = pool.starmap(logged_call,
                                     ((logq, opt_fn, mem) for
                                      mem in self))

        logq.put(None)
        log_thread.join()

        # If anything failed, raise an error.
        for result in optimized:
            if isinstance(result, Exception):
                raise result

        # Update the structures in the population.
        sorted_opt = sorted(optimized, key=lambda m: m.key)
        sorted_pop = sorted(self, key=lambda m: m.key)
        for old, new in zip(sorted_pop, sorted_opt):
            old.__dict__ = dict(vars(new))
            if optimizer.use_cache:
                optimizer.cache.add((old.key, -1))

        # Make sure the cache is updated with the optimized versions.
        if OPTIONS['cache']:
            for member in optimized:
                member.update_cache()

    def _optimize_serial(self, optimizer):
        for member in self:
            optimizer.optimize(member)

    def optimize(self, optimizer, processes=psutil.cpu_count()):
        """
        Optimizes the structures of molecules in the population.

        The molecules are optimized serially or in parallel depending
        if `processes` is ``1`` or more. The serial version may be
        faster in cases where all molecules have already been
        optimized, for example if optimized molecules are skipped.
        In this case creating a parallel process pool creates
        unncessary overhead.

        Notes
        -----
        This function modifies the structures of molecules held by the
        population. This means their :attr:`.Molecule.mol` attributes
        are modified.

        Parameters
        ----------
        optimizer : :class:`.Optimizer`
            The optimizer used to carry out the optimizations.

        processes : :class:`int`
            The number of parallel processes to create. Optimization
            will run serially if ``1``.

        Returns
        -------
        None : :class:`NoneType`

        """

        if processes == 1:
            self._optimize_serial(optimizer)
        else:
            self._optimize_parallel(optimizer, processes)

    def remove_duplicates(self,
                          between_subpops=True,
                          key=id,
                          top_seen=None):
        """
        Removes duplicates from the population and preserves structure.

        The question of which molcule is preserved when duplicates are
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
        between_subpops : :class:`bool`, optional
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

        # Whether duplicates are being removed from within a single
        # subpopulation or from different subpopulations, the duplicate
        # must be removed from the `members` attribute of some
        # ``Population`` instance. This means ``dedupe`` can be run
        # on the `members` attribute of every population or
        # subpopulation. The only difference is that when duplicates
        # must be removed between different subpopulations a global
        # ``seen`` set must be defined for the entire top level
        # ``Population`` instance. This can be passed each time dedupe
        # is being called on a subpopulation's `members` attribute.
        if between_subpops:
            if top_seen is None:
                seen = set()
            if isinstance(top_seen, set):
                seen = top_seen

            self.members = list(dedupe(self.members, seen, key))
            for subpop in self.populations:
                subpop.remove_duplicates(True, key, seen)

        # If duplicates are only removed from within the same
        # subpopulation, only the `members` attribute of each
        # subpopulation needs to be cleared of duplicates. To do this,
        # each `members` attribute is deduped recursively.
        if not between_subpops:
            self.members = list(dedupe(self.members, key=key))
            for subpop in self.populations:
                subpop.remove_duplicates(False, key)

    def remove_members(self, key):
        """
        Removes all members where ``key(member)`` is ``True``.

        The structure of the population is preserved.

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

        self.members = [ind for ind in self.members if not key(ind)]
        for subpop in self.populations:
            subpop.remove_members(key)

    def to_list(self, include_attrs=None):
        """
        Converts the population to a list representation.

        The population and any subpopulations are represented as lists
        (and sublists), while members are represented by their JSON
        dictionaries (as strings).

        Parameters
        ----------
        include_attrs : :class:`list` of :class:`str`, optional
            The names of attributes of the molecules to be added to
            the JSON. Each attribute is saved as a string using
            :func:`repr`.

        Returns
        -------
        :class:`str`
            A JSON string representing the population.

        """

        if include_attrs is None:
            include_attrs = []

        include_attrs = set(include_attrs)
        pop = [
            m.json(list(include_attrs & m.__dict__.keys()))
            for m in self.members
        ]
        for sp in self.populations:
            pop.append(sp.to_list(include_attrs))
        return pop

    def write(self, dir_path, use_name=False):
        """
        Writes the ``.mol`` files of members to a directory.

        Parameters
        ----------
        dir_path : :class:`str`
            The full path of the directory into which the ``.mol`` file
            is written.

        use_name : :class:`bool`, optional
            When ``True`` the :attr:`.Molecule.name` attribute of the
            members is used to make the name of the ``.mol`` file. If
            ``False``, the files are just named after the member's
            index in the population.

        Returns
        -------
        None : :class:`NoneType`

        """

        # If the directory does not exist, create it.
        if not os.path.exists(dir_path):
            os.mkdir(dir_path)

        for i, member in enumerate(self):
            if use_name:
                fname = join(dir_path, '{}.mol'.format(
                                                        member.name))
            else:
                fname = join(dir_path, '{}.mol'.format(i))

            member.write(fname)

    def __iter__(self):
        """
        Allows populations to be iterated through,

        When :class:`Population` instances are iterated through, they
        yield :class:`.Molecule` instances from the :meth:`all_members`
        generator.

        Returns
        -------
        :class:`generator`
            The :meth:`all_members` generator.

        """

        return self.all_members()

    def __getitem__(self, key):
        """
        Allows indexing to select members of the population.

        Molecules held by the :class:`Population` instance can be
        accessed by their index. Slices are also supported. These return
        a new :class:`Population` instance holding the members with the
        requested indices. Using slices will return a flat
        :class:`Population` instance, meaing no nesting is preserved.

        The index corresponds to the order of members in
        :meth:`all_members`.

        Parameters
        ----------
        key : :class:`int`, :class:`slice`
            An :class:`int` or :class:`slice` can be used depending on
            if a single members needs to be returned or a cllection of
            them.

        Returns
        -------
        :class:`.Molecule`
            If the supplied `key` is an :class:`int`. Returns the
            :class:`.Molecule` instance with the corresponding index
            from :meth:`all_members`.

        :class:`Population`
            If the supplied `key` is a :class:`slice`. The returned
            :class:`Population` instance holds members at the given
            indices in :meth:`all_members`.

        Raises
        ------
        :class:`TypeError`
            If the supplied `key` is not an :class:`int` or
            :class:`slice`.

        """

        # Determine if provided key was an ``int`` or a ``slice``.
        # If ``int``, return the corresponding ``Molecule``
        # instance from the `all_members` generator.
        if isinstance(key, int):
            return list(self.all_members())[key]

        # If ``slice`` return a ``Population`` of the corresponding
        # ``Molecule`` instances.
        if isinstance(key, slice):
            mols = it.islice(self.all_members(),
                             key.start, key.stop, key.step)
            pop = self.__class__(*mols)
            return pop

        # If `key` is not ``int`` or ``slice`` raise ``TypeError``.
        raise TypeError("Index must be an integer or slice, not"
                        " {}.".format(type(key).__name__))

    def __len__(self):
        """
        Returns the total number of members in the population.

        Returns
        -------
        :class:`int`
            The number of members held by the population, including
            those held within its subpopulations.

        """

        size = len(self.members)
        stack = list(self.populations)
        while stack:
            subpop = stack.pop()
            stack.extend(subpop.populations)
            size += len(subpop.members)
        return size

    def __sub__(self, other):
        """
        Removes members of `other` from the population.

        Subtracting one population from another,

        .. code-block:: python

            pop3 = pop1 - pop2

        returns a new population, ``pop3``. The returned population
        contains all molecules in ``pop1`` excep those also found in
        ``pop2``. This refers to all molcules, including those held
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
        new_pop.add_members(mol for mol in self if mol not in other)
        return new_pop

    def __add__(self, other):
        """
        Joins two populations.

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
        return any(item is mol for mol in self.all_members())

    def __str__(self):
        output_string = (" Population " + str(id(self)) + "\n" +
                         "--------------------------\n" +
                         "\tMembers\n" + "   ---------\n")

        for mol in self.members:
            output_string += "\t" + str(mol) + "\n"

        if len(self.members) == 0:
            output_string += "\tNone\n\n"

        output_string += (("\tSub-populations\n" +
                           "   -----------------\n\t"))

        for pop in self.populations:
            output_string += str(id(pop)) + ", "

        if len(self.populations) == 0:
            output_string += "None\n\n"

        output_string += "\n\n"

        for pop in self.populations:
            output_string += str(pop)

        return output_string

    def __repr__(self):
        return str(self)


class _Guard:
    """
    A decorator for optimization functions.

    This decorator should be applied to all functions which are to
    be used with :mod:`multiprocessing`. It prevents functions from
    raising if they fail, which prevents the multiprocessing pool
    from hanging.

    """

    def __init__(self, calculator, fn):
        self.calc_name = calculator.__class__.__name__
        wraps(fn)(self)

    def __call__(self, mol):
        """
        Decorates and calls the function.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        Returns
        -------
        :class:`.Molecule`
            The optimized molecule.

        """

        fn = self.__wrapped__.__name__
        cls = self.calc_name
        try:
            logger.info(f'Running "{cls}.{fn}()" on "{mol.name}"')
            self.__wrapped__(mol)
            return mol

        except Exception as ex:
            errormsg = (
                f'"{cls}.{fn}()" failed on molecule "{mol.name}"'
            )
            logger.error(errormsg, exc_info=True)
            return ex
