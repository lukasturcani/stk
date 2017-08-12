"""
Defines the Population class.

"""

import itertools as it
import os
from os.path import join
import numpy as np
from collections import Counter, defaultdict
import json
from glob import iglob, glob
import multiprocessing as mp

from .fitness import _calc_fitness, _calc_fitness_serial
from .plotting import plot_counter
from .ga_tools import GATools
from ..convenience_tools import dedupe
from ..molecular import (Molecule, Cage,
                         StructUnit, StructUnit2, StructUnit3)
from ..molecular.topologies.cage.base import _VertexOnlyCageTopology
from ..molecular.optimization.optimization import (
                                   _optimize_all_serial, _optimize_all)


class Population:
    """
    A container for  :class:`.Molecule` objects.

    :class:`Population` instances can be nested.

    In addtion to holding :class:`.Molecule` objects, the
    :class:`Population` class can be used to create large numbers of
    these instances through the :meth:`init\_` class methods.

    :class:`.Molecule` instances held by a :class:`Population` can have
    their structures optimized in parallel through the
    :meth:`optimize_population` method.

    The EA is invoked by calling a number of methods of this class,
    such as :meth:`gen_offspring`, :meth:`gen_mutants` and
    :meth:`select`. However, this class only implements container
    related functionality. It delegates EA operations to the
    :class:`.Crossover`, :class:`.Mutation` and :class:`.Selection`
    classes.

    These classes are organised as follows. Each :class:`Population`
    instance has a :attr:`ga_tools` attribute. This holds a
    :class:`.GATools` instance. The :class:`.GATools` instance is just
    a container. It holds a :class:`.Crossover`, :class:`.Mutation`
    and :class:`.Selection` instance. These instances deal with the
    :class:`Population` instance they are held by and perform the
    various crossover, mutation and selection operations on it.
    Any functionality related to the EA should be delegated to these
    instances. The :meth:`gen_offspring` and :meth:`gen_mutants`
    methods can serve as a guide to how this should be done.

    The only operations directly implemented by this class are those
    relevant to its role as a container. It supports all expected and
    necessary container operations such as iteration, indexing,
    membership checks (via the ``is in`` operator) as would be
    expected. Details of the various implementations and a full list of
    supported operations can be found by examining the included
    methods.

    Parameters
    ----------
    *args : :class:`.Molecule`, :class:`Population`, :class:`.GATools`
        A population is initialized with the :class:`.Molecule` and
        :class:`Population` instances it should hold. These are placed
        into the :attr:`members` or :attr:`populations` attributes,
        respectively. A :class:`.GATools` instance may be included
        and will be placed into the :attr:`ga_tools` attribute.

    Raises
    ------
    TypeError
        If initialized with something other than :class:`.Molecule`,
        :class:`Population` or :class:`.GATools` instances.

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

    ga_tools : :class:`.GATools`
        An instance of the :class:`.GATools` class. Calls to preform EA
        operations on the :class:`.Population` instance are delegated
        to this attribute.

    """

    def __init__(self, *args):
        # Generate `populations`, `members` and `ga_tools` attributes.
        self.populations = []
        self.members = []
        self.ga_tools = GATools.init_empty()

        # Determine type of supplied arguments and place in the
        # appropriate attribute.  ``Population`` types added to
        # `populations` attribute, ``Molecule`` into `members` and
        # if ``GATools`` is supplied it is placed into `ga_tools`.
        # Raise a ``TypeError``  if an argument was not ``GATools``,
        # ``Molecule`` or ``Population`` type.
        for arg in args:
            if isinstance(arg, Population):
                self.populations.append(arg)
                continue

            if isinstance(arg, Molecule):
                self.members.append(arg)
                continue

            if isinstance(arg, GATools):
                self.ga_tools = arg
                continue

            raise TypeError(
                    ("Population can only be"
                     " initialized with ``Population``,"
                     " ``Molecule`` and ``GATools`` types."), arg)

    @classmethod
    def init_all(cls, databases, bb_classes,
                 topologies, macromol_class,
                 ga_tools=GATools.init_empty(), duplicates=False):
        """
        Creates all possible molecules from a given set of databases.

        All possible combinations of building blocks from `databases`
        and topologies from `topologies` are built. The initialization
        of a macromolecule in the population has the form

        .. code-block:: python

            mol = macromol_class([bb1, bb2, ...], topology)

        where ``bb1`` is initialized from a  molecule in
        ``databases[0]``, ``bb2`` from a molecule in ``databases[1]``
        and so on.

        Parameters
        ----------
        databases : :class:`list` of :class:`str`
            List of paths to directories, which hold molecular
            structure files of the building blocks.

        bb_classes : :class:`list` of :class:`type`
            This list must be equal in length to `databases`. For each
            database provided in `databases`, a class used to
            initialize building blocks from that database is provided
            here. For example, if

            .. code-block:: python

                databases = ['/path/to/amines2f',
                             '/path/to/aldehydes3f']

            then a valid `bb_classes` list would be

            .. code-block:: python

                bb_classes = [StructUnit2, StructUnit3]

            This means all molecules in the ``amines2f`` database are
            initialized as :class:`.StructUnit2` objects, while all
            molecules in ``aldehydes3f`` are initialized as
            :class:`.StructUnit3` objects.

        topologies : :class:`list` of :class:`.Topology`
            The topologies of macromolecules being made.

        macromol_class : :class:`type`
            The class of the :class:`.MacroMolecule` objects being
            built.

        duplicates : :class:`bool`, optional
            If ``True`` duplicate structures are not removed from
            the population.

        Returns
        -------
        :class:`Population`
            A population holding all possible macromolecules from
            assembled from `databases`.

        Example
        -------
        If the name of the functional group needs to be provided to the
        building blocks, a lambda function can be used.

        .. code-block:: python

            dbs = ['/path/to/db1', 'path/to/db2']
            bb_classes = [lambda x: StructUnit2(x, 'aldehyde'),
                          lambda x: StructUnit3(x, 'amine')]
            tops = [Linear("AB", [0, 0], 6)]
            pop = Population.init_all(dbs, bb_classes, tops, Polymer)

        """

        databases = [glob(join(db, '*')) for db in databases]
        args = []
        for *bb_files, topology in it.product(*databases, topologies):
            bbs = [su(f) for su, f in zip(bb_classes, bb_files)]
            args.append((bbs, topology))

        with mp.Pool() as pool:
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

        p = Population(*mols, ga_tools)
        if not duplicates:
            p.remove_duplicates()
        return p

    @classmethod
    def init_cage_isomers(cls, lk_file, bb_file, topology,
                          ga_tools=GATools.init_empty(),
                          lk_fg=None, bb_fg=None):
        """
        Creates a population holding all structural isomers of a cage.

        Structural isomers here means that the building blocks are
        rotated in position so that every possible bond combination
        with linkers is formed.

        This only works for cage topologies which have both building
        blocks and linkers. It will not with topologies where all
        building blocks have the same number of functional groups.

        Parameters
        ----------
        lk_file : str
            The full path of the file holding the linker of the cage.

        bb_file : str
            The full path to the file holding the building block of the
            cage.

        topology : type
            A _CageTopology child class. Exluding child classes of
            _NoLinkerCageTopology.

        ga_tools : GATools (default = GATools.init_empty())
            The GATools instance to be used by created population.

        lk_fg : str (default = None)
            The name of the linker's functional group. If ``None`` then
            `lk_file` is checked for a name.

        bb_fg : str (default = None)
            The name of the building block's functional group. If
            ``None`` then `bb_file` is checked for a name.

        Returns
        -------
        Population
            A population filled with isomers of a cage.

        """

        n_A = len(topology.positions_A)
        A_alignments = set()

        for x in it.combinations_with_replacement([0, 1, 2], n_A):
            for y in it.permutations(x, n_A):
                A_alignments.add(y)

        n_B = len(topology.positions_B)
        B_alignments = set()
        orientations = ([0, 1, 2] if
                        issubclass(topology, _VertexOnlyCageTopology)
                        else [1, -1])
        for x in it.combinations_with_replacement(orientations, n_B):
            for y in it.permutations(x, n_B):
                B_alignments.add(y)

        lk = StructUnit(lk_file, lk_fg)
        if len(lk.functional_group_atoms()) > 2:
            lk = StructUnit3(lk_file, lk_fg)
        else:
            lk = StructUnit2(lk_file, lk_fg)

        bb = StructUnit3(bb_file, bb_fg)

        pop = cls(ga_tools)
        for A_align in A_alignments:
            for B_align in B_alignments:
                c = Cage([bb, lk], topology(A_align, B_align))
                pop.members.append(c)

        return pop

    @classmethod
    def init_diverse_cages(cls, bb_db, lk_db,
                           topologies, size, ga_tools,
                           bb_fg=None, lk_fg=None):
        """
        Creates a population of cages built from provided databases.

        All cages are held in the population's `members` attribute.

        From the supplied databases a random linker and building block
        molecule is selected to form a cage. The next linker and
        building block selected are those which have the most different
        Morgan fingerprints to the first pair. The next pair random
        again and so on. This is done until `size` cages have been
        formed.

        Parameters
        ----------
        bb_db : str
            The full path of the database of building-block* molecules.

        lk_db : str
            The full path of the database of linker molecules.

        topolgies : iterable of Topology objects
            An iterable holding topologies which should be randomly
            selected for cage initialization.

        size : int
            The size of the population to be initialized.

        ga_tools : GATools
            The GATools instance to be used by created population.

        bb_fg : str (default = None)
            The name of the functional group present in molecules in
            `bb_db`. It is the name of the functional group used to
            build the macromolecules. If ``None`` it is assumed that
            the name is present in `bb_db`.

        lk_fg : str (default = None)
            The name of the functional group present in molecules in
            `lk_db`. It is the name of the functional group used to
            build the macromolecules. If ``None`` it is assumed that
            the name is present in `lk_db`.

        Returns
        -------
        Population
            A population filled with random cages.

        """

        pop = cls(ga_tools)
        bb_files = glob(join(bb_db, '*'))
        # Remove any files which are not valid structure files.
        bb_files = [x for x in bb_files if
                    os.path.splitext(x)[1] in StructUnit.init_funcs]
        lk_files = glob(join(lk_db, '*'))
        lk_files = [x for x in lk_files if
                    os.path.splitext(x)[1] in StructUnit.init_funcs]

        pairs = defaultdict(list)
        bbindices = list(range(len(bb_files)))
        i = -1
        while bbindices:
            i += 1

            topology = np.random.choice(topologies)
            if i % 2 == 0:
                # First pick the index of a building block file.
                bbi = np.random.choice(bbindices)
                # Next get the indices of all linker files which are
                # not already used together with `bbi`.
                lkindices = list(range(len(lk_files)))
                for pairedi in pairs[bbi]:
                    lkindices.remove(pairedi)
                # If `bbi` has been paired with all linkers already,
                # remove it from the list of possible indices and try
                # again.
                if not lkindices:
                    bbindices.remove(bbi)
                    i -= 1
                    continue
                # Pick a linker index and note the pairing.
                lki = np.random.choice(lkindices)
                pairs[bbi].append(lki)

                bb = StructUnit3(bb_files[bbi], bb_fg)
                lk = StructUnit(lk_files[lki], lk_fg)

            else:
                bb_file = bb.similar_molecules(bb_db)[-1][1]
                bb = StructUnit3(bb_file, bb_fg)

                lk_file = lk.similar_molecules(lk_db)[-1][1]
                lk = StructUnit(lk_file, lk_fg)

            if len(lk.bonder_ids) >= 3:
                lk = StructUnit3(lk.file, lk_fg)
            else:
                lk = StructUnit2(lk.file, lk_fg)

            cage = Cage([bb, lk], topology)
            if cage not in pop:
                pop.members.append(cage)

            if len(pop) >= size:
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
            example, :class:`~.StructUnit` or :class:`.StructUnit2`.
            If `folder` contains ``.json`` dump files of
            :class:`MacroMolecule` then :meth:`.Molecule.load` could
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
    def init_random_cages(cls, bb_db, lk_db,
                          topologies, size, ga_tools,
                          bb_fg=None, lk_fg=None):
        """
        Creates a population of cages built from provided databases.

        All cages are held in the population's `members` attribute.

        From the supplied databases a random linker and building-block*
        molecule is selected to form a cage. This is done until `size`
        cages have been formed. After this, all of them are returned
        together in a ``Population`` instance.

        Parameters
        ----------
        bb_db : str
            The full path of the database of building-block* molecules.

        lk_db : str
            The full path of the database of linker molecules.

        topolgies : iterable of Topology objects
            An iterable holding topologies which should be randomly
            selected for cage initialization.

        size : int
            The size of the population to be initialized.

        ga_tools : GATools
            The GATools instance to be used by created population.

        bb_fg : str (default = None)
            The name of the functional group present in molecules in
            `bb_db`. It is the name of the functional group used to
            build the macromolecules. If ``None`` it is assumed that
            the name is present in `bb_db`.

        lk_fg : str (default = None)
            The name of the functional group present in molecules in
            `lk_db`. It is the name of the functional group used to
            build the macromolecules. If ``None`` it is assumed that
            the name is present in `lk_db`.

        Returns
        -------
        Population
            A population filled with random cages.

        """

        pop = cls(ga_tools)
        bb_files = glob(join(bb_db, '*'))
        # Remove any files which are not valid structure files.
        bb_files = [x for x in bb_files if
                    os.path.splitext(x)[1] in StructUnit.init_funcs]
        lk_files = glob(join(lk_db, '*'))
        lk_files = [x for x in lk_files if
                    os.path.splitext(x)[1] in StructUnit.init_funcs]

        pairs = defaultdict(list)
        bbindices = list(range(len(bb_files)))
        while bbindices:
            # First pick the index of a building block file.
            bbi = np.random.choice(bbindices)
            # Next get the indices of all linker files which are not
            # already used together with `bbi`.
            lkindices = list(range(len(lk_files)))
            for pairedi in pairs[bbi]:
                lkindices.remove(pairedi)
            # If `bbi` has been paired with all linkers already, remove
            # it from the list of possible indices and try again.
            if not lkindices:
                bbindices.remove(bbi)
                continue
            # Pick a linker index and note the pairing.
            lki = np.random.choice(lkindices)
            pairs[bbi].append(lki)

            topology = np.random.choice(topologies)
            bb = StructUnit3(bb_files[bbi], bb_fg)
            lk = StructUnit(lk_files[lki], lk_fg)

            if len(lk.bonder_ids) >= 3:
                lk = StructUnit3(lk.file, lk_fg)
            else:
                lk = StructUnit2(lk.file, lk_fg)

            cage = Cage([bb, lk], topology)
            if cage not in pop:
                pop.members.append(cage)

            if len(pop) >= size:
                break

        assert len(pop) == size
        return pop

    def add_members(self, population, duplicates=False):
        """
        Adds ``Molecule`` instances into `members`.

        The |Molecule| instances held within the supplied
        ``Population`` instance, `population`, are added into the
        `members` attribute of `self`. The supplied `population` itself
        is not added. This means that any information the `population`
        instance had about subpopulations is lost. This is because all
        of its ``Molecule`` instances are added into the `members`
        attribute, regardless of which subpopulation they were
        originally in.

        The `duplicates` parameter indicates whether multiple instances
        of the same molecule are allowed to be added into the
        population. Note that the sameness of a molecule is judged
        by the `same` method of the ``Molecule`` class, which is
        invoked by the ``in`` operator within this method. See the
        `__contains__` method of the ``Population`` class for details
        on how the ``in`` operator uses the `same` method.

        Parameters
        ----------
        population : iterable of :class:`~mmead.molecular.molecules.Molecule`
            ``Molecule`` instances to be added to the `members`
            attribute and/or ``Population`` instances who's members, as
            generated by `all_members`, will be added to the `members`
            attribute.

        duplicates : bool (default = False)
            When ``False`` only molecules which are not already
            held by the population will be added. ``True`` allows more
            than one instance of the same molecule to be added.
            Whether two molecules are the same is defined by the
            `same()` method of the ``Molecule`` class.

        Returns
        -------
        None : NoneType

        """

        if duplicates:
            self.members.extend(mol for mol in population)
        else:
            self.members.extend(mol for mol in population if
                                mol not in self)

    def add_subpopulation(self, population):
        """
        Appends a population into the `populations` attribute.

        The `population` instance itself is not added, only a copy.
        However the items it holds are not copied.

        Parameters
        ----------
        population : Population
            The population to be added as a subpopulation.

        Modifies
        --------
        populations : list of Populations
            The `populations` attribute of `self` has ``Population``
            instaces added to it.

        Returns
        -------
        None : NoneType

        """

        pop = Population(*population.members)
        for sp in population.populations:
            pop.add_subpopulation(sp)

        self.populations.append(pop)

    def all_members(self):
        """
        Yields all members in the population and its subpopulations.

        Yields
        ------
        Molecule
            The next ``Molecule`` instance held within the
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

        `n` is a number which is incremented. Each population member
        will have a unique name as a result.

        Parameters
        ----------
        n : int
            A number. Members of this population are given a unique
            number as a name, starting from this number.

        overwrite : bool (default = False)
            If ``True`` existing names are replaced.


        Modifies
        --------
        name : str
            The `name` attribute of the population's members is
            modified.

        Returns
        -------
        int
            The final value of `n`.

        """

        for mem in self:
            if not mem.name or overwrite:
                mem.name = str(n)
                n += 1

        return n

    def calculate_member_fitness(self):
        """
        Applies the fitness function to all members.

        The calculation will be performed serially or in parallel
        depending on the flag `ga_tools.parallel`. The serial version
        may be faster in cases where all molecules have already had
        their fitness calcluated. This is because all calculations will
        be skipped. In this case creating a parallel process pool
        creates unncessary overhead.

        Modifies
        --------
        MarcroMolecule
            The MacroMolecule instances held by the population have
            fitness values calculated and placed in the
            `unscaled_fitness` attribute.

        Returns
        -------
        None : NoneType

        """

        if self.ga_tools.parallel:
            _calc_fitness(self.ga_tools.fitness, self)
        else:
            _calc_fitness_serial(self.ga_tools.fitness, self)

    def dump(self, path):
        """
        Write the population to a file.

        The population is written in the JSON format in the following
        way. The population is represented as a list,

            [mem1.json(), mem2.json(), [mem3.json(), [mem4.json()]]]

        where each member of the population held directly in the
        `members` attribute is placed an an element in the list. Any
        subpopulations are held as  sublists.

        Parameters
        ----------
        path : str
            The full path of the file to which the population should
            be written.

        Returns
        -------
        None : NoneType

        """

        with open(path, 'w') as f:
            json.dump(self.tolist(), f, indent=4)

    def exit(self):
        """
        Checks the if the exit criterion has been satisfied.

        Returns
        -------
        bool
            ``True`` if the exit criterion is satisfied, else
            ``False``.

        """

        return self.ga_tools.exit(self)

    @classmethod
    def fromlist(cls, pop_list, load_names=True):
        """
        Initializes a population from a list representation of one.

        Parameters
        ----------
        pop_list : list of str and lists
            A list which represents a population. Like the ones created
            by `tolist()`.

        load_names : bool (default = True)
            If ``True`` then the `name` attribute stored in the JSON
            objects is loaded. If ``False`` then it's not.

        Returns
        -------
        Population
            The population represented by `pop_list`.

        """

        pop = cls()
        for item in pop_list:
            if isinstance(item, dict):
                pop.members.append(
                    Molecule.fromdict(item, load_names=load_names))
            elif isinstance(item, list):
                pop.populations.append(
                    cls.fromlist(item, load_names=load_names))

            else:
                raise TypeError(('Population list must consist only'
                                 ' of strings and lists.'))
        return pop

    def gen_mutants(self, counter_name='mutation_counter.png'):
        """
        Returns a population of mutant ``MacroMolecule`` instances.

        This is a GA operation and as a result this method merely
        delegates the request to the ``Mutation`` instance held in the
        `ga_tools` attribute.

        Parameters
        ----------
        counter_name : str (default='mutation_counter.png')
            The name of the .png file showing which members were
            selected for mutation.

        Returns
        -------
        Population
            A population holding mutants created by mutating contained
            ``MacroMolecule`` instances.

        """

        return self.ga_tools.mutation(self, counter_name)

    def gen_next_gen(self, pop_size, counter_path=''):
        """
        Returns a population hodling the next generation of structures.

        Parameters
        ----------
        pop_size : int
            The size of the next generation.

        counter_path : str (default= '')
            The name of the .png file showing which members were
            selected for the next generation. If '' then no file
            is made.

        Returns
        -------
        Population
            A population holding the next generation of individuals.

        """

        new_gen = Population(self.ga_tools)
        counter = Counter()
        for member in self.select('generational'):
            counter.update([member])
            new_gen.members.append(member)
            if len(new_gen) == pop_size:
                break

        if counter_path:
            for member in self:
                if member not in counter.keys():
                    counter.update({member: 0})
            plot_counter(counter, counter_path)

        return new_gen

    def gen_offspring(self, counter_name='crossover_counter.png'):
        """
        Returns a population of offspring ``MacroMolecule`` instances.

        This is a GA operation and as a result this method merely
        delegates the request to the ``Crossover`` instance held in the
        `ga_tools` attribute. The ``Crossover`` instance takes care of
        selecting parents and combining them to form offspring. The
        ``Crossover`` instance delegates the selection to the
        ``Selection`` instance as would be expected. The request to
        perform crossovers is done by calling the ``Crossover``
        instance with the population as the argument. Calling of the
        ``Crossover``instance returns a ``Population`` instance holding
        the generated offspring. All details regarding the crossover
        procedure are handled by the ``Crossover`` instance.

        For more details about how crossover is implemented see the
        ``Crossover`` class documentation.

        Parameters
        ----------
        counter_name : str (default='crossover_counter.png')
            The name of the .png file showing which members were
            selected for crossover.

        Returns
        -------
        Population
            A population holding offspring created by crossing
            contained the ``MacroMolecule`` instances.

        """

        return self.ga_tools.crossover(self, counter_name)

    def has_structure(self, mol):
        """
        Returns ``True`` if molecule with `mol` structure is held.

        Parameters
        ----------
        mol : Molecule
            A molecule whose structure is being evaluated for presence
            in the population.

        Returns
        -------
        bool
            ``True`` if a molecule with the same structure as `mol`
            is held by the population.

        """

        return any(x.same(mol) for x in self)

    @classmethod
    def load(cls, path, load_names=True):
        """
        Initializes a Population from one dumped to a file.

        Parameters
        ----------
        path : str
            The full path of the file holding the dumped population.

        load_names : bool (default = True)
            If ``True`` then the `name` attribute stored in the JSON
            objects is loaded. If ``False`` then it's not.

        Returns
        -------
        Population
            The population stored in the dump file.

        """

        with open(path, 'r') as f:
            pop_list = json.load(f)

        pop = cls.fromlist(pop_list, load_names)
        pop.ga_tools = GATools.init_empty()
        return pop

    def max(self, key):
        """
        Calculates the max given a key.

        This method applies key(member) on every member of the
        population and returns the max of returned values.

        For example, if the max value of the attribute `cavity_size`
        was desired:

            population.max(
                 lambda macro_mol : macro_mol.topology.cavity_size())

        Parameters
        ----------
        key : function
            A function which should take a Molecule instance as
            its argument and return a value.

        Returns
        -------
        float
            The max of the values returned by the function `key` when
            its applied to all members of the population.

        """

        return np.max([key(member) for member in self], axis=0)

    def mean(self, key):
        """
        Calculates the mean given a key.

        This method applies key(member) on every member of the
        population and returns the mean of returned values.

        For example, if the mean value of the attribute `cavity_size`
        was desired:

            population.mean(
                 lambda macro_mol : macro_mol.topology.cavity_size())

        Parameters
        ----------
        key : function
            A function which should take a Molecule instance as
            its argument and return a value.

        Returns
        -------
        float
            The mean of the values returned by the function `key` when
            its applied to all members of the population.

        """

        return np.mean([key(member) for member in self], axis=0)

    def min(self, key):
        """
        Calculates the min given a key.

        This method applies key(member) on every member of the
        population and returns the min of returned values.

        For example, if the min value of the attribute `cavity_size`
        was desired:

            population.min(
                 lambda macro_mol : macro_mol.topology.cavity_size())

        Parameters
        ----------
        key : function
            A function which should take a Molecule instance as
            its argument and return a value.

        Returns
        -------
        float
            The min of the values returned by the function `key` when
            its applied to all members of the population.

        """

        return np.min([key(member) for member in self], axis=0)

    def normalize_fitness_values(self):
        """
        Applies the normalization function.

        Returns
        -------
        None : NoneType

        """

        return self.ga_tools.normalization(self)

    def optimize_population(self):
        """
        Optimizes all the members of the population.

        The population is optimized serially or in parallel depending
        on the flag `ga_tools.parallel`. The serial version may be
        faster in cases where all molecules have already been
        optimized. This is because all optimizations will be skipped.
        In this case creating a parallel process pool creates
        unncessary overhead.

        Modifies
        --------
        Molecule
            The Molecule instances held by the population have their
            structures optimized.

        Returns
        -------
        None : NoneType

        """

        if self.ga_tools.parallel:
            _optimize_all(self.ga_tools.optimization, self)
        else:
            _optimize_all_serial(self.ga_tools.optimization, self)

    def remove_duplicates(self, between_subpops=True,
                          key=lambda x: id(x), top_seen=None):
        """
        Removes duplicates from a population and preserves structure.

        The question of which ``Molecule`` instance is preserved
        from a choice of two is difficult to answer. The iteration
        through a population is depth-first, so a rule such as ``the
        molecule in the topmost population is preserved`` is not
        the case here. Rather, the first ``Molecule`` instance
        iterated through is preserved.

        However, this question is only relevant if duplicates in
        different subpopulations are being removed. In this case it is
        assumed that it is more important to have a single instance
        than to worry about which subpopulation it is in.

        If the duplicates are being removed from within subpopulations,
        each subpopulation will end up with a single instance of all
        molecules held before. There is no ``choice``.

        Parameters
        ----------
        between_subpops : bool (default = False)
            When ``False`` duplicates are only removed from within a
            given subpopulation. If ``True`` all duplicates are
            removed, regardless of which subpopulation they are in.

        key : callable (default = lambda x : id)
            Duplicates are removed by on the value returned by this
            function.

        Modifies
        --------
        members
            Duplicate instances are removed from the `members`
            attribute of the population or subpopulations.

        Returns
        -------
        None : NoneType

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
        Removes all members where key(member) is ``True``.

        The structure of the population is preserved.

        Parameters
        ----------
        key : callable
            A callable which takes 1 argument. Each member of the
            population is passed as the argument to `key` in turn. If
            the result is ``True`` then the member is removed from the
            population.

        Modifies
        --------
        self : Population
            All members of the population which have their `failed`
            attribute set to ``True`` are removed.

        Returns
        -------
        None : NoneType

        """

        self.members = [ind for ind in self.members if not key(ind)]
        for subpop in self.populations:
            subpop.remove_members(key)

    def select(self, type_='generational'):
        """
        Returns a generator field yielding selected members of `self`.

        Selection is a GA procedure and as a result this method merely
        delegates the selection request to the ``Selection`` instance
        held within the `ga_tools` attribute. The ``Selection``
        instance then returns a generator which yields
        ``MacroMolecule`` instances held within the population. Which
        macromolecules are yielded depends on the selection algorithm
        which was chosen during initialization and when calling this
        method. The selection instance (`self.ga_tools.selection`)
        returns the generator when it is called. See ``Selection``
        class documentation for more information.

        Because selection is required in a number of different ways,
        such as selecting the parents, ``MacroMolecule`` instances for
        mutation and ``MacroMolecule`` instances for the next
        generation, the type of selection must be specificed with the
        `type_` parameter. The valid values for `type_` will correspond
        to one of the attribute names of the ``Selection`` instance.

        For example, if `type_` is set to 'crossover' a selection
        algorithm which yields a parents will be invoked. If the
        `type_` is set to 'generational' an algorithm which yields the
        next generation will be invoked.

        The information regarding which generational, parent pool, etc.
        algorithm is used is held by the ``Selection`` instance. This
        method merely requests that the ``Selection`` instance performs
        the selection algorithm of the relevant type. The ``Selection``
        instance takes care of all the details to do with selection.

        Parameters
        ----------
        type_ : str (default = 'generational')
            A string specifying the type of selection to be performed.
            Valid values will correspond to names of attributes of the
            ``Selection`` class. Check ``Selection`` class
            documentation for details.

            Valid values include:
                'generational' - selects the next generation
                'crossover' - selects parents
                'mutation' - selects ``MacroMolecule`` instances for
                             mutation

        Returns
        -------
        generator
           A generator which yields ``MacroMolecule`` instances or
           tuples of them. Which instances are yielded depends on the
           selection algorithm used by the generator. This will depend
           on the `type_` provided.

        """

        return self.ga_tools.selection(self, type_)

    def tolist(self):
        """
        Converts the population to a list representation.

        The population and any subpopulations are represented as lists
        (and sublists), while members are represented by their JSON
        dictionaries (as strings).

        Returns
        -------
        str
            A JSON string representing the population.

        """

        pop = [x.json() for x in self.members]
        for sp in self.populations:
            pop.append(sp.tolist())
        return pop

    def write(self, dir_path, use_name=False):
        """
        Writes the ``.mol`` files of members to a directory.

        Parameters
        ----------
        dir_path : str
            The full path of the directory into which the ``.mol`` file
            is written.

        use_name : bool (default = False)
            When ``True`` the `name` attribute of the population's
            members is used to make the name of the .mol file. If
            ``False`` the files are just named after the member's
            index in the population.

        Returns
        -------
        None : NoneType

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
        Allows the use of ``for`` loops, ``*`` and ``iter`` function.

        When ``Population`` instances are iterated through they yield
        ``Molecule`` instances generated by the `all_members`
        method. It also means that a ``Population`` instance can be
        unpacked with the ``*`` operator. This will produce the
        ``Molecule`` instances yielded by the `all_members`
        method.

        Returns
        -------
        Generator
            The `all_members` generator. See `all_members` method
            documentation for more information.

        """

        return self.all_members()

    def __getitem__(self, key):
        """
        Allows the use of ``[]`` operator.

        Molecules held by the ``Population`` instance can be
        accesed by their index. Slices are also supported. These return
        a new ``Population`` instance holding the ``Molecule``
        instances with the requested indices. Using slices will return
        a flat ``Population`` instance meaing no subpopulation
        information is preserved. All of the ``Molecule``
        instances are placed into the `members` attribute of the
        returned ``Population`` instance.

        The index corresponds to the ``Molecule`` yielded by the
        `all_members` method.

        This can be exploited if one desired to remove all
        subpopulations and transfer all the ``Molecules``
        instances into the members attribute. For example,

        >>> pop2 = pop[:]

        ``pop2`` is a ``Population`` instance with all the same
        ``Molecule`` instances as ``pop``, however all
        ``Molecule`` instances are held within its `members`
        attribute and its `populations` attribute is empty. This may or
        may not be the case for the ``pop`` instance.

        Parameters
        ----------
        key : int, slice
            An int or slice can be used depending on if a single
            ``Molecule`` instance needs to be returned or a
            collection of ``Molecule`` instances.

        Returns
        -------
        Molecule
            If the supplied `key` is an ``int``. Returns the
            ``Molecule`` instance with the corresponding index
            from the `all_members` generator.

        Population
            If the supplied `key` is a ``slice``. The returned
            ``Population`` instance holds ``Molecule`` instances
            in its `members` attribute. The ``Molecule`` instances
            correspond to indices defined by the slice. The slice is
            implemented on the `all_members` generator.

        Raises
        ------
        TypeError
            If the supplied `key` is not an ``int`` or ``slice`` type.

        """

        # Determine if provided key was an ``int`` or a ``slice``.
        # If ``int``, return the corresponding ``Molecule``
        # instance from the `all_members` generator.
        if isinstance(key, int):
            return list(self.all_members())[key]

        # If ``slice`` return a ``Population`` of the corresponding
        # ``Molecule`` instances. The returned ``Population`` will
        # have the same `ga_tools` attribute as original ``Population``
        # instance.
        if isinstance(key, slice):
            mols = it.islice(self.all_members(),
                             key.start, key.stop, key.step)
            pop = Population(*mols)
            pop.ga_tools = self.ga_tools
            return pop

        # If `key` is not ``int`` or ``slice`` raise ``TypeError``.
        raise TypeError("Index must be an integer or slice, not"
                        " {}.".format(type(key).__name__))

    def __len__(self):
        """
        Returns the number of members yielded by `all_members`.

        Returns
        -------
        int
            The number of members held by the population, including
            those held within its subpopulations.

        """

        return len(list(self.all_members()))

    def __sub__(self, other):
        """
        Allows use of the ``-`` operator.

        Subtracting one from another,

            pop3 = pop1 - pop2,

        returns a new population, pop3. The returned population
        contains all the ``Molecule`` instances in pop1 except
        those also in pop2. This refers to all of the ``Molecule``
        instances, including those held within any subpopulations. The
        returned population is flat. This means all information about
        subpopulations in pop1 is lost as all the ``Molecule``
        instances are held in the `members` attribute of pop3.

        The resulting population, pop3, will inherit the `ga_tools`
        attribute from pop1.

        Parameters
        ----------
        other : Population
            A collection of ``Molecule`` instances to be removed
            from `self`, if held by it.

        Returns
        -------
        Population
            A flat population of ``Molecule`` instances which are
            not also held in `other`.

        """

        new_pop = Population(self.ga_tools)
        new_pop.add_members(mol for mol in self if mol not in other)
        return new_pop

    def __add__(self, other):
        """
        Allows use fo the ``+`` operator.

        Creates a new ``Population`` instance which holds two
        subpopulations and no direct members. The two subpopulations
        are the two ``Population`` instances on which the ``+``
        operator was applied.

        Parameters
        ----------
        other : Population

        Returns
        -------
        Population


        """

        return Population(self, other, self.ga_tools)

    def __contains__(self, item):
        """
        Allows use of the ``in`` operator.

        Parameters
        ----------
        item : Molecule

        Returns
        -------
        bool

        """

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
