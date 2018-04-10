"""
Defines :class:`GAPopulation`.

:class:`GAPopulation` extends :class:`.Population` by addings tools for
carrying out GA operations.

"""

from collections import Counter
import psutil

from .fitness import _calc_fitness, _calc_fitness_serial
from .plotting import plot_counter
from .ga_tools import GATools
from ..population import Population


class GAPopulation(Population):
    """
    Used for applying GA operations to molecules.

    The GA is invoked by calling a number of methods of this class,
    such as :meth:`gen_offspring`, :meth:`gen_mutants` and
    :meth:`select`. However, this class only implements container
    related functionality. It delegates GA operations to the
    :class:`.Crossover`, :class:`.Mutation` and :class:`.Selection`
    classes.

    These classes are organised as follows. Each :class:`GAPopulation`
    instance has a :attr:`ga_tools` attribute. This holds a
    :class:`.GATools` instance. The :class:`.GATools` instance is just
    a container. It holds a :class:`.Crossover`, :class:`.Mutation`
    and :class:`.Selection` instance. These instances deal with the
    :class:`GAPopulation` instance they are held by and perform the
    various crossover, mutation and selection operations on it.
    Any functionality related to the GA should be delegated to these
    instances. The :meth:`gen_offspring` and :meth:`gen_mutants`
    methods can serve as a guide to how this should be done.

    Attributes
    ----------
    ga_tools : :class:`.GATools`
        An instance of the :class:`.GATools` class. It stores
        instances of classes such as :class:`.Selection`,
        :class:`.Mutation` and :class:`.Crossover`, which carry out GA
        operations on the population.

    """

    def __init__(self, *args, ga_tools=GATools.init_empty()):
        """
        Intializes a population.

        Parameters
        ----------
        *args : :class:`.Molecule`, :class:`Population`
            A population is initialized with the :class:`.Molecule` and
            :class:`Population` instances it should hold. These are
            placed into the :attr:`members` or :attr:`populations`
            attributes, respectively.

        ga_tools : :class:`.GATools`, optional
            The :class:`.GATools` object holding the GA settings via
            the other GA related classes.

        """

        super().__init__(*args)
        self.ga_tools = ga_tools

    @classmethod
    def init_all(cls,
                 macromol_class,
                 building_blocks,
                 topologies,
                 processes=None,
                 duplicates=False,
                 ga_tools=GATools.init_empty()):
        """
        See :meth:`.Population.init_all`.

        Parameters
        ----------
        ga_tools : :class:`.GATools`, optional
            The :class:`.GATools` instance to use.

        """

        p = super().init_all(macromol_class,
                             building_blocks,
                             topologies,
                             processes,
                             duplicates)
        p.ga_tools = ga_tools
        return p

    @classmethod
    def init_diverse(cls,
                     macromol_class,
                     building_blocks,
                     topologies,
                     size,
                     ga_tools=GATools.init_empty()):
        """
        See :meth:`.Population.init_diverse`.

        Parameters
        ----------
        ga_tools : :class:`.GATools`, optional
            The :class:`.GATools` instance to use.

        """

        p = super().init_diverse(macromol_class,
                                 building_blocks,
                                 topologies,
                                 size)
        p.ga_tools = ga_tools
        return p

    @classmethod
    def init_from_files(cls,
                        folder,
                        moltype,
                        glob_pattern='*',
                        ga_tools=GATools.init_empty()):
        """
        See :meth:`.Population.init_from_files`.

        Parameters
        ----------
        ga_tools : :class:`.GATools`, optional
            The :class:`.GATools` instance to use.

        """

        p = super().init_from_files(folder, moltype, glob_pattern)
        p.ga_tools = ga_tools
        return p

    @classmethod
    def init_random(cls,
                    macromol_class,
                    building_blocks,
                    topologies,
                    size,
                    ga_tools):
        """
        See :meth:`.Population.init_random`

        Parameters
        ----------
        ga_tools : :class:`.GATools`, optional
            The :class:`.GATools` instance to use.

        """

        p = super().init_random(macromol_class,
                                building_blocks,
                                topologies,
                                size)
        p.ga_tools = ga_tools
        return p

    def calculate_member_fitness(self, processes=psutil.cpu_count()):
        """
        Applies the fitness function to all members.

        The fitness function is defined in the attribute
        :attr:`.GATools.fitness` of the :class:`.GATools` instance
        held in the :attr:`ga_tools` attribute of the population.

        The calculation will be performed serially or in parallel
        depending on the flag :attr:`.GATools.parallel`. The serial
        version may be faster in cases where all molecules have already
        had their fitness values calcluated. This is because all
        calculations will be skipped. In this case creating a parallel
        process pool creates unncessary overhead.

        Notes
        -----
        This method will result in the modification of the
        :attr:`.MacroMolecule.unscaled_fitness` attribute of molecules
        held by the population.

        Parameters
        ----------
        processes : :class:`int`
            The number of parallel processes to create.

        Returns
        -------
        None : :class:`NoneType`

        """

        if processes == 1:
            _calc_fitness_serial(self.ga_tools.fitness, self)
        else:
            _calc_fitness(self.ga_tools.fitness, self, processes)

    def exit(self, progress):
        """
        Checks the if the EA exit criterion has been satisfied.

        Parameters
        ----------
        progress : :class:`GAPopulation`
            population where each subpopulation is a previous generation.

        Returns
        -------
        :class:`bool`
            ``True`` if the exit criterion is satisfied, else
            ``False``.

        """

        return self.ga_tools.exit(self, progress)

    def gen_mutants(self, counter_name='mutation_counter.png'):
        """
        Returns a :class:`GAPopulation` of mutants.

        The selection function which decides which molecules are
        selected for mutation is defined in the
        :attr:`.Selection.mutation` attribute of the
        :class:`.Selection` instance held by the population via its
        :attr:`ga_tools` attribute.

        The mutation function(s) to be used are defined in the
        :class:`.Mutation` instance held by the population via its
        :attr:`ga_tools` attribute.

        Parameters
        ----------
        counter_name : :class:`str`, optional
            The name of the ``.png`` file which holds a graph showing
            which members were selected for mutation.

        Returns
        -------
        :class:`GAPopulation`
            A population holding mutants created by mutating the
            :class:`.MacroMolecule` instances held by the population.

        """

        return self.ga_tools.mutation(self, counter_name)

    def gen_next_gen(self, pop_size, counter_path=''):
        """
        Returns a population holding the next generation of structures.

        The selection function to be used for selecting the next
        generation of molecules is defined in the :class:`.Selection`
        instance held by the population via its :attr:`ga_tools`
        attribute.

        Parameters
        ----------
        pop_size : :class:`int`
            The size of the next generation.

        counter_path : :class:`str`, optional
            The name of the ``.png`` file holding a graph showing which
            members were selected for the next generation. If ``''``
            then no file is made.

        Returns
        -------
        :class:`GAPopulation`
            A population holding the next generation of molecules.

        """

        new_gen = GAPopulation(ga_tools=self.ga_tools)
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
        Returns a :class:`GAPopulation` of offspring molecules.

        The selection function which decides which molecules are
        selected for crossover is defined in the
        :attr:`.Selection.crossover` attribute of the
        :class:`.Selection` instance held by the population via its
        :attr:`ga_tools` attribute.

        The crossover function(s) to be used are defined in the
        :class:`.Crossover` instance held by the population via its
        :attr:`ga_tools` attribute.

        Parameters
        ----------
        counter_name : :class:`str`, optional
            The name of the ``.png`` file showing which members were
            selected for crossover.

        Returns
        -------
        :class:`GAPopulation`
            A population of offspring, created by crossing
            :class:`.MacroMolecule` instances contained in the
            population.

        """

        return self.ga_tools.crossover(self, counter_name)

    def normalize_fitness_values(self):
        """
        Applies the normalization functions on the population.

        Normalization functions scale or modify the fitness values of
        molecules in the population.

        The normalization functions which are applied on the
        population, along with their order, are defined in the
        :class:`.Normalization` instance held by the population via the
        :attr:`ga_tools` attribute.

        Notes
        -----
        This method modifies the :attr:`.MacroMolecule.fitness`
        attribute of molecules held by the population.

        Returns
        -------
        None : :class:`NoneType`

        """

        return self.ga_tools.normalization(self)

    def select(self, type='generational'):
        """
        Returns a generator for yielding members of the population.

        Members are yielded based on the selection function defined in
        the :class:`.Selection` instance held by the population via the
        :attr:`ga_tools` attribute. The instance defines 3 selection
        functions, 1 for each `type`.

        Parameters
        ----------
        type : :class:`str`, optional
            A string specifying the type of selection to be performed.
            Valid values will correspond to names of attributes of the
            :class:`.Selection` class.

            Valid values are:

                * ``'generational'`` - selects the next generation
                * ``'crossover'`` - selects molecules for crossover
                * ``'mutation'`` - selects molecules for mutation

        Returns
        -------
        :class:`generator`
           A generator which yields molecules or tuples of them.

        """

        return self.ga_tools.selection(self, type)

    def __getitem__(self, key):
        p = super().__getitem__(key)
        p.ga_tools = self.ga_tools
        return p

    def __sub__(self, other):
        p = super().__sub__(other)
        p.ga_tools = self.ga_tools
        return p

    def __add__(self, other):
        p = super().__add__(other)
        p.ga_tools = self.ga_tools
        return p
