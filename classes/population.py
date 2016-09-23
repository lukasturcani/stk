import itertools
import shutil
import os
from statistics import mean
import pickle

from .molecular import MacroMolecule, Cage
from .ga import GATools
from ..convenience_functions import dedupe
from ..optimization import optimize_all, optimize_all_serial
from ..fitness import calc_fitness

class Population:
    """
    A container for instances of ``MacroMolecule`` and ``Population``.

    This is the central class of MMEA. The GA is invoked by calling the
    ``gen_offspring``, ``gen_mutants`` and ``select`` methods of this 
    class on a given instance. However, this class is a container of 
    ``MacroMolecule`` and other ``Population`` instances first and 
    foremost. It delegates GA operations to its `ga_tools` attribute. 
    Any functionality related to the GA should be delegated to this 
    attribute. The ``gen_offspring`` and ``gen_mutants`` methods can 
    serve as a guide to how this should be done. A comphrehensive 
    account of how the interaction between these two classes is provided 
    in the developer's guide.
    
    For consistency and maintainability, collections of 
    ``MacroMolecule`` or ``Population`` instances should always be 
    placed in a ``Population`` instance. As a result, any function which 
    should return multiple ``MacroMolecule`` or ``Population`` instances 
    can be expected to return a single ``Population`` instance holding 
    the desired instances. Some functions may have to return the 
    population organized in a specific way, depending on use.
    
    The only operations directly addressed by this class and definined 
    within it are those relevant to its role as a container. It supports
    all expected and necessary container operations such as iteration, 
    indexing, membership checks (via the ``is in`` operator) as would be 
    expected. Additional operations such as comparison via the ``==``, 
    ``>``, etc. operators is also supported. Details of the various 
    implementations and a full list of supported operations can be found 
    by examining the included methods. Note that all comparison 
    operations are accounted for with the ``total_ordering`` decorator, 
    even if they are not explicity defined.

    Attributes
    ----------
    populations : list of ``Population`` instances
        A list of other instances of the ``Population`` class. This
        allows the implementation of subpopulations or evolutionary 
        islands. This attribute is also used for grouping macromolecules 
        within a given population for organizational purposes if need 
        be.
        
    members : list of ``MacroMolecule`` instances
        A list of ``MacroMolecule`` instances. These are the members of 
        the population which are not held within any subpopulations. 
        This means that not all members of a population are stored here. 
        To access all members of a population the generator method 
        ``all_members`` should be used.
    
    ga_tools : GATools, optional
        An instance of the ``GATools`` class. Calls to preform GA
        operations on the ``Population`` instance are delegated 
        to this attribute.
    
    """
    
    def __init__(self, *args):
        """
        Initializer for ``Population`` class.
        
        This initializer creates a new population from the 
        ``Population`` and ``MacroMolecule`` instances given as 
        arguments. It also accepts a single, optional, ``GATools`` 
        instance if the population is to have GA operations performed on 
        it.
        
        The arguments can be provided in any order regardless of type.
        
        Parameters
        ----------
        *args : MacroMolecule, Population, GATools
            A population is initialized with as many ``MacroMolecule`` 
            or ``Population`` arguments as required. These are placed 
            into the `members` or `populations` attributes, 
            respectively. A single ``GATools`` instance may be included 
            and will be placed into the `ga_tools` attribute. 
        
        Raises
        ------
        TypeError
            If the instance is initialized with something other than
            ``MacroMolecule``, ``Population`` or more than 1 ``GATools`` 
            object. 
            
        """    
        
        # Generate `populations`, `members` and `ga_tools` attributes.
        self.populations = []
        self.members = []
        self.ga_tools = None
    
        # Determine type of supplied arguments and place in the
        # appropriate attribute.  ``Population`` types added to
        # `populations` attribute, ``MacroMolecule`` into `members` and 
        # if ``GATools`` is supplied it is placed into `ga_tools`.
        # Raise a ``TypeError`` if more than 1 ``GATools`` argument
        # was supplied or if an argument was not ``GATools``, 
        # ``MacroMolecule`` or ``Population`` type.
        for arg in args:
            if isinstance(arg, Population):
                self.populations.append(arg)
                continue
            
            if isinstance(arg, MacroMolecule):
                self.members.append(arg)
                continue           
            
            if isinstance(arg, GATools) and self.ga_tools is None:
                self.ga_tools = arg
                continue

            # Some methods create a new ``Population`` instance by using 
            # another as a template. In these cases the `ga_tools` 
            # attribute of the template population is passed to the 
            # initializer. If the template instance did not have a 
            # defined ``GATools`` instance, the default ``None`` 
            # argument is passed. The following 2 lines prevent 
            # exceptions being raised in such a scenario. A consequence
            # of this is that any number of ``None`` arguments can be 
            # passed to the initializer. 
            if arg is None:
                continue

            raise TypeError(("Population can only be"
                             " initialized with 'Population',"
                             " 'MacroMolecule' and 1 'GATools' type."), 
                                                                    arg)
                    
    @classmethod
    def init_random_cages(cls, bb_db, lk_db, 
                          topologies, size, ga_tools):
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
            
        topolgies : iterable of ``Topology`` child classes
            An iterable holding topologies which should be randomly 
            selected for cage initialization.            
            
        size : int
            The size of the population to be initialized.
        
        ga_tools : GATools
            The GATools instance to be used by created population.
        
        """
        
        cage_gen = iter(Cage.init_random(bb_db, lk_db, topologies,
                    os.path.join(os.getcwd(),"init_{}.mol".format(x))) 
                        for x in range(size))
        
        return cls(*cage_gen, ga_tools)

    @classmethod
    def init_fixed_bb_cages(cls, bb_file, lk_db, 
                            topologies, size, ga_tools):
        
        cage_gen = iter(Cage.init_fixed_bb(bb_file, lk_db, topologies,
                    os.path.join(os.getcwd(),"init_{}.mol".format(x))) 
                        for x in range(size))
        return cls(*cage_gen, ga_tools)

    @classmethod
    def load(cls, file_name):
        """
        Initializes a Population from one dumped to a file with pickle.
        
        Parameters
        ----------
        file_name : str
            The full path of the file holding the dumped population.
            
        Returns
        -------
        Population
            The population stored in the dump file.
            
        """
        
        with open(file_name, 'rb') as dump_file:
            return pickle.load(dump_file)

    def all_members(self):
        """
        Yields all members in the population and its subpopulations.

        Yields
        ------
        MacroMolecule
            The next ``MacroMolecule`` instance held within the 
            population or its subpopulations.
        
        """
        
        # Go through `members` attribute and yield ``MacroMolecule`` 
        # instances held within one by one.
        for ind in self.members:
            yield ind
        
        # Go thorugh `populations` attribute and for each ``Population``
        # instance within, yield ``MacroMolecule`` instances from its
        # `all_members` generator.
        for pop in self.populations:
            yield from pop.all_members()
            
    def add_members(self, population, duplicates=False):
        """
        Adds ``MacroMolecule`` instances into `members`.        
        
        The ``MacroMolecule`` instances held within the supplied 
        ``Population`` instance, `population`, are added into the 
        `members` attribute of `self`. The supplied `population` itself 
        is not added. This means that any information the `population` 
        instance had about subpopulations is lost. This is because all 
        of its ``MacroMolecule`` instances are added into the `members` 
        attribute, regardless of which subpopulation they were 
        originally in.

        The `duplicates` parameter indicates whether multiple instances
        of the same macromolecule are allowed to be added into the 
        population. Note that the sameness of a macromolecule is judged 
        by the `same` method of the ``MacroMolecule`` class, which is 
        invoked by the ``in`` operator within this method. See the 
        `__contains__` method of the ``Population`` class for details on 
        how the ``in`` operator uses the `same` method.
        
        Parameters
        ----------
        population : Population (or iterable of ``MacroMolecule``s)
            ``MacroMolecule`` instances to be added to the `members` 
            attribute and/or ``Population`` instances who's members, as 
            generated by `all_members`, will be added to the `members` 
            attribute.
        
        duplicates : bool (default = False)
            When ``False`` only macromolecules which are not already 
            held by the population will be added. ``True`` allows more 
            than one instance of the same macromolecule to be added. 
            Whether two macromolecules are the same is defined by the 
            `same` method of the ``MacroMolecule`` class.
        
        Modifies
        --------
        members
            Adds instances into the `members` attribute of `self`.        
        
        Returns
        -------
        None : NoneType        
        
        """
        
        if duplicates:
            self.members.extend(mol for mol in population)
        else:
            self.members.extend(mol for mol in population
                                                    if mol not in self)
    def add_subpopulation(self, population):
        """
        Appends a population into the `populations` attribute.        
        
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
        
        self.populations.append(population)
        
    def remove_duplicates(self, between_subpops=True, top_seen=None):
        """
        Removes duplicates from a population while preserving structure.        
        
        The question of which ``MacroMolecule`` instance is preserved 
        from a choice of two is difficult to answer. The iteration 
        through a population is depth-first, so a rule such as ``the 
        macromolecule in the topmost population is preserved`` is not 
        the case here. Rather, the first ``MacroMolecule`` instance 
        iterated through is preserved.
        
        However, this question is only relevant if duplicates in 
        different subpopulations are being removed. In this case it is 
        assumed that it is more important to have a single instance than
        to worry about which subpopulation it is in.
        
        If the duplicates are being removed from within subpopulations,
        each subpopulation will end up with a single instance of all
        macromolecules held before. There is no ``choice``.
        
        Parameters
        ----------
        between_subpops : bool (default = False)
            When ``False`` duplicates are only removed from within a
            given subpopulation. If ``True`` all duplicates are removed,
            regardless of which subpopulation they are in.
        
        Modifies
        --------
        members
            Duplicate instances are removed from the `members` attribute
            of the population or subpopulations.
        
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
            if type(top_seen) == set:
                seen = top_seen
                
            self.members = list(dedupe(self.members, seen=seen))
            for subpop in self.populations:
                subpop.remove_duplicates(between_subpops, top_seen=seen)
        
        # If duplicates are only removed from within the same
        # subpopulation, only the `members` attribute of each
        # subpopulation needs to be cleared of duplicates. To do this,
        # each `members` attribute is deduped recursively.
        if not between_subpops:
            self.members = list(dedupe(self.members))
            for subpop in self.populations:
                subpop.remove_duplicates(between_subpops=False)
        
    def select(self, type_='generational'):
        """
        Returns a generator field yielding selected members of `self`.
        
        Selection is a GA procedure and as a result this method merely 
        delegates the selection request to the ``Selection`` instance 
        held within the `ga_tools` attribute. The ``Selection`` instance
        then returns a generator which yields ``MacroMolecule``
        instances held within the population. Which macromolecules are
        yielded depends on the selection algorithm which was chosen
        during initialization and when calling this method. The 
        selection instance (`self.ga_tools.selection`) returns the 
        generator when it is called. See ``Selection`` class 
        documentation for more information.
        
        Because selection is required in a number of different ways,
        such as selecting the parents, ``MacroMolecule`` instances for 
        mutation and ``MacroMolecule`` instances for the next 
        generation, the type of selection must be specificed with the 
        `type_` parameter. The valid values for `type_` will correspond 
        to one of the attribute names of the ``Selection`` instance.

        For example, if `type_` is set to 'mating' a selection 
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
            ``Selection`` class. Check ``Selection`` class documentation
            for details. 
            
            Valid values include:              
                'generational' - selects the next generation
                'mating' - selects parents
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
        
    def gen_offspring(self):
        """
        Returns a population of offspring ``MacroMolecule`` instances.        
        
        This is a GA operation and as a result this method merely
        delegates the request to the ``Mating`` instance held in the 
        `ga_tools` attribute. The ``Mating`` instance takes care of 
        selecting parents and combining them to form offspring. The
        ``Mating`` instace delegates the selection to the ``Selection`` 
        instance as would be expected. The request to perform mating is 
        done by calling the ``Mating`` instance with the population as 
        the argument. Calling of the ``Mating``instance returns a 
        ``Population`` instance holding the generated offspring. All 
        details regarding the mating procedure are handled by the 
        ``Mating`` instance.

        For more details about how mating is implemented see the
        ``Mating`` class documentation.
        
        Returns
        -------
        Population
            A population holding offspring created by mating contained 
            the ``MacroMolecule`` instances.
        
        """
        
        return self.ga_tools.mating(self)
        
    def gen_mutants(self):
        """
        Returns a population of mutant ``MacroMolecule`` instances.        

        This is a GA operation and as a result this method merely
        delegates the request to the ``Mutation`` instance held in the 
        `ga_tools` attribute. 

        Returns
        -------
        Population
            A population holding mutants created by mutating contained
            ``MacroMolecule`` instances.
        
        """

        return self.ga_tools.mutation(self)
        
    def write_population_to_dir(self, dir_path):
        """
        Writes the ``.mol`` files of members to a directory.
        
        This writes the pristine version of the ``.mol`` file.        
        
        Parameters
        ----------
        dir_path : str
            The full path of the directory into which the ``.mol`` file
            is copied.
        
        Returns
        -------
        None : NoneType
        
        """
        
        # If the directory does not exist, create it.        
        if not os.path.exists(dir_path):
            os.mkdir(dir_path)
        
        for member in self:
            member.write_mol_file('prist', dir_path)

    def optimize_population(self):
        """
        Optimizes all the members of the population.

        This function should invoke either the ``optimize_all`` or
        ``optimize_all_serial`` functions. ``optimize_all`` optimizes
        all members of the population in parallel. 
        ``optimize_all_serial`` does them serially. Probably best not to
        use the serial version unless debugging.
        
        The parallel optimization creates cloned instances of the
        population's members. It is these that are optimized. This means
        that the ``.mol`` files are changed but any instance attributes
        are not. See ``optimize_all`` function documentation in 
        ``optimization.py`` for more details.

        Modifies
        --------
        MacroMolecule
            This function replaces the pristine rdkit molecule instances 
            with optimizes versions. It also replaces the content of the
            pristine ``.mol`` files with pristine structures.

        Returns
        -------
        iterator of MacroMolecule objects
            If a parallel optimization was chosen, this iterator yields 
            the ``MacroMolecule`` objects that have had their attributes 
            changed as a result of the optimization. They are modified 
            clones of the original population's macromolecules.       
            
            If a serial optimization is done the iterator does not yield
            clones.
        
        """

        return optimize_all(self.ga_tools.optimization, self)
            
    def calculate_member_fitness(self):
        calc_fitness(self.ga_tools.fitness, self)     
     
    def mean(self, attr_name, delegate=False):
        """
        Calculates the mean value of some attribute of the members.
        
        Parameters
        ----------
        attr_name : str
            The name of the attribute whose mean value should be
            calculated.
            
        delegate : False or string (default = False)
            If ``False`` the attr_name is assumed to correspond to an
            attribute held by the members directly. If a string is
            provided it should hold the name of the attribute which
            holds the desired attribute. For example, if `cavity_size`
            is requested: attr_name=`cavity_size`, delegate=`topology`.
            This is because `cavity_size` is not held by members
            directly. It is held in the `topology` attribute of the 
            members.

        Returns
        -------
        float
            The mean of a given attribute held by members of the 
            population.

        """
        
        if delegate:
            return mean(getattr(getattr(x, delegate), attr_name) 
                                                          for x in self)    
        else:
            return mean(getattr(x, attr_name) for x in self)

    def dump(self, file_name):
        """
        Write the population object to a file.
        
        Parameters
        ----------
        file_name : str
            The full path of the file to which the population should
            be written.
            
        Returns
        -------
        None : NoneType

        """
        
        with open(file_name, 'wb') as dump_file:    
            pickle.dump(self, dump_file)        
        
     
    def __iter__(self):
        """
        Allows the use of ``for`` loops, ``*`` and ``iter`` function.        

        When ``Population`` instances are iterated through they yield 
        ``MacroMolecule`` instances generated by the `all_members` 
        method. It also means that a ``Population`` instance can be 
        unpacked with the ``*`` operator. This will produce the 
        ``MacroMolecule`` instances yielded by the `all_members` method.

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

        Macromolecules held by the ``Population`` instance can be 
        accesed by their index. Slices are also supported. These return 
        a new ``Population`` instance holding the ``MacroMolecule`` 
        instances with the requested indices. Using slices will return 
        a flat ``Population`` instance meaing no subpopulation
        information is preserved. All of the ``MacroMolecule`` instances 
        are placed into the `members` attribute of the returned 
        ``Population`` instance.

        The index corresponds to the ``MacroMolecule`` yielded by the 
        `all_members` method.

        This can be exploited if one desired to remove all
        subpopulations and transfer all the ``MacroMolecules`` instances 
        into the members attribute. For example, 
        
        >>> pop2 = pop[:]
        
        ``pop2`` is a ``Population`` instance with all the same
        ``MacroMolecule`` instances as ``pop``, however all 
        ``MacroMolecule`` instances are held within its `members` 
        attribute and its `populations` attribute is empty. This may or 
        may not be the case for the ``pop`` instance.   
        
        Parameters
        ----------
        key : int, slice
            An int or slice can be used depending on if a single 
            ``MacroMolecule`` instance needs to be returned or a 
            collection of ``MacroMolecule`` instances.
        
        Returns
        -------
        MacroMolecule
            If the supplied `key` is an ``int``. Returns the 
            ``MacroMolecule`` instance with the corresponding index from 
            the `all_members` generator.
        
        Population
            If the supplied `key` is a ``slice``. The returned 
            ``Population`` instance holds ``MacroMolecule`` instances in 
            its `members` attribute. The ``MacroMolecule`` instances 
            correspond to indices defined by the slice. The slice is 
            implemented on the `all_members` generator.
        
        Raises
        ------
        TypeError
            If the supplied `key` is not an ``int`` or ``slice`` type.        
        
        """
        
        # Determine if provided key was an ``int`` or a ``slice``.
        # If ``int``, return the corresponding ``MacroMolecule`` 
        # instance from the `all_members` generator.
        if isinstance(key, int):
            return list(self.all_members())[key]
        
        # If ``slice`` return a ``Population`` of the corresponding 
        # ``MacroMolecule`` instances. The returned ``Population`` will 
        # have the same `ga_tools` attribute as original ``Population`` 
        # instance.        
        if isinstance(key, slice):
            mols = itertools.islice(self.all_members(), 
                                     key.start, key.stop, key.step)
            return Population(*mols, self.ga_tools)

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
        
        returns a new population, pop3. The returned population contains 
        all the ``MacroMolecule`` instances in pop1 except those also in 
        pop2. This refers to all of the ``MacroMolecule`` instances, 
        including those held within any subpopulations. The returned 
        population is flat. This means all information about 
        subpopulations in pop1 is lost as all the ``MacroMolecule`` 
        instances are held in the `members` attribute of pop3.

        The resulting population, pop3, will inherit the `ga_tools` 
        attribute from pop1.

        Parameters
        ----------
        other : Population
            A collection of ``MacroMolecule`` instances to be removed 
            from `self`, if held by it.
            
        Returns
        -------
        Population
            A flat population of ``MacroMolecule`` instances which are 
            not also held in `other`.

        """

        new_pop = Population(self.ga_tools)        
        new_pop.add_members(mol for mol in self 
                                                if mol not in other)
        return new_pop
        
    def __add__(self, other):
        """
        Allows use fo the ``+`` operator.
        
        Creates a new ``Population`` instance which holds two 
        subpopulations and no direct members. The two subpopulations
        are the two ``Population`` instances on which the ``+`` operator
        was applied.
        
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
        item : MacroMolecule

        Returns
        -------
        bool
        
        """

        return any(item.same(mol) for mol in self.all_members())

    def __str__(self):
        output_string = (" Population " + str(id(self)) + "\n" + 
                            "--------------------------\n" + 
                            "\tMembers\n" + "   ---------\n")
        
        for mol in self.members:
            output_string += "\t"  + str(mol) + "\n"
        
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

    """
    The following methods are inteded for convenience while 
    debugging or testing and should not be used during typical 
    execution of the program.
    
    """

    @classmethod
    def init_empty(cls):
        pops = []
        
        for x in range(0,8):
            pop = cls(*iter(Cage.init_empty() for x in range(0,3)), 
                      GATools.default())
            pops.append(pop)
        
        pops[1].populations.append(pops[3])
        pops[1].populations.append(pops[4])
        pops[1].populations.append(pops[5])
        
        pops[2].populations.append(pops[6])
        pops[2].populations.append(pops[7])
        
        pops[0].populations.append(pops[1])
        pops[0].populations.append(pops[2])        
        
        return pops[0]