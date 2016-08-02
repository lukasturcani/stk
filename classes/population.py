import itertools



from .molecular import MacroMolecule, Cage
from .ga import GATools
from ..convenience_functions import dedupe, flatten

class Population:
    """
    A container for instances of ``Cage`` and ``Population``.

    This is the central class of MMEA. The GA is invoked by calling the
    ``gen_offspring``, ``gen_mutants`` and ``select`` methods of this 
    class on a given instance. However, this class is a container of 
    ``Cage`` and other ``Population`` instances first and foremost. It 
    delegates GA operations to its `ga_tools` attribute. Any 
    functionality related to the GA should be delegated to this 
    attribute. The ``gen_offspring`` and ``gen_mutants`` methods can 
    serve as a guide to how this should be done. A comphrehensive 
    account of how the interaction between these two classes is provided 
    in the developer's guide.
    
    For consistency and maintainability, collections of ``Cage`` or 
    ``Population`` instances should always be placed in a ``Population`` 
    instance. As a result, any function which should return multiple 
    ``Cage`` or ``Population`` instances can be expected to return a 
    single ``Population`` instance holding the desired instances. Some 
    functions will have to return the population organized in a specific 
    way. For example, functions generating the parent pool will generate 
    a population with no direct members but multiple subpopulations of 2 
    members each. More specific guidelines are provided within the 
    ``Mating`` class.
    
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
    populations : list
        A list of other instances of the ``Population`` class. This
        allows the implementation of subpopulations or 'evolutionary 
        islands'. This attribute is also used for grouping cages within
        a given population such as when grouping parents together in a 
        parent pool.
        
    members : list
        A list of ``Cage`` instances. These are the members of the
        population which are not held within any subpopulations. This 
        means that not all members of a population are stored here. 
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
        ``Population`` and ``Cage`` instances given as arguments.
        It also accepts a single, optional, ``GATools`` instance if the 
        population is to have GA operations performed on it.
        
        The arguments can be provided in any order regardless of type.
        
        Parameters
        ----------
        *args : Cage, Population, GATools
            A population is initialized with as many ``Cage`` or 
            ``Population`` arguments as required. These are placed into 
            the `members` or `populations` attributes, respectively. 
            A single ``GATools`` instance may be included and will be
            placed into the `ga_tools` attribute. 
        
        Raises
        ------
        TypeError
            If the instance is initialized with something other than
            ``Cage``, ``Population`` or more than 1 ``GATools`` object. 
            
        """    
        
        # Generate `populations`, `members` and `ga_tools` attributes.
        self.populations = []
        self.members = []
        self.ga_tools = None
    
        # Determine type of supplied arguments and place in the
        # appropriate attribute.  ``Population`` types added to
        # `populations` attribute, ``Cage`` into `members` and if
        # ``GATools`` is supplied it is placed into `ga_tools`.
        # Raise a ``TypeError`` if more than 1 ``GATools`` argument
        # was supplied or if an argument was not ``GATools``, ``Cage`` 
        # or ``Population`` type.
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
                             " 'Cage' and 1 'GATools' type."), arg)
                                    
    def all_members(self):
        """
        Yields all members in the population and its subpopulations.

        Yields
        ------
        Cage
            The next ``Cage`` instance held within the population or its
            subpopulations.
        
        """
        
        # Go through `members` attribute and yield ``Cage`` instances
        # held within one by one.
        for ind in self.members:
            yield ind
        
        # Go thorugh `populations` attribute and for each ``Population``
        # instance within, yield ``Cage`` instances from its
        # `all_members()` generator.
        for pop in self.populations:
            yield from pop.all_members()
            
    def add_members(self, population, duplicates=False):
        """
        Adds ``Cage`` instances into `members`.        
        
        The ``Cage`` instances held within the supplied ``Population`` 
        instance, `population`, are added into the `members` attribute 
        of `self`. The supplied `population` itself is not added. This 
        means that any information the `population` instance had about 
        subpopulations is lost. This is because all of its ``Cage`` 
        instances are added into the `members` attribute, regardless of 
        which subpopulation they were originally in.

        The `duplicates` parameter indicates whether multiple instances
        of the same cage are allowed to be added into the population.
        Note that the sameness of a cage is judged by the `same_cage`
        method of the ``Cage`` class, which is invoked by the ``in``
        operator within this method. See the `__contains__` method of 
        the ``Population`` class for details on how the ``in`` operator 
        uses the `same_cage` method.
        
        Parameters
        ----------
        population : Population
            ``Cage`` instances to be added to the `members` attribute
            and/or ``Population`` instances who's members, as generated 
            by `all_members`, will be added to the `members` attribute.
        
        duplicates : bool (default = False)
            When `False` only cages which are not already held by
            the population will be added. `True` allows more than one
            instance of the same cage to be added. Whether two cages
            are the same is defined by the `same_cage` method of the
            ``Cage`` class.
        
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
            The population to added as a subpopulation.
            
        Returns
        -------
        None : NoneType
        
        """
        
        self.populations.append(population)
        
    def remove_duplicates(self, between_subpops=True, top_seen=None):
        """
        Removes duplicates from a population while preserving structure.        
        
        The question of which ``Cage`` instance is preserved from a 
        choice of two is difficult to answer. The iteration through a 
        population is depth-first so a rule such as ``the cage in the
        topmost population is preserved`` is not the case here. Rather,
        the first ``Cage`` instance iterated through is preserved.
        
        However, this question is only relevant if duplicates in 
        different subpopulations are being removed. In this case it is 
        assumed that it is more important to have a single instance than
        to worry about which subpopulation it is in.
        
        If the duplicates are being removed from within subpopulations,
        each subpopulation will end up with a single instance of all
        cages held before. There is no ``choice``.
        
        Parameters
        ----------
        between_subpops : bool (default = False)
            When ``False`` duplicates are only removed from within a
            given subpopulation. If ``True`` all duplicates are removed,
            regardless of which subpopulation they are in.
            
        Returns
        -------
        None : NoneType
        
        """
        
        if between_subpops:
            if top_seen is None:
                seen = set()
            if type(top_seen) == set:
                seen = top_seen
                
            self.members = list(dedupe(self.members, seen=seen))
            for subpop in self.populations:
                subpop.remove_duplicates(between_subpops, top_seen=seen)
        
        if not between_subpops:
            self.members = list(dedupe(self.members))
            for subpop in self.populations:
                subpop.remove_duplicates(between_subpops=False)
        
    def select(self, type_='generational'):
        """
        Selects some members to form a new ``Population`` instance.
        
        Selection is a GA procedure and as a result this method merely 
        delegates the selection request to the ``Selection`` instance 
        held within the `ga_tools` attribute. The ``Selection`` instance
        then returns the new ``Population`` instance. This 
        ``Population`` instance is then returned by this method.
        The selection instance (`self.ga_tools.selection`) performs the
        selection by being called. See ``Selection`` class documentation
        for more information.
        
        Because selection is required in a number of different ways,
        such as selecting the parents, ``Cage`` instances for mutation
        and ``Cage`` instances for the next generation, the type of 
        selection must be specificed with the `type_` parameter. The
        valid values for `type_` will correspond to one of the
        attribute names of the ``Selection`` instance.

        For example, if `type_` is set to 'mating' a selection 
        algorithm which creates a parent pool will be invoked. If the 
        `type_` is set to 'generational' an algorithm which selects the 
        next generation will be invoked. It should be noted that a 
        ``Population`` instance representing a parent pool will be
        organized differently to one representing a generation. See
        ``Selection`` class documentation for more details.
        
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
                'mutation' - selects ``Cage`` instances for mutation

        Returns
        -------        
        Population
            A population generated by using applying a selection
            algorithm. Can represent a generation, a parent pool, etc.
        
        """
        
        return self.ga_tools.selection(self, type_)
        
    def gen_offspring(self):
        """
        Returns a population of offspring ``Cage`` instances.        
        
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
            the ``Cage`` instances.
        
        """
        
        return self.ga_tools.mating(self)
        
    def gen_mutants(self):
        """
        Returns a population of mutant ``Cage`` instances.        
        
        
        
        Returns
        -------
        Population
            A population holding mutants created by mutating contained
            ``Cage`` instances.
        
        """

        return self.ga_tools.mutation(self)
        
    def __iter__(self):
        """
        Allows the use of ``for`` loops, ``*`` and ``iter`` function.        

        When ``Population`` instances are iterated through they yield 
        ``Cage`` instances generated by the `all_members` method. It 
        also means that a ``Population`` instance can be unpacked with
        the ``*`` operator. This will produce the ``Cage`` instances
        yielded by the `all_members` method.

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

        Cages held by the ``Population`` instance can be accesed by
        their index. Slices are also supported. These return a new
        ``Population`` instance holding the ``Cage`` instances with
        the requested indices. Using slices will return a flat 
        ``Population`` instance meaing no subpopulation
        information is preserved. All of the ``Cage`` instances are
        placed into the `members` attribute of the returned 
        ``Population`` instance.

        The index corresponds to the ``Cages`` yielded by the 
        `all_members` method.

        This can be exploited if one desired to remove all
        subpopulations and transfer all the ``Cage`` instances into the 
        members attribute. For example, 
        
        >>> pop2 = pop[:]
        
        ``pop2`` is a ``Population`` instance with all the same
        ``Cage`` instances as ``pop``, however all ``Cage`` 
        instances are held within its `members` attribute and its 
        `populations` attribute is empty. This may or may not be the 
        case for the ``pop`` instance.   
        
        Parameters
        ----------
        key : int, slice
            An int or slice can be used depending on if a single 
            ``Cage`` instance needs to be returned or a collection of 
            ``Cage`` instances.
        
        Returns
        -------
        Cage
            If the supplied `key` is an ``int``. Returns the ``Cage`` 
            instance with the corresponding index from the `all_members` 
            generator.
        
        Population
            If the supplied `key` is a ``slice``. The returned 
            ``Population`` instance holds ``Cage`` instances in its 
            `members` attribute. The ``Cage`` instances correspond to 
            indices defined by the slice. The slice is implemented 
            on the `all_members` generator.
        
        Raises
        ------
        TypeError
            If the supplied `key` is not an ``int`` or ``slice`` type.        
        
        """
        
        # Determine if provided key was an ``int`` or a ``slice``.
        # If ``int``, return the corresponding ``Cage`` instance from
        # the `all_members` generator.
        if type(key) is int:
            return list(self.all_members())[key]
        
        # If ``slice`` return a ``Population`` of the corresponding 
        # ``Cage`` instances. The returned ``Population`` will have the 
        # same `ga_tools` attribute as original ``Population`` instance.        
        if type(key) is slice:
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
        all the ``Cage`` instances in pop1 except those also in pop2.
        This refers to all of the ``Cage`` instances, including those
        held within any subpopulations. The returned population is 
        flat. This means all information about subpopulations in pop1 is 
        lost as all the ``Cage`` instances are held in the `members` 
        attribute of pop3.

        The resulting population, pop3, will inherit the `ga_tools` 
        attribute from pop1.

        Parameters
        ----------
        other : Population
            A collection of ``Cage`` instances to be removed from 
            `self`, if held by it.
            
        Returns
        -------
        Population
            A flat population of ``Cage`` instances which are not also
            held in `other`.

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
        item : Cage

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