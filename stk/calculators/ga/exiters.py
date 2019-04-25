"""
Defines :class:`Exiter` classes.

:class:`Exiter` objects check to see if an exit condition for the GA
has been fulfilled.

.. _`adding exiters`:

Extending stk: Making new exiters.
----------------------------------

A new :class:`Exiter` class should inherit :class:`Exiter` and
define an :meth:`~Exiter.exit` method, which takes progress population
as its only argument and returns ``True`` or ``False`` depending on if
an exit condition has been satisfied. The progress population is a
:class:`.Population` instance which holds each GA generation as a
subpopulation.

"""


class Exiter:
    """
    Checks if the exit criterion for the GA has been satisfied.

    """

    def exit(self, progress):
        """
        Checks to see the the GA should stop.

        Parameters
        ----------
        progress : :class:`.Population`
            A population where every generation is a
            subpopulation.

        Returns
        -------
        :class:`bool`
            ``True`` if the GA should stop and ``False`` otherwise.

        """

        raise NotImplementedError()


class AnyExiter(Exiter):
    """
    Checks if any :class:`Exiter` has satisfied its exit condition.

    Attributes
    ----------
    exiters : :class:`tuple` of :class:`Exiter`
        :class:`Exiter` objects which are checked to see if their exit
        conditions have been satisfied.

    """

    def __init__(self, *exiters):
        """
        Initializes a :class:`AnyExiter` instance.

        Parameters
        ----------
        *exiters : :class:`Exiter`
            :class:`Exiter` objects which are checked to see if their
            exit conditions have been satisfied.

        """

        self.exiters = exiters

    def exit(self, progress):
        """
        Checks to see if any exit condition has been satisfied.

        Parameters
        ----------
        progress : :class:`.Population`
            A population where every generation is a subpopulation.

        Returns
        -------
        :class:`bool`
            ``True`` if any :class:`Exiter` in :attr:`exiters` has
            satisfied its exit condition.

        """

        return any(exiter.exit(progress) for exiter in self.exiters)


class AllExiters(Exiter):
    """
    Checks if all :class:`Exiter` objects return ``True``.

    Attributes
    ----------
    exiters : :class:`tuple` of :class:`Exiter`
        :class:`Exiter` objects which are checked to see if their exit
        conditions have been satisfied.

    """

    def __init__(self, *exiters):
        """
        Initializes a :class:`AllExiter` instance.

        Parameters
        ----------
        *exiters : :class:`Exiter`
            :class:`Exiter` objects which are checked to see if their
            exit conditions have been satisfied.

        """

        self.exiters = exiters

    def exit(self, progress):
        """
        Checks to see if all exit conditions have been satisfied.

        Parameters
        ----------
        progress : :class:`.Population`
            A population where every generation is a subpopulation.

        Returns
        -------
        :class:`bool`
            ``True`` if all :class:`Exiter` objects in :attr:`exiters`
            have satisfied its exit condition.

        """

        return all(exiter.exit(progress) for exiter in self.exiters)


class NumGenerations(Exiter):
    """
    Stop the GA after a certain number of generations.

    Attributes
    ----------
    num_generations : :class:`int`
        The number of generations after which the GA should stop.

    """

    def __init__(self, num_generations):
        """
        Initializes a :class:`NumGenerations` instance.

        Parameters
        ----------
        num_generations : :class:`int`
            The number of generations after which the GA should stop.

        """

        self.num_generations = num_generations

    def exit(self, progress):
        """
        Checks if a number of generations has passed.

        Parameters
        ----------
        progress : :class:`.Population`
            A population where every generation is a subpopulation.

        Returns
        -------
        ``True`` if :attr:`num_generations` or more has passed.

        """

        return len(progress.populations) >= self.num_generations


class MolPresent(Exiter):
    """
    Stops the GA if a specific molecule has been found.

    Attributes
    ----------
    mol : :class:`.Molecule`
        A molecule which if present in any of the GA's generations
        causes it to stop running.

    """

    def __init__(self, mol):
        """
        Initializes a :class:`MolPresent` instance.

        Parameters
        ----------
        mol : :class:`.Molecule`
            A molecule which if present in any of the GA's generations
            causes it to stop running.

        """

        self.mol = mol

    def exit(self, progress):
        """
        ``True`` if :attr:`mol` is in `progress`.

        Parameters
        ----------
        progress : :class:`.Population`
            A population where every generation is a subpopulation.

        Returns
        -------
        :class:`bool`
            ``True`` if :attr:`mol` in `progress`, ``False`` otherwise.

        """

        # Check for the presence of the molecule, starting with the
        # newest generation first.
        for pop in progress.populations:
            if any(mol.same(self.mol) for mol in pop):
                return True
        return False


class FitnessPlateau(Exiter):
    """
    Checks if the fittest molecules remain the same.

    Attributes
    ----------
    num_generations : :class:`int`
        Number of generations in which the fittest molecules did
        not change.

    top_members : :class:`int`
        The number of fittest molecules which are checked. This
        number needs to be smaller than the population size.

    """

    def __init__(self, num_generations, top_members=1):
        """
        Initializes a :class:`FitnessPlateau` instance.

        Parameters
        ----------
        num_generations : :class:`int`
            Number of generations in which the fittest molecules did
            not change.

        top_members : :class:`int`, optional
            The number of fittest molecules which are checked. This
            number needs to be smaller than the population size.

        """

        self.num_generations = num_generations
        self.top_members = top_members

    def exit(self, progress):
        """
        Checks if the fittest molecules changed between generations.

        Parameters
        ----------
        progress : :class:`.Population`
            A population where every generation is a subpopulation.

        Returns
        -------
        :class:`bool`
            ``True`` if the :attr:`top_members` have not changed over
            the last :attr:`num_gens` generations.

        """

        # Check that the GA has run for more than num_gens generations.
        if len(progress.populations) >= self.num_generations:
            gens = set()
            for i in range(self.num_generations):
                gen = sorted(progress.populations[-i-1],
                             reverse=True,
                             key=lambda mol: mol.fitness)
                # Get the top members of the generation.
                keys = frozenset(
                    mol.key for mol in gen[:self.top_members]
                )
                gens.add(keys)
            unique_gens = len(gens)
            if unique_gens == 1:
                return True
            else:
                return False

        return False
