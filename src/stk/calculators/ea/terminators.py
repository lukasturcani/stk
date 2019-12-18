"""
Termination
===========


#. :class:`.NumGenerations`
#. :class:`.MoleculePresent`
#. :class:`.FitnessPlateau`
#. :class:`.AnyTerminator`
#. :class:`.AllTerminators`



:class:`Terminator` objects check to see if an exit condition for the
EA has been fulfilled.

.. _`adding Terminators`:

Making New Terminators
----------------------

A new :class:`Terminator` class should inherit :class:`Terminator` and
define an :meth:`~Terminator.terminate` method, which takes a progress
population as its only argument and returns ``True`` or ``False``
depending on if an exit condition has been satisfied. The progress
population is a :class:`.Population` instance which holds each EA
generation as a subpopulation.

"""


import rdkit.Chem.AllChem as rdkit
from functools import partial


class Terminator:
    """
    Checks if the exit criterion for the EA has been satisfied.

    """

    def terminate(self, progress):
        """
        Check to see the the EA should stop.

        Parameters
        ----------
        progress : :class:`.Population`
            A population where every generation is a
            subpopulation.

        Returns
        -------
        :class:`bool`
            ``True`` if the EA should stop and ``False`` otherwise.

        """

        raise NotImplementedError()


class AnyTerminator(Terminator):
    """
    Checks if any :class:`Terminator` has satisfied its exit condition.

    """

    def __init__(self, *Terminators):
        """
        Initialize a :class:`AnyTerminator` instance.

        Parameters
        ----------
        *Terminators : :class:`Terminator`
            :class:`Terminator` objects which are checked to see if their
            exit conditions have been satisfied.

        """

        self._Terminators = Terminators

    def terminate(self, progress):
        """
        Check to see if any exit condition has been satisfied.

        Parameters
        ----------
        progress : :class:`.Population`
            A population where every generation is a subpopulation.

        Returns
        -------
        :class:`bool`
            ``True`` if any :class:`Terminator` in :attr:`Terminators` has
            satisfied its exit condition.

        """

        return any(
            Terminator.terminate(progress)
            for Terminator in self._Terminators
        )


class AllTerminators(Terminator):
    """
    Checks if all :class:`Terminator` objects return ``True``.

    """

    def __init__(self, *Terminators):
        """
        Initialize a :class:`AllTerminator` instance.

        Parameters
        ----------
        *Terminators : :class:`Terminator`
            :class:`Terminator` objects which are checked to see if their
            exit conditions have been satisfied.

        """

        self._Terminators = Terminators

    def terminate(self, progress):
        """
        Checks to see if all exit conditions have been satisfied.

        Parameters
        ----------
        progress : :class:`.Population`
            A population where every generation is a subpopulation.

        Returns
        -------
        :class:`bool`
            ``True`` if all :class:`Terminator` objects in :attr:`Terminators`
            have satisfied its exit condition.

        """

        return all(
            Terminator.terminate(progress)
            for Terminator in self._Terminators
        )


class NumGenerations(Terminator):
    """
    Stop the EA after a certain number of generations.

    """

    def __init__(self, num_generations):
        """
        Initialize a :class:`NumGenerations` instance.

        Parameters
        ----------
        num_generations : :class:`int`
            The number of generations after which the EA should stop.

        """

        self._num_generations = num_generations

    def terminate(self, progress):
        """
        Check if a number of generations has passed.

        Parameters
        ----------
        progress : :class:`.Population`
            A population where every generation is a subpopulation.

        Returns
        -------
        ``True`` if :attr:`num_generations` or more has passed.

        """

        return len(progress.subpopulations) >= self._num_generations


class MoleculePresent(Terminator):
    """
    Stops the EA if a specific molecule has been found.

    """

    def __init__(self, mol):
        """
        Initialize a :class:`MoleculePresent` instance.

        Parameters
        ----------
        mol : :class:`.Molecule`
            A molecule which if present in any of the EA's generations
            causes it to stop running.

        """

        self._mol = mol
        self._is_same_molecule = partial(
            self._is_same_molecule,
            mol.to_rdkit_mol()
        )

    def terminate(self, progress):
        """
        Return ``True`` if :attr:`mol` is in `progress`.

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
        for pop in progress.subpopulations:
            if any(
                self._is_same_molecule(mol.to_rdkit_mol())
                for mol in pop
            ):
                return True
        return False

    @staticmethod
    def _is_same_molecule(mol1, mol2):
        return rdkit.MolToInchi(mol1) == rdkit.MolToInchi(mol2)


class FitnessPlateau(Terminator):
    """
    Checks if the fittest molecules remain the same.

    """

    def __init__(self, num_generations, top_members=1):
        """
        Initialize a :class:`FitnessPlateau` instance.

        Parameters
        ----------
        num_generations : :class:`int`
            Number of generations in which the fittest molecules did
            not change.

        top_members : :class:`int`, optional
            The number of fittest molecules which are checked. This
            number needs to be smaller than the population size.

        """

        self._num_generations = num_generations
        self._top_members = top_members

    def terminate(self, progress):
        """
        Check if the fittest molecules changed between generations.

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

        # Check that the EA has run for more than num_gens generations.
        if len(progress.subpopulations) >= self._num_generations:
            gens = set()
            for i in range(self._num_generations):
                gen = sorted(
                    progress.subpopulations[-i-1],
                    reverse=True,
                    key=lambda mol: mol.fitness
                )
                # Get the top members of the generation.
                keys = frozenset(
                    mol.get_identity_key() for mol in gen[:self._top_members]
                )
                gens.add(keys)
            unique_gens = len(gens)
            if unique_gens == 1:
                return True
            else:
                return False

        return False
