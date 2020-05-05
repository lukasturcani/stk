"""
Reaction Factory
================

.. toctree::
    :maxdepth: 2

    Generic Reaction Factory <\
stk.molecular.reactions.factories.generic_reaction_factory\
>
    Dative Reaction Factory <\
stk.molecular.reactions.factories.dative_reaction_factory\
>

"""


class ReactionFactory:
    """
    An abstract base class for reaction factories.

    Reaction factories are responsible for creating :class:`.Reaction`
    instances. Different subclasses of this class will provide
    different options for creating :class:`.Reaction` instances.

    If you want to change which reactions are used to create and
    delete atoms and bonds during :class:`.ConstructedMolecule`
    construction, you want to subclass this abstract base class
    and implement :meth:`.get_reaction`. Your implementation can
    then pick which reaction to use for a particular
    :class:`.EdgeGroup`. You will then pass an instance of your
    :class:`.ReactionFactory` subclass to the :class:`.TopologyGraph`
    initializer, so that it knows to use it.

    Examples
    --------
    *Subclass Implementation*

    The source code of the subclasses, listed in
    :mod:`.reaction_factory`, can serve as good examples.

    """

    def get_reaction(self, construction_state, edge_group):
        """
        Get a reaction to use on the `edge_group`.

        Parameters
        ----------
        construction_state : :class:`.ConstructionState`
            The state of the current construction.

        edge_group : :class:`.EdgeGroup`
            The edge group for which a reaction should be found.

        Returns
        -------
        :class:`.Reaction`
            The reaction to use on the `edge_group`.

        """

        raise NotImplementedError()
