"""
Optimizer
=========

.. toctree::
    :maxdepth: 2

    Collapser <\
stk.molecular.topology_graphs.topology_graph.optimizers.collapser\
>

    Periodic Collapser <\
stk.molecular.topology_graphs.topology_graph.optimizers.\
periodic_collapser\
>

    MCHammer <\
stk.molecular.topology_graphs.topology_graph.optimizers.mchammer\
>
    Spinner <\
stk.molecular.topology_graphs.topology_graph.optimizers.spinner\
>

    NullOptimizer <\
stk.molecular.topology_graphs.topology_graph.optimizers.null\
>

"""

from ....construction_state import ConstructionState


class Optimizer:
    """
    An abstract base class for optimizers.

    An optimizer is used to change the structure of the molecule under
    construction to be more realistic.

    """

    def optimize(
        self,
        state: ConstructionState,
    ) -> ConstructionState:
        """
        Optimize the structure of a molecule under construction.

        Parameters:

            state:
                The molecule being constructed.

        Returns:

            The optimized construction state.

        """

        raise NotImplementedError()
