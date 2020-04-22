"""
Host Guest Complex
==================

"""

from .vertices import _HostVertex, _GuestVertex
from ..topology_graph import TopologyGraph


class Complex(TopologyGraph):
    """
    Represents a host-guest complex topology graph.

    Examples
    --------
    *Construction*

    You can use :class:`.ConstructedMolecule` instances as the host,
    but you should turn them into a :class:`.BuildingBlock` first

    .. code-block:: python

        import stk

        host = stk.ConstructedMolecule(
            topology_graph=stk.cage.FourPlusSix(
                building_blocks=(
                    stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
                    stk.BuildingBlock(
                        smiles='BrCC(Br)CBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
            ),
        )
        complex1 = stk.ConstructedMolecule(
            topology_graph=stk.host_guest.Complex(
                host=stk.BuildingBlock.init_from_molecule(host),
                guest=stk.BuildingBlock('[Br][Br]'),
            ),
        )

    *Changing the Position of the Guest*

    You can change the position and orientation of the guest, as well
    as its displacement

    .. code-block:: python

        complex2 = stk.ConstructedMolecule(
            topology_graph=stk.host_guest.Complex(
                host=stk.BuildingBlock.init_from_molecule(host),
                guest=stk.BuildingBlock('[Br][Br]'),
                # Apply a rotation onto the guest molecule such that
                # the vector returned by get_direction() has the same
                # direction as [1, 1, 1].
                guest_start=guest.get_direction(),
                guest_target=[1, 1, 1],
                # Change the displacement of the guest.
                displacement=[5.3, 2.1, 7.1],
            ),
        )

    """

    def __init__(
        self,
        host,
        guest,
        guest_start=None,
        guest_target=None,
        displacement=(0, 0, 0),
        num_processes=1,
    ):
        """
        Initialize an instance of :class:`.Complex`.

        Parameters
        ----------
        host : :class:`.BuildingBlock`
            The host molecule.

        guest : :class:`.BuildingBlock`
            The guest molecule.

        guest_start : :class:`tuple` of :class:`float`, optional
            A direction vector which gets aligned with `guest_target`.

        guest_target : :class:`tuple` of :class:`float`, optional
            A direction vector which determines the rotation applied to
            the `guest`. A rotation such that `guest_start` is
            transformed into `guest_target` is applied.

        displacement : :class:`tuple` of :class:`float`, optional
            The translational offset of the guest.

        num_processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        Raises
        ------
        :class:`TypeError`
            If `guest_start` or `guest_target` is defined but the other
            is not.

        """

        num_nones = sum(
            1 for vector in (guest_start, guest_target)
            if vector is None
        )
        if num_nones == 1:
            raise TypeError(
                'If guest_start or guest_target is defined, '
                'the other must be too.'
            )

        if guest_start is None:
            start = target = (1., 0., 0.)
        else:
            start = guest_start = tuple(guest_start)
            target = guest_target = tuple(guest_target)

        # Save the values as None, for __repr__.
        self._guest_start = guest_start
        self._guest_target = guest_target
        self._displacement = displacement

        super().__init__(
            building_block_vertices={
                host: (_HostVertex(0, [0, 0, 0]), ),
                guest: (
                    _GuestVertex(1, displacement, start, target),
                ),
            },
            edges=(),
            reaction_factory=None,
            construction_stages=(),
            num_processes=num_processes,
            edge_groups=(),
        )

    def clone(self):
        clone = super().clone()
        clone._guest_start = self._guest_start
        clone._guest_target = self._guest_target
        clone._displacement = self._displacement
        return clone

    def _run_reactions(self, state):
        return state

    def _get_scale(self, building_block_vertices):
        return 1

    def __repr__(self):
        return (
            f'host_guest.Complex('
            f'guest_start={self._guest_start!r}, '
            f'guest_target={self._guest_target!r}, '
            f'displacement={self._displacement!r})'
        )
