"""
Cage Key Maker
==============

"""

from ...building_block import BuildingBlockKeyMaker


class CageKeyMaker:
    """
    Makes keys for :class:`.Cage` instances.

    """

    def __init__(
        self,
        building_block_key_maker=BuildingBlockKeyMaker(),
    ):
        """
        Initialize a :class:`.CageKeyMaker` instance.

        Parameters
        ----------
        building_block_key_maker : \
                :class:`.BuildingBlockKeyMaker`, optional
            Used to create the keys of the cage's building blocks.

        """

        return

    def get_key(self, cage):
        """
        Get the key of `cage`.

        Parameters
        ----------
        cage : :class:`.Cage`
            The cage for which a key is wanted.

        Returns
        -------
        :class:`object`
            The key.

        """

        building_blocks = (
            building_block.with_canonical_atom_ordering()
            for building_block in cage.get_building_blocks()
        )
        # Cage does not provide a public interface to for these
        # attributes, so private ones must be used.

        key = (
            cage.__class__.__name__,
            vertex_alignments,
            tuple(sorted(map(
                self._building_block_key_maker.get_key,
                building_blocks,
            ))),
        )
        return str(key)
