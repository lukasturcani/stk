"""
COF Key Maker
=============

"""

from ...building_block import BuildingBlockKeyMaker


class CofKeyMaker:
    """
    Makes keys for :class:`.Cof` instances.

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
            Used to create the keys of the COF's building blocks.

        """

        return

    def get_key(self, cof):
        """
        Get the key of `cof`.

        Parameters
        ----------
        cof : :class:`.Cof`
            The cof for which a key is wanted.

        Returns
        -------
        :class:`object`
            The key.

        """

        building_blocks = (
            building_block.with_canonical_atom_ordering()
            for building_block in cof.get_building_blocks()
        )
        key = (
            cof.__class__.__name__,
            tuple(sorted(map(
                self._building_block_key_maker.get_key,
                building_blocks,
            ))),
        )
        return str(key)
