"""
Ring Amine Key
==============

"""


class RingAmineKey:
    """
    Generates keys for :class:`.RingAmineKey` instances.

    """

    # Keep the empty __init__() method for the docstring.
    def __init__(self):
        """
        Initialize a :class:`.RingAmineKey` instance.

        """

        return

    def get_key(self, ring_amine):
        """
        Get the key of `ring_amine`.

        Parameters
        ----------
        ringAmine : :class:`.Ring_amine`
            The ring amine functional group for which a key is needed.

        Returns
        -------
        :class:`object`
            The key.

        """

        return (
            'ring_amine',
            ring_amine.get_nitrogen().get_id(),
            ring_amine.get_hydrogen1().get_id(),
            ring_amine.get_hydrogen2().get_id(),
            ring_amine.get_carbon1().get_id(),
            ring_amine.get_carbon2().get_id(),
            ring_amine.get_hydrogen3().get_id(),
            ring_amine.get_carbon3().get_id(),
        )
