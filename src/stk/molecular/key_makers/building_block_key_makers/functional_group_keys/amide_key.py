"""
Amide Key
=========

"""


class AmideKey:
    """
    Generates keys for :class:`.Amide` instances.

    """

    # Keep the empty __init__() method for the docstring.
    def __init__(self):
        """
        Initialize an :class:`.AmideKey` instance.

        """

        return

    def get_key(self, amide):
        """
        Get the key of `amide`.

        Parameters
        ----------
        amide : :class:`.Amide`
            The amide for which a key is needed.

        Returns
        -------
        :class:`object`
            The key.

        """

        return (
            'amide',
            amide.get_carbon().get_id(),
            amide.get_oxygen().get_id(),
            amide.get_nitrogen().get_id(),
            amide.get_hydrogen1().get_id(),
            amide.get_hydrogen2().get_id(),
            amide.get_atom().get_id(),
            tuple(sorted(amide.get_bonder_ids())),
            tuple(sorted(amide.get_deleter_ids())),
        )
