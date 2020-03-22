"""
Dibromo Key
===========

"""


class DibromoKey:
    """
    Generates keys for :class:`.DibromoKey` instances.

    """

    # Keep the empty __init__() method for the docstring.
    def __init__(self):
        """
        Initialize a :class:`.DibromoKey` instance.

        """

        return

    def get_key(self, dibromo):
        """
        Get the key of `dibromo`.

        Parameters
        ----------
        dibromo : :class:`.Dibromo`
            The dibromo for which a key is needed.

        Returns
        -------
        :class:`object`
            The key.

        """

        return (
            'dibromo',
            dibromo.get_bromine1().get_id(),
            dibromo.get_atom1().get_id(),
            dibromo.get_bromine2().get_id(),
            dibromo.get_atom2().get_id(),
            tuple(sorted(dibromo.get_bonder_ids())),
            tuple(sorted(dibromo.get_deleter_ids())),
        )
