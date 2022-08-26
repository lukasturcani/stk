def assert_database_state(state1, state2):
    assert state1 == state2


class HashableDict(dict):
    def __hash__(self):
        return hash((frozenset(self), frozenset(self.values())))

    def __eq__(self, other):
        return super().__eq__(other)


def to_hashable(items):
    for key, value in items:
        if key == "a":
            yield key, tuple(tuple(v) for v in value)

        elif key == "b":
            yield (
                key,
                tuple(
                    (atom1, atom2, order, tuple(periodicity))
                    for atom1, atom2, order, periodicity in value
                ),
            )

        elif key == "m":
            yield key, tuple(tuple(v) for v in value)

        elif key == "BB":
            yield key, tuple(map(HashableDict, value))

        elif key == "aI":
            yield key, tuple(tuple(v) for v in value)

        elif key == "bI":
            yield key, tuple(tuple(v) for v in value)

        elif key == "nBB":
            yield key, tuple(value)

        else:
            yield key, value


def get_entry(item):
    """
    Create a :class:`.DatabaseEntry` from a :class:`dict`.

    Parameters
    ----------
    item : :class:`dict`
        A entry retrieved from a MongoDB collection.

    Returns
    -------
    :class:`.DatabaseEntry`
        The entry created from `item`.

    """

    return DatabaseEntry(
        **{
            key: value
            for key, value in to_hashable(item.items())
            if key != "_id"
        }
    )


class DatabaseEntry:
    """
    Represents an entry in a MongoDB collection.

    """

    def __init__(self, **kwargs):
        """
        Initialize a :class:`.DatabaseEntry`.

        Parameters
        ----------
        kwargs : :class:`dict`
            The data held by the entry.

        """

        self._items = kwargs

    def __eq__(self, other):
        return self._items == other._items

    def __hash__(self):
        return hash(
            (
                frozenset(self._items),
                frozenset(self._items.values()),
            )
        )

    def __repr__(self):
        return f"DatabaseEntry({self._items})"

    def __str__(self):
        return repr(self)


class DatabaseState:
    """
    Represents the state of a MongoDB database.

    """

    def __init__(self, entries):
        """
        Initialize a :class:`.DatabaseState`.

        Parameters
        ----------
        entries : :class:`dict`
            Maps a :class:`.DatabaseEntry` instance to the number
            of times an entry with that data appears in the database.

        """

        self._entries = dict(entries)

    def __eq__(self, other):
        return self._entries == other._entries

    def __repr__(self):
        return f"DatabaseState({self._entries})"

    def __str__(self):
        return repr(self)
