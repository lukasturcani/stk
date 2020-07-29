from collections import Counter

from ...utilities import DatabaseState, get_entry


def get_database_state(database):
    """
    Get the state of a :class:`.ValueMongoDb`.

    Parameters
    ----------
    database : :class:`.ValueMongoDb`
        The database whose state is wanted.

    Returns
    -------
    :class:`.DatabaseState`
        The current state of `database`.

    """

    entries = Counter()
    entries.update(map(get_entry, database._values.find({})))
    return DatabaseState(entries)
