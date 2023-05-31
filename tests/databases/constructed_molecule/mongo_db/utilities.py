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
    entries.update(map(get_entry, database._molecules.find({})))
    entries.update(map(get_entry, database._position_matrices.find({})))
    entries.update(
        map(
            get_entry,
            database._constructed_molecules.find({}),
        )
    )
    entries.update(
        map(
            get_entry,
            database._building_block_position_matrices.find({}),
        )
    )
    return DatabaseState(entries)
