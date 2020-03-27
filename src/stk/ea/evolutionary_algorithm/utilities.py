import os
import logging

from stk.molecular import Inchi


def get_logger():
    # Define the formatter for logging messages.
    try:
        f = '\n' + '='*os.get_terminal_size().columns + '\n\n'
    except OSError:
        # os.get_terminal_size() may fail because stdout is
        # not connected to a terminal.
        f = '\n' + '='*100 + '\n\n'
    formatter = logging.Formatter(
        fmt=f+('%(asctime)s - %(levelname)s - %(name)s - %(message)s'),
        datefmt='%H:%M:%S'
    )

    # Define logging handlers.
    errorhandler = logging.FileHandler(
        'ea_errors.log',
        delay=True
    )
    errorhandler.setLevel(logging.ERROR)
    streamhandler = logging.StreamHandler()

    errorhandler.setFormatter(formatter)
    streamhandler.setFormatter(formatter)

    # Get the loggers.
    rootlogger = logging.getLogger()
    rootlogger.addHandler(errorhandler)
    rootlogger.addHandler(streamhandler)
    return logging.getLogger(__name__)


def get_inchi(record):
    """
    Get the InChI of a :class:`.MoleculeRecord`.

    Parameters
    ----------
    record : :class:`.MoleculeRecord`
        Holds a molecule.

    Returns
    -------
    :class:`str`
        The InChI of the molecule in `record`.

    """

    return Inchi().get_key(record.get_molecule())
