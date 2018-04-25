"""
Defines tools to make ``logging`` comptaible with ``multiprocessing``.

Based largely on: https://gist.github.com/schlamar/7003737

"""

import logging
import os

# Define the formatter for logging messages.
try:
    f = '\n' + '='*os.get_terminal_size().columns + '\n\n'
except OSError as ex:
    # When testing os.get_terminal_size() will fail because stdout is
    # not connceted to a terminal.
    f = '\n' + '='*100 + '\n\n'
formatter = logging.Formatter(fmt=f+('%(asctime)s - %(levelname)s - '
                                     '%(name)s - %(message)s'),
                              datefmt='%H:%M:%S')


# Define logging handlers.
errorhandler = logging.FileHandler('output/scratch/errors.log',
                                   delay=True)
errorhandler.setLevel(logging.ERROR)

streamhandler = logging.StreamHandler()

errorhandler.setFormatter(formatter)
streamhandler.setFormatter(formatter)


class MPLogger(logging.Logger):
    """

    """

    log_queue = None

    def isEnabledFor(self, level):
        return True

    def handle(self, record):
        """

        """

        ei = record.exc_info
        if ei:
            # To get traceback text into record.exc_text
            logging._defaultFormatter.format(record)
            record.exc_info = None  # Not needed any more.
        d = dict(record.__dict__)
        d['msg'] = record.getMessage()
        d['args'] = None
        self.log_queue.put(d)


def logged_call(log_queue, func, *args, **kwargs):
    """

    """

    MPLogger.log_queue = log_queue
    logging.setLoggerClass(MPLogger)
    # Monkey patch root logger and already defined loggers.
    logging.root.__class__ = MPLogger
    for logger in logging.Logger.manager.loggerDict.values():
        if not isinstance(logger, logging.PlaceHolder):
            logger.__class__ = MPLogger
    return func(*args, **kwargs)


def daemon_logger(log_queue):
    """
    Logs messages from subprocesses in the main process.

    Parameters
    ----------
    log_queue : :class:`multiprocessing.Queue`
        A queue into which log records are sent by subprocesses.

    Returns
    -------
    None : :class:`NoneType`

    """

    while True:
        try:
            record_data = log_queue.get()
            if record_data is None:
                break
            record = logging.makeLogRecord(record_data)

            logger = logging.getLogger(record.name)
            if logger.isEnabledFor(record.levelno):
                logger.handle(record)
        except (KeyboardInterrupt, SystemExit):
            raise
        except EOFError:
            break
        except Exception:
            logging.exception('Error in log handler.')
