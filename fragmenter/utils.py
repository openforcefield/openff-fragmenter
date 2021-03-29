import logging
import sys

BOHR_2_ANGSTROM = 0.529177210
ANGSROM_2_BOHR = 1. / BOHR_2_ANGSTROM
HARTREE_2_KJMOL = 2625.50
KJMOL_2_HARTREE = 1. / HARTREE_2_KJMOL


def logger(name='fragmenter', pattern='%(asctime)s %(levelname)s %(name)s: %(message)s',
           date_format='%H:%M:%S', handler=logging.StreamHandler(sys.stdout)):
    """
    Retrieves the logger instance associated to the given name

    :param name: The name of the logger instance
    :param pattern: The associated pattern
    :param date_format: The date format to be used in the pattern
    :param handler: The logging handler
    :return: The logger
    """
    _logger = logging.getLogger(name)
    _logger.setLevel(log_level(False))

    if not _logger.handlers:
        formatter = logging.Formatter(pattern, date_format)
        handler.setFormatter(formatter)
        handler.setLevel(log_level(False))
        _logger.addHandler(handler)
        _logger.propagate = False
    return _logger


def log_level(verbose=False):
    if verbose:
        return logging.DEBUG
    else:
        return logging.INFO
