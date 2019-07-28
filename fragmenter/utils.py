import logging
import sys
import json

"""
~~~~~~~~~~~~~~~~~~~
Conversion factors
~~~~~~~~~~~~~~~~~~~
"""
BOHR_2_ANGSTROM = 0.529177210
ANGSROM_2_BOHR = 1. / BOHR_2_ANGSTROM
HARTREE_2_KJMOL = 2625.50
KJMOL_2_HARTREE = 1. / HARTREE_2_KJMOL

"""
~~~~~~~
Logger
~~~~~~~
"""


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
    _logger.setLevel(log_level(verbose))

    if not _logger.handlers:
        formatter = logging.Formatter(pattern, date_format)
        handler.setFormatter(formatter)
        handler.setLevel(log_level(verbose))
        _logger.addHandler(handler)
        _logger.propagate = False
    return _logger

verbose = False


def log_level(verbose=verbose):
    if verbose:
        return logging.DEBUG
    else:
        return logging.INFO


"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Movie making functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""


def merge_images(file1, file2, filename):
    """Merge two images into one, displayed side by side
    :param file1: path to first image file
    :param file2: path to second image file
    :return: the merged Image object
    """
    try:
        from PIL import Image
    except ImportError:
        raise ImportError("Must have PIL installed to use this function")

    image1 = Image.open(file1)
    image2 = Image.open(file2)

    (width1, height1) = image1.size
    (width2, height2) = image2.size

    result_width = width1 + width2 - 500 # play with the number here
    result_height = max(height1, height2)

    result = Image.new('RGB', (result_width, result_height), color='white')
    result.paste(im=image2, box=(width1-250, 0)) # play with the number here
    result.paste(im=image1, box=(0, int(result_height/2) - int(height1/2)))
    result.save(filename)


def make_movie(files, fps=15, filename='movie'):
    """
    Make movie from list of files

    Parameters
    ----------
    files: list of paths to files
    fps: files per second, int, optional
        defualt 15
    filename: str, optional
        default 'movie

    """
    try:
        import moviepy.editor as mpy
    except ImportError:
        raise ImportError("You need to install moviepy to use this function")

    im = mpy.ImageSequenceClip(files, fps=fps)
    im.write_videofile('{}.mp4'.format(filename), fps=fps)

def serialize_bond(key):
    """
    Serialize bond tuple for JSON
    Parameters
    ----------
    key : tuple of ints

    Returns
    -------
    serialized tuple

    """
    if isinstance(key, (int, float)):
        key = (int(key), )

    return json.dumps(key)

def deserialize_bond(key):
    """
    Convert string tuple to tuple
    Parameters
    ----------
    key :

    Returns
    -------

    """
    return tuple(int(x) for x in key.split('[')[-1].split(']')[0].split(','))



