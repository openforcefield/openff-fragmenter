import logging
import sys
import re
import codecs
import copy
import numpy as np
from PIL import Image


"""
~~~~~~~~~~~~~~~~~~~
Conversion factors
~~~~~~~~~~~~~~~~~~~
"""
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


def sort_nicely(l):
    """

    Parameters
    ----------
    l

    Returns
    -------

    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    l.sort(key=alphanum_key)


def flatten(l, ltypes=(list, tuple)):
    """
    Flatten list of lists
    Parameters
    ----------
    l: list to flatten
    ltypes: tuple of types

    Returns
    -------
    flattened list
    """
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -=1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)


def make_python_identifier(string, namespace=None, reserved_words=None,
                           convert='hex', handle='force'):
    """
    Taken from https://gist.github.com/JamesPHoughton
    Takes an arbitrary string and creates a valid Python identifier.
    If the python identifier created is already in the namespace,
    or if the identifier is a reserved word in the reserved_words
    list, or is a python default reserved word,
    adds _1, or if _1 is in the namespace, _2, etc.
    Parameters
    ----------
    string : <basestring>
        The text to be converted into a valid python identifier
    namespace : <dictionary>
        Map of existing translations into python safe identifiers.
        This is to ensure that two strings are not translated into
        the same python identifier
    reserved_words : <list of strings>
        List of words that are reserved (because they have other meanings
        in this particular program, such as also being the names of
        libraries, etc.
    convert : <string>
        Tells the function what to do with characters that are not
        valid in python identifiers
        - 'hex' implies that they will be converted to their hexidecimal
                representation. This is handy if you have variables that
                have a lot of reserved characters
        - 'drop' implies that they will just be dropped altogether
    handle : <string>
        Tells the function how to deal with namespace conflicts
        - 'force' will create a representation which is not in conflict
                  by appending _n to the resulting variable where n is
                  the lowest number necessary to avoid a conflict
        - 'throw' will raise an exception
    Returns
    -------
    identifier : <string>
        A vaild python identifier based on the input string
    namespace : <dictionary>
        An updated map of the translations of words to python identifiers,
        including the passed in 'string'.
    Examples
    --------
    >>> make_python_identifier('Capital')
    ('capital', {'Capital': 'capital'})
    >>> make_python_identifier('multiple words')
    ('multiple_words', {'multiple words': 'multiple_words'})
    >>> make_python_identifier('multiple     spaces')
    ('multiple_spaces', {'multiple     spaces': 'multiple_spaces'})
    When the name is a python keyword, add '_1' to differentiate it
    >>> make_python_identifier('for')
    ('for_1', {'for': 'for_1'})
    Remove leading and trailing whitespace
    >>> make_python_identifier('  whitespace  ')
    ('whitespace', {'  whitespace  ': 'whitespace'})
    Remove most special characters outright:
    >>> make_python_identifier('H@t tr!ck')
    ('ht_trck', {'H@t tr!ck': 'ht_trck'})
    Replace special characters with their hex representations
    >>> make_python_identifier('H@t tr!ck', convert='hex')
    ('h40t_tr21ck', {'H@t tr!ck': 'h40t_tr21ck'})
    remove leading digits
    >>> make_python_identifier('123abc')
    ('abc', {'123abc': 'abc'})
    namespace conflicts
    >>> make_python_identifier('Variable$', namespace={'Variable@':'variable'})
    ('variable_1', {'Variable@': 'variable', 'Variable$': 'variable_1'})
    >>> make_python_identifier('Variable$', namespace={'Variable@':'variable', 'Variable%':'variable_1'})
    ('variable_2', {'Variable@': 'variable', 'Variable%': 'variable_1', 'Variable$': 'variable_2'})
    throw exception instead
    >>> make_python_identifier('Variable$', namespace={'Variable@':'variable'}, handle='throw')
    Traceback (most recent call last):
     ...
    NameError: variable already exists in namespace or is a reserved word
    References
    ----------
    Identifiers must follow the convention outlined here:
        https://docs.python.org/2/reference/lexical_analysis.html#identifiers
    """
    if namespace is None:
        namespace = {}

    if reserved_words is None:
        reserved_words = []

    s = copy.deepcopy(string)

    # remove leading and trailing whitespace
    s = s.strip()

    # Make spaces into underscores
    s = re.sub('[\\s\\t\\n]+', '_', s)

    if convert == 'hex':
        # Convert invalid characters to hex
        hexlify = codecs.getencoder('hex')
        s = ''.join(['_' + (hexlify(c.encode('utf-8'))[0]).decode('utf-8') + '_'
                     if re.findall('[^0-9a-zA-Z_]', c) else c for c in s])

    elif convert == 'drop':
        # Remove invalid characters
        s = re.sub('[^0-9a-zA-Z_]', '', s)

    # Remove leading characters until we find a letter or underscore
    s = re.sub('^[^a-zA-Z_]+', '', s)

    # Check that the string is not a python identifier
    while (#s in keyword.kwlist or
           s in namespace.values() or
           s in reserved_words):
        if handle == 'throw':
            raise NameError(s + ' already exists in namespace or is a reserved word')
        if handle == 'force':
            if re.match(".*?_\d+$", s):
                i = re.match(".*?_(\d+)$", s).groups()[0]
                s = s.strip('_'+i) + '_'+str(int(i)+1)
            else:
                s += '_1'

    namespace[string] = s

    return s, namespace

"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Util functions for manipulating QCFractal output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""


def grid_id_from_str(grid_id_str):
    """
    Only works for 1D grids
    Parameters
    ----------
    grid_id_str

    Returns
    -------

    """
    return int(grid_id_str.split('[')[-1].split(']')[0])


def sort_energies(final_energies):
    for frag in final_energies:
        for job in final_energies[frag]:
            angles = []
            energies = []
            for angle in final_energies[frag][job]:
                energy = final_energies[frag][job][angle]
                angle = grid_id_from_str(angle)
                angles.append(angle)
                energies.append(energy)
            energies = np.asarray(energies)
            energies = energies * HARTREE_2_KJMOL
            rel_energies = energies - energies.min()
            sorted_energies = [x for _, x in sorted(zip(angles, rel_energies))]
            sorted_angles = sorted(angles)
            final_energies[frag][job] = (sorted_angles, sorted_energies)


def deserialze_molecules(final_molecules):
    deserialized = {}
    for frag in final_molecules:
        deserialized[frag] = {}
        for job in final_molecules[frag]:
            deserialized[frag][job] = {}
            for angle in final_molecules[frag][job]:
                molecule = final_molecules[frag][job][angle]
                angle_int = grid_id_from_str(angle)
                deserialized[frag][job][angle_int] = molecule
    return deserialized

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
