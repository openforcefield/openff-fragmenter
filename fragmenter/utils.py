import os
import time
import logging
import sys
import numpy as np
import itertools
import re
import codecs
import warnings
import copy

# def write_oedatabase(moldb, ofs, mlist, size):
#     """
#     This function writes out a new oedatabase from an existing database
#
#     Parameters
#     ----------
#     moldb: OpenEye database object
#     ofs: output stream
#     mlist: list of indices of molecules in moldb for new database
#     size: size of new database
#
#     """
#     for molidx in mlist[:size:]:
#         moldb.WriteMolecule(ofs, molidx)


# def create_oedatabase_idxfile(ifname):
#     """
#     This function creates an index file associated with a given molecule filename. It write out the file in the same
#     directory the parent molecule file is and adds an .idx extension to parent molecule filename.
#
#     From OpenEye's documentation:
#     The index file of a molecule dtabase stores file position offsets of the molecules in the file. Generating an index
#     file can be expensive, but it can be created only once and then it can speed up the handling of large molecule files
#     significantly
#
#     Parameters
#     ----------
#     ifname: str
#         absolute path to molecule file
#     """
#     idx_fname = oechem.OEGetMolDatabaseIdxFileName(ifname)
#
#     if os.path.exists(idx_fname):
#         oechem.OEThrow.Warning("{} index file already exists".format(idx_fname))
#     elif not oechem.OECreateMolDatabaseIdx(ifname):
#         oechem.OEThrow.Warning("Unable to create {} molecule index file".format(idx_fname))
#


# def to_psi4_input(fragment, molecule_name, memory='500 mb', command='gradient', symmetry='C1', crank_job='crank_job_1', filename=None):
#     """
#     Write out psi4 input for crank-launch input
#     Parameters
#     ----------
#     fragment: dict
#         JSON crank jobs for fragment
#     molecule_name: str
#         molecule name
#     memory: str
#         memory for psi4 job
#     command: str
#         gradient or optimize. Defualt is gradient. Must be gradient if using geomeTRIC for geometry optimization
#     symmetry: str
#         point group symmetry. Default is C1
#     crank_job: str
#         key for crank job in fragment dictionary. Default is crank_job_1
#     filename: str
#         Name for psi4 input file. Default is None. If None, function will return input string.
#
#     Returns
#     -------
#     psi4_input: str
#
#     """
#     job = fragment['crank_torsion_drives'][crank_job]
#     psi4_input = ""
#
#     if memory is not None:
#         psi4_input += "memory {}\n".format(memory)
#
#     psi4_input += "\nmolecule {}".format(molecule_name)
#     psi4_input += " {\n"
#     psi4_input += "symmetry {}\n".format(symmetry)
#     charge = fragment['molecule']['molecular_charge']
#     multiplicity = fragment['molecule']['molecular_multiplicity']
#     psi4_input += "{}  {} \n".format(charge, multiplicity)
#
#     # Convert 1-D Geom to 2-D geometry
#     n_atoms = len(fragment['molecule']['symbols'])
#     coords_3d = np.array(fragment['molecule']['geometry'], dtype=float).reshape(n_atoms, 3)
#
#     for element, coord in zip(fragment['molecule']['symbols'], coords_3d):
#         psi4_input += "%-7s %13.7f %13.7f %13.7f\n" % (element, coord[0], coord[1], coord[2])
#
#     psi4_input += "units Angstrom\n"
#     psi4_input += "}\n"
#
#     psi4_input += "\nset {\n"
#     psi4_input += " basis {}\n ".format(job['crank_specs']['model']['basis'])
#     try:
#         options = job['crank_specs']['options']
#         for option in options:
#             psi4_input += "{} {}\n".format(option, options[option])
#     except KeyError:
#         # No options were given
#         logger().info('no options found')
#         pass
#
#     psi4_input += "}\n"
#
#     psi4_input += "\n{}('{}')\n".format(command, job['crank_specs']['model']['method'])
#
#     if filename:
#         f = open(filename, 'w')
#         f.write(psi4_input)
#         f.close()
#     else:
#         return psi4_input
#
#
# def to_dihedraltxt(fragment, crank_job='crank_job_1', filename=None):
#     """
#     Generate dihedral txt file defining dihedrals that should be restraint in crank job
#
#     Parameters
#     ----------
#     fragment: dict
#         JSON crank jobs for fragment
#     crank_job: str
#         key for crank job to run
#     filename: str
#         filename of dihedral file. If None, returns input string.
#
#     Returns
#     -------
#     dih_input: str
#
#     """
#     needed_torsions = fragment['needed_torsion_drives']
#     crank_jobs = fragment['crank_torsion_drives'][crank_job]
#     dihedrals = []
#     for torsion in needed_torsions:
#         if torsion in crank_jobs:
#             dihedrals.append([i-1 for i in needed_torsions[torsion]])
#
#     dih_input = ""
#     dih_input += "# dihedral definition by atom indices starting from 0\n"
#     dih_input += "# i     j     k     l\n"
#     for dihedral in dihedrals:
#         dih_input += "  {}     {}     {}     {}\n".format(dihedral[0], dihedral[1], dihedral[2], dihedral[3])
#     if filename:
#         f = open(filename, 'w')
#         f.write(dih_input)
#         f.close()
#     else:
#         return dih_input
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