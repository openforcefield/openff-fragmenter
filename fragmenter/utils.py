import os
import time
import logging
import sys
import json
from uuid import UUID
import numpy as np

from openeye import oechem, oeiupac, oedepict
from openmoltools import openeye


def write_oedatabase(moldb, ofs, mlist, size):
    """
    This function writes out a new oedatabase from an existing database

    Parameters
    ----------
    moldb: OpenEye database object
    ofs: output stream
    mlist: list of indices of molecules in moldb for new database
    size: size of new database

    """
    for molidx in mlist[:size:]:
        moldb.WriteMolecule(ofs, molidx)


def to_smi(smiles, path, base, return_fname=False):
    """
    This function writes out an .smi file for a list of SMILES
    Parameters
    ----------
    smiles: list of SMILES.
        The list can also contain strings that include name for SMILES separated by a space. ("SMILES Name")
    path: str
        path to output file
    base: str
        base name for output file

    """
    fname = os.path.join(path, base + '.smi')
    outf = open(fname, 'w')
    smiles_list = map(lambda x: x+"\n", list(smiles))
    outf.writelines(smiles_list)
    outf.close()
    if return_fname:
        return fname


def create_oedatabase_idxfile(ifname):
    """
    This function creates an index file associated with a given molecule filename. It write out the file in the same
    directory the parent molecule file is and adds an .idx extension to parent molecule filename.

    From OpenEye's documentation:
    The index file of a molecule dtabase stores file position offsets of the molecules in the file. Generating an index
    file can be expensive, but it can be created only once and then it can speed up the handling of large molecule files
    significantly

    Parameters
    ----------
    ifname: str
        absolute path to molecule file
    """
    idx_fname = oechem.OEGetMolDatabaseIdxFileName(ifname)

    if os.path.exists(idx_fname):
        oechem.OEThrow.Warning("{} index file already exists".format(idx_fname))
    elif not oechem.OECreateMolDatabaseIdx(ifname):
        oechem.OEThrow.Warning("Unable to create {} molecule index file".format(idx_fname))


def new_output_stream(outname):
    """
    This function creates a new oechem.oemolostream.
    Parameters
    ----------
    outname: str
        name of outputfile.

    Returns
    -------
    ofs: oechem.oemolostream

    """
    ofs = oechem.oemolostream()
    if not ofs.open(outname):
        oechem.OEThrow.Fatal("Unable to open {} for writing".format(outname))
    return ofs


def to_oemol(filename, title=True, verbose=True):
    """Create OEMol from file. If more than one mol in file, return list of OEMols.

    Parameters
    ----------
    filename: str
        absolute path to
    title: str, title
        title for molecule. If None, IUPAC name will be given as title.

    Returns
    -------
    mollist: list
        list of OEMol for multiple molecules. OEMol if file only has one molecule.
    """

    if not os.path.exists(filename):
        raise Exception("File {} not found".format(filename))
    if verbose:
        logger().info("Loading molecules from {}".format(filename))

    ifs = oechem.oemolistream(filename)
    #moldb = oechem.OEMolDatabase(ifs)
    mollist = []

    molecule = oechem.OECreateOEGraphMol()
    while oechem.OEReadMolecule(ifs, molecule):
        molecule_copy = oechem.OEMol(molecule)
        if title:
            title = molecule_copy.GetTitle()
            if verbose:
                logger().info("Reading molecule {}".format(title))

        mollist.append(normalize_molecule(molecule_copy, title))

    if len(mollist) <= 1:
        mollist = mollist[0]

    ifs.close()

    return mollist


def normalize_molecule(molecule, title=''):
    """Normalize a copy of the molecule by checking aromaticity, adding explicit hydrogens and renaming by IUPAC name
    or given title

    Parameters
    ----------
    molecule: OEMol
        The molecule to be normalized:
    title: str
        Name of molecule. If the string is empty, will use IUPAC name

    Returns
    -------
    molcopy: OEMol
        A (copied) version of the normalized molecule

    """
    molcopy = oechem.OEMol(molecule)

    # Assign aromaticity.
    oechem.OEAssignAromaticFlags(molcopy, oechem.OEAroModelOpenEye)

    # Add hydrogens.
    oechem.OEAddExplicitHydrogens(molcopy)

    # Set title to IUPAC name.
    name = title
    if not name:
        name = oeiupac.OECreateIUPACName(molcopy)
    molcopy.SetTitle(name)

    # Check for any missing atom names, if found reassign all of them.
    if any([atom.GetName() == '' for atom in molcopy.GetAtoms()]):
        oechem.OETriposAtomNames(molcopy)
    return molcopy


def png_atoms_labeled(smiles, fname):
    """Write out png file of molecule with atoms labeled with their map index.

    Parameters
    ----------
    smiles: str
        SMILES
    fname: str
        absolute path and filename for png

    """

    mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol, smiles)
    oedepict.OEPrepareDepiction(mol)

    width, height = 300, 200

    opts = oedepict.OE2DMolDisplayOptions(width, height, oedepict.OEScale_AutoScale)
    opts.SetAtomPropertyFunctor(oedepict.OEDisplayAtomMapIdx())
    opts.SetAtomPropLabelFont(oedepict.OEFont(oechem.OEDarkGreen))

    disp = oedepict.OE2DMolDisplay(mol, opts)
    return oedepict.OERenderMolecule(fname, disp)


def png_wiberg_labels(mol, fname, width=600, height=400):
    """
    Generate png figure of molecule. Bonds are labeled with Wiberg bond order

    Parameters
    ----------
    mol: OpenEye OEMol
    fname: str
        filename for png
    width: int
    height: int

    Returns
    -------
    bool:
    """

    oedepict.OEPrepareDepiction(mol)


    opts = oedepict.OE2DMolDisplayOptions(width, height, oedepict.OEScale_AutoScale)
    # opts.SetAtomPropertyFunctor(oedepict.OEDisplayAtomIdx())
    # opts.SetAtomPropLabelFont(oedepict.OEFont(oechem.OEDarkGreen))

    bondlabel = LabelBondOrder()
    opts.SetBondPropertyFunctor(bondlabel)

    disp = oedepict.OE2DMolDisplay(mol, opts)
    return oedepict.OERenderMolecule(fname, disp)


class LabelBondOrder(oedepict.OEDisplayBondPropBase):
    def __init__(self):
        oedepict.OEDisplayBondPropBase.__init__(self)

    def __call__(self, bond):
        bondOrder = bond.GetData('WibergBondOrder')
        label = "{:.2f}".format(bondOrder)
        return label

    def CreateCopy(self):
        copy = LabelBondOrder()
        return copy.__disown__()


def mol2_to_psi4json(infile):
    """

    Parameters
    ----------
    infile

    Returns
    -------

    """
    pass


def mol_to_tagged_smiles(infile, outfile):
    """
    Generate .smi from input mol with index-tagged explicit hydrogen SMILES
    Parameters
    ----------
    infile: str
        input molecule file
    outfile: str
        output smi file. Must be smi or ism

    """
    ifs = oechem.oemolistream()
    if not ifs.open(infile):
        oechem.OEThrow.Fatal("Unable to open {} for reading".format(infile))

    ofs = oechem.oemolostream()
    if not ofs.open(outfile):
        oechem.OEThrow.Fatal("Unable to open {} for writing".format(outfile))
    if ofs.GetFormat() not in [oechem.OEFormat_ISM, oechem.OEFormat_SMI]:
        oechem.OEThrow.Fatal("Output format must be SMILES")

    for mol in ifs.GetOEMols():
        smiles = create_mapped_smiles(mol)
        oechem.OEWriteMolecule(ofs, mol)


def create_mapped_smiles(molecule):
    """
    Generate an index-tagged explicit hydrogen SMILES.
    Exmaple:
    SMILES string for carbon monoxide "CO"
    With index-tagged explicit hydrogen SMILES this becomes
    '[H:3][C:1]([H:4])([H:5])[O:2][H:6]'

    Parameters
    ----------
    molecule: OEMOl

    Returns
    -------
    index-tagged explicit hydrogen SMILES str

    """
    #ToDo check if tags already exist raise warning about overwritting existing tags. Maybe also add an option to override existing tags
    oechem.OEAddExplicitHydrogens(molecule)

    for atom in molecule.GetAtoms():
        atom.SetMapIdx(atom.GetIdx() + 1)

    # add tag to data
    tag = oechem.OEGetTag("has_map")
    molecule.SetData(tag, bool(True))

    return oechem.OEMolToSmiles(molecule)


def get_atom_map(tagged_smiles, molecule=None, is_mapped=False):
    """
    Returns a dictionary that maps tag on SMILES to atom index in molecule.
    Parameters
    ----------
    tagged_smiles: str
        index-tagged explicit hydrogen SMILES string
    molecule: OEMol
        molecule to generate map for. If None, a new OEMol will be generated from the tagged SMILES, the map will map to
        this molecule and it will be returned.
    is_mapped: bool
        Default: False
        When an OEMol is generated from SMART string with tags, the tag will be stored in atom.GetMapIdx(). The index-tagged
        explicit-hydrogen SMILES are tagges SMARTS. Therefore, if a molecule was generated with the tagged SMILES, there is
        no reason to do a substructure search to get the right order of the atoms. If is_mapped is True, the atom map will be
        generated from atom.GetMapIdx().


    Returns
    -------
    molecule: OEMol
        The molecule for the atom_map. If a molecule was provided, it's that molecule.
    atom_map: dict
        a dictionary that maps tag to atom index {tag:idx}
    """
    if molecule is None:
        molecule = openeye.smiles_to_oemol(tagged_smiles)
        # Since the molecule was generated from the tagged smiles, the mapping is already in the molecule.
        is_mapped = True

    # Check if conformer was generated. The atom indices can get reordered when generating conformers and then the atom
    # map won't be correct
    if molecule.GetMaxConfIdx() <= 1:
        for conf in molecule.GetConfs():
            values = np.asarray([conf.GetCoords().__getitem__(i) == (0.0, 0.0, 0.0) for i in
                                range(conf.GetCoords().__len__())])
        if values.all():
            # Generate on Omega conformer
            molecule = openeye.generate_conformers(molecule, max_confs=1)

    if is_mapped:
        atom_map = {}
        for atom in molecule.GetAtoms():
            atom_map[atom.GetMapIdx()] = atom.GetIdx()
        return molecule, atom_map

    ss = oechem.OESubSearch(tagged_smiles)
    oechem.OEPrepareSearch(molecule, ss)
    ss.SetMaxMatches(1)

    atom_map = {}
    t1 = time.time()
    matches = [m for m in ss.Match(molecule)]
    t2 = time.time()
    seconds = t2-t1
    logger().info("Substructure search took {} seconds".format(seconds))
    if not matches:
        logger().info("MCSS failed for {}, smiles: {}".format(molecule.GetTitle(), tagged_smiles))
        return False
    for match in matches:
        for ma in match.GetAtoms():
            atom_map[ma.pattern.GetMapIdx()] = ma.target.GetIdx()

    # sanity check
    mol = oechem.OEGraphMol()
    oechem.OESubsetMol(mol, match, True)
    logger().info("Match SMILES: {}".format(oechem.OEMolToSmiles(mol)))

    return molecule, atom_map


def to_mapped_xyz(molecule, atom_map, conformer=None, xyz_format=False, filename=None):
    """
    Generate xyz coordinates for molecule in the order given by the atom_map. atom_map is a dictionary that maps the
    tag on the SMILES to the atom idex in OEMol.
    Parameters
    ----------
    molecule: OEMol with conformers
    atom_map: dict
        maps tag in SMILES to atom index
    conformer: int
        Which conformer to write xyz file for. If None, write out all conformers. Default is None
    xyz_format: bool
        If True, will write out number of atoms and molecule name. If false, will only write out elements and coordinates
    filename: str
        Name of file to save to. If None, only returns a string.

    Returns
    -------
    str: elements and xyz coordinates in order of tagged SMILES

    """
    xyz = ""
    for k, mol in enumerate(molecule.GetConfs()):
        if k == conformer or conformer is None:
            if xyz_format:
                xyz += "{}\n".format(mol.GetMaxAtomIdx())
                xyz += "{}\n".format(mol.GetTitle())
            coords = oechem.OEFloatArray(mol.GetMaxAtomIdx() * 3)
            mol.GetCoords(coords)
            if k != 0 and not xyz_format:
                    xyz += "*"
            for mapping in range(1, len(atom_map)+1):
                idx = atom_map[mapping]
                atom = mol.GetAtom(oechem.OEHasAtomIdx(idx))
                syb = oechem.OEGetAtomicSymbol(atom.GetAtomicNum())
                xyz += "  {}      {:05.3f}   {:05.3f}   {:05.3f}\n".format(syb,
                                                                           coords[idx * 3],
                                                                           coords[idx * 3 + 1],
                                                                           coords[idx * 3 + 2])

    if filename:
        file = open("{}.xyz".format(filename, 'w'))
        file.write(xyz)
        file.close()
    return xyz


def to_mapped_QC_JSON_geometry(molecule, atom_map):
    """
    Generate xyz coordinates for molecule in the order given by the atom_map. atom_map is a dictionary that maps the
    tag on the SMILES to the atom idex in OEMol.
    Parameters
    ----------
    molecule: OEMol with conformers
    atom_map: dict
        maps tag in SMILES to atom index

    Returns
    -------
    dict: QC_JSON Molecule spec {symbols: [], geometry: []}

    """
    symbols = []
    geometry = []
    if molecule.GetMaxConfIdx() != 1:
        raise Warning("The molecule must have at least and at most 1 conformation")

    for mapping in range(1, len(atom_map)+1):
        idx = atom_map[mapping]
        atom = molecule.GetAtom(oechem.OEHasAtomIdx(idx))
        syb = oechem.OEGetAtomicSymbol(atom.GetAtomicNum())
        symbols.append(syb)
        for i in range(3):
            geometry.append(molecule.GetCoords()[idx][i])

    return {'symbols': symbols, 'geometry': geometry}


# def run_psi4_json(tagged_smiles, molecule, driver, method, basis, properties=None, return_output=True, xyz_file=True):
#     """
#
#     Parameters
#     ----------
#     tagged_smiles
#     xyz
#     driver
#     method
#     basis
#     properties
#     return_output
#     xyz_file
#
#     Returns
#     -------
#
#     """
#     json_data = {}
#     json_data["tagged_smiles"] = tagged_smiles
#     json_data["molecule"] = molecule
#     json_data["driver"] = "property"
#     json_data["kwargs"] = {"properties": properties}
#     json_data["method"] = method
#     json_data["options"] = {"BASIS": basis}
#     json_data["return_output"] = return_output
#
#     name = molecule.split("\n")[2]
#     if xyz_file:
#         file = open("{}.xyz".format(name), 'w')
#         file.write(molecule)
#         file.close()
#
#     j = json.dumb(json_data, indent=4, sort_keys=True)
#     f = open("{}.input.json".format(name), 'w')
#     f.close()
#
#     psi4.json_wrapper.run_json(json_data)
#     j = json.dump(json_data, indent=4, sort_keys=True)
#
#     f = open("{}.output.json".format(name))


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

class UUIDEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, UUID):
            # if the obj is uuid, we simply return the value of uuid
            return obj.hex
        return json.JSONEncoder.default(self, obj)