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

from fragmenter import chemi
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


def to_smi(smiles, filename, return_fname=False):
    """
    This function writes out an .smi file for a list of SMILES
    Parameters
    ----------
    smiles: list of SMILES.
        The list can also contain strings that include name for SMILES separated by a space. ("SMILES Name")
    filename: str
        name of output file
    return_fname: bool, optional, default=False
        If True, returns absolute path to filename.

    """

    smiles_list = map(lambda x: x+"\n", list(smiles))
    with open(filename, 'w') as outf:
        outf.writelines(smiles_list)

    if return_fname:
        filenmae = os.path.join(os.getcwd(), filename)
        return filename


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


def file_to_oemols(filename, title=True, verbose=True):
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

    # if len(mollist) <= 1:
    #     mollist = mollist[0]

    ifs.close()

    return mollist


def smiles_to_oemol(smiles, name='', normalize=True):
    """Create a OEMolBuilder from a smiles string.
    Parameters
    ----------
    smiles : str
        SMILES representation of desired molecule.
    Returns
    -------
    molecule : OEMol
        A normalized molecule with desired smiles string.
    """

    molecule = oechem.OEMol()
    if not oechem.OEParseSmiles(molecule, smiles):
        raise ValueError("The supplied SMILES '%s' could not be parsed." % smiles)

    if normalize:
        molecule = normalize_molecule(molecule, name)

    return molecule


def oemols_to_smiles_list(OEMols, isomeric=True):

    if not isinstance(OEMols, list):
        OEMols = [OEMols]

    SMILES = []
    for mol in OEMols:
        SMILES.append(create_mapped_smiles(mol, tagged=False, explicit_hydrogen=False, isomeric=isomeric))

    return SMILES


def file_to_smiles_list(filename, return_titles=True, isomeric=True):

    oemols = file_to_oemols(filename)
    # Check if oemols have names
    if return_titles:
        names = []
        for mol in oemols:
            title = mol.GetTitle()
            if not title:
                logger().warning("an oemol does not have a name. Adding an empty str to the titles list")
            names.append(title)
    smiles_list = oemols_to_smiles_list(oemols, isomeric=isomeric)

    if return_titles:
        return smiles_list, names
    return smiles_list


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


def check_molecule(molecule, title=''):

    if isinstance(molecule, oechem.OEMol):
        return molecule

    mol = oechem.OEMol()
    # First try reading as smiles
    if not oechem.OESmilesToMol(mol, molecule):
        # Try reading as input file
        ifs = oechem.oemolistream()
        if not ifs.open(molecule):
            raise Warning('Could not parse molecule.')

    # normalize molecule
    if not title:
        title = mol.GetTitle()
    molecule = normalize_molecule(mol, title=title)
    return molecule


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


def create_mapped_smiles(molecule, tagged=True, explicit_hydrogen=True, isomeric=True):
    """
    Generate an index-tagged explicit hydrogen SMILES.
    Exmaple:
    SMILES string for carbon monoxide "CO"
    With index-tagged explicit hydrogen SMILES this becomes
    '[H:3][C:1]([H:4])([H:5])[O:2][H:6]'

    **Warning**
    The map that is added to the SMILES is not guaranteed to correspond to the indices to the molecule provided since a
    new conformer is generated to ensure explicit stereochemistry in the SMILES.

    Parameters
    ----------
    molecule: OEMOl
    tagged: Bool, optional, default=True
        If True, will add index tags to each atom. If True, explicit_hydrogen should also be tru so that hydrogens are
        numbered as well.
    explicit_hydrogen: Bool, optional, default=True
        If True, will write explicit hydrogens in SMILES
    isomeric: Bool, optiona, default=True
        If True, will specify cis/trans and R/S isomerism

    Returns
    -------
    index-tagged explicit hydrogen SMILES str

    """

    # Generate conformer if not conformers exist
    if isomeric:
        if not has_conformer(molecule):
            # Check if conformer already exists
            try:
                # Set strict_stereo so that we get a conformer to perceive chirality for isomeric SMILES
                # Set strict type to False so conformer will be generated even if exact MMFF atom type does not exist for
                # those atoms. This is probably a bigger problem for charging than conformer generation ?
                # Set copy to False so that input molecule gets tagged and mapped consistent with SMILES
                molecule = chemi.generate_conformers(molecule, max_confs=1, strict_stereo=False, strict_types=False, copy=False)
            except RuntimeError:
                        logger().warning("Omega failed to generate a conformer for {}. Steroechemistry might not be specified"
                                      "in tagged SMILES".format(oechem.OEMolToSmiles(molecule)))

        oechem.OEPerceiveChiral(molecule)
        oechem.OE3DToAtomStereo(molecule)
        oechem.OE3DToBondStereo(molecule)

    #iso_smiles = oechem.OEMolToSmiles(molecule)
    # Create a new molecule with isomeric SMILES
    # molecule = oechem.OEMol()
    # oechem.OESmilesToMol(molecule, iso_smiles)
    # if molecule.GetMaxConfIdx() <= 1:
    #     for conf in molecule.GetConfs():
    #         values = np.asarray([conf.GetCoords().__getitem__(i) == (0.0, 0.0, 0.0) for i in
    #                             range(conf.GetCoords().__len__())])
    #     if values.all():
    #         # No conformers were generated so generate one
    #         molcopy = copy.deepcopy(molecule)
    #         molecule = chemi.generate_conformers(molcopy, max_confs=1, strict_stereo=False)
    #         iso_smiles = oechem.OEMolToSmiles(molecule)
    #
    #         molecule = oechem.OEMol()
    #         oechem.OESmilesToMol(molecule, iso_smiles)

    if not explicit_hydrogen and not tagged and isomeric:
        return oechem.OEMolToSmiles(molecule)
    if not explicit_hydrogen and not tagged and not isomeric:
         return oechem.OECreateSmiString(molecule, oechem.OESMILESFlag_Canonical | oechem.OESMILESFlag_RGroups)

    oechem.OEAddExplicitHydrogens(molecule)
    if not tagged and explicit_hydrogen and isomeric:
        return oechem.OECreateSmiString(molecule, oechem.OESMILESFlag_Hydrogens | oechem.OESMILESFlag_ISOMERIC)

    if not tagged and explicit_hydrogen and not isomeric:
        return oechem.OECreateSmiString(molecule, oechem.OESMILESFlag_Hydrogens | oechem.OESMILESFlag_Canonical |
                                        oechem.OESMILESFlag_RGroups)

    # Add tags to molecule
    for atom in molecule.GetAtoms():
        atom.SetMapIdx(atom.GetIdx() + 1)

    if tagged and not explicit_hydrogen:
        raise Warning("Tagged SMILES must include hydrogens to retain order")

    for atom in molecule.GetAtoms():
        atom.SetMapIdx(atom.GetIdx() + 1)

    if tagged and not isomeric:
        raise Warning("Tagged SMILES must include stereochemistry ")

    # add tag to data
    tag = oechem.OEGetTag("has_map")
    molecule.SetData(tag, bool(True))

    return oechem.OEMolToSmiles(molecule)


def has_conformer(molecule, check_two_dimension=False):
    """
    Check if conformer exists for molecule. Return True or False
    Parameters
    ----------
    molecule
    check_two_dimension: bool, optional. Default False
        If True, will also check if conformation is a 2D conformation (all z coordinates are zero) and return False if
        conformation is 2D

    Returns
    -------

    """
    conformer_bool = True
    try:
        if molecule.NumConfs() <= 1:
            # Check if xyz coordinates are not zero
            for conf in molecule.GetConfs():
                # print(conf.GetCoords().__len__())
                # coords = molecule.GetCoords()
                # values = np.asarray(list(coords.values()))
                # print(values)
                # print(values.all())
                # if not values.all():
                #     conformer_bool = False
                #for i in range(conf.GetCoords().__len__()):
                values = np.asarray([conf.GetCoords().__getitem__(i) == (0.0, 0.0, 0.0) for i in
                                    conf.GetCoords()])
            if values.all():
                conformer_bool = False
    except AttributeError:
        conformer_bool = False

    if conformer_bool and check_two_dimension:
        for conf in molecule.GetConfs():
            values = np.asarray([conf.GetCoords().__getitem__(i)[-1] == 0.0 for i in conf.GetCoords()])
            if values.all():
                conformer_bool = False
    return conformer_bool


def is_mapped(molecule):
    """
    Check if atoms are mapped. If any atom is missing a tag, this will return False
    Parameters
    ----------
    molecule: OEMol

    Returns
    -------
    Bool: True if mapped. False otherwise
    """
    IS_MAPPED = True
    for atom in molecule.GetAtoms():
        if atom.GetMapIdx() == 0:
            IS_MAPPED = False
    return IS_MAPPED


def get_atom_map(tagged_smiles, molecule=None, StrictStereo=True, verbose=True):
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

    # Check if conformer was generated. The atom indices can get reordered when generating conformers and then the atom
    # map won't be correct
    if not has_conformer(molecule):
    # if molecule.GetMaxConfIdx() <= 1:
    #     for conf in molecule.GetConfs():
    #         values = np.asarray([conf.GetCoords().__getitem__(i) == (0.0, 0.0, 0.0) for i in
    #                             range(conf.GetCoords().__len__())])
    #     if values.all():
        # Generate an Omega conformer
        try:
            molecule = openeye.generate_conformers(molecule, max_confs=1, strictStereo=StrictStereo, strictTypes=False)
            # Omega can change the ordering so whatever map existed is not there anymore
        except RuntimeError:
            logger().warning("Omega failed to generate a conformer for {}. Mapping can change when a new conformer is "
                          "generated".format(molecule.GetTitle()))

    # Check if molecule is mapped
    mapped = is_mapped(molecule)
    if mapped:
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
    logger().debug("Substructure search took {} seconds".format(seconds))
    if not matches:
        logger().info("MCSS failed for {}, smiles: {}".format(molecule.GetTitle(), tagged_smiles))
        return False
    for match in matches:
        for ma in match.GetAtoms():
            atom_map[ma.pattern.GetMapIdx()] = ma.target.GetIdx()

    # sanity check
    mol = oechem.OEGraphMol()
    oechem.OESubsetMol(mol, match, True)
    logger().debug("Match SMILES: {}".format(oechem.OEMolToSmiles(mol)))

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
        file = open("{}.xyz".format(filename), 'w')
        file.write(xyz)
        file.close()
    return xyz


def to_psi4_input(fragment, molecule_name, memory='500 mb', command='gradient', symmetry='C1', crank_job='crank_job_1', filename=None):
    """
    Write out psi4 input for crank-launch input
    Parameters
    ----------
    fragment: dict
        JSON crank jobs for fragment
    molecule_name: str
        molecule name
    memory: str
        memory for psi4 job
    command: str
        gradient or optimize. Defualt is gradient. Must be gradient if using geomeTRIC for geometry optimization
    symmetry: str
        point group symmetry. Default is C1
    crank_job: str
        key for crank job in fragment dictionary. Default is crank_job_1
    filename: str
        Name for psi4 input file. Default is None. If None, function will return input string.

    Returns
    -------
    psi4_input: str

    """
    job = fragment['crank_torsion_drives'][crank_job]
    psi4_input = ""

    if memory is not None:
        psi4_input += "memory {}\n".format(memory)

    psi4_input += "\nmolecule {}".format(molecule_name)
    psi4_input += " {\n"
    psi4_input += "symmetry {}\n".format(symmetry)
    charge = fragment['molecule']['molecular_charge']
    multiplicity = fragment['molecule']['molecular_multiplicity']
    psi4_input += "{}  {} \n".format(charge, multiplicity)

    # Convert 1-D Geom to 2-D geometry
    n_atoms = len(fragment['molecule']['symbols'])
    coords_3d = np.array(fragment['molecule']['geometry'], dtype=float).reshape(n_atoms, 3)

    for element, coord in zip(fragment['molecule']['symbols'], coords_3d):
        psi4_input += "%-7s %13.7f %13.7f %13.7f\n" % (element, coord[0], coord[1], coord[2])

    psi4_input += "units Angstrom\n"
    psi4_input += "}\n"

    psi4_input += "\nset {\n"
    psi4_input += " basis {}\n ".format(job['crank_specs']['model']['basis'])
    try:
        options = job['crank_specs']['options']
        for option in options:
            psi4_input += "{} {}\n".format(option, options[option])
    except KeyError:
        # No options were given
        logger().info('no options found')
        pass

    psi4_input += "}\n"

    psi4_input += "\n{}('{}')\n".format(command, job['crank_specs']['model']['method'])

    if filename:
        f = open(filename, 'w')
        f.write(psi4_input)
        f.close()
    else:
        return psi4_input


def to_dihedraltxt(fragment, crank_job='crank_job_1', filename=None):
    """
    Generate dihedral txt file defining dihedrals that should be restraint in crank job

    Parameters
    ----------
    fragment: dict
        JSON crank jobs for fragment
    crank_job: str
        key for crank job to run
    filename: str
        filename of dihedral file. If None, returns input string.

    Returns
    -------
    dih_input: str

    """
    needed_torsions = fragment['needed_torsion_drives']
    crank_jobs = fragment['crank_torsion_drives'][crank_job]
    dihedrals = []
    for torsion in needed_torsions:
        if torsion in crank_jobs:
            dihedrals.append([i-1 for i in needed_torsions[torsion]])

    dih_input = ""
    dih_input += "# dihedral definition by atom indices starting from 0\n"
    dih_input += "# i     j     k     l\n"
    for dihedral in dihedrals:
        dih_input += "  {}     {}     {}     {}\n".format(dihedral[0], dihedral[1], dihedral[2], dihedral[3])
    if filename:
        f = open(filename, 'w')
        f.write(dih_input)
        f.close()
    else:
        return dih_input


def to_mapped_QC_JSON_geometry(molecule, atom_map, charge=0, multiplicity=1):
    """
    Generate xyz coordinates for molecule in the order given by the atom_map. atom_map is a dictionary that maps the
    tag on the SMILES to the atom idex in OEMol.
    Parameters
    ----------
    molecule: OEMol with conformers
    atom_map: dict
        maps tag in SMILES to atom index
    charge: int
        charge of molecule. Default is 0 (neural)
    multiplicity: int
        spin multiplicity of molecule (2S +1). Default is 1 (all electrons are paired)

    Returns
    -------
    dict: QC_JSON Molecule spec {symbols: [], geometry: [], 'molecular_charge': int, 'molecular_multiplicity': int}

    """
    symbols = []
    geometry = []
    if molecule.GetMaxConfIdx() != 1:
        raise Warning("The molecule must have at least and at most 1 conformation")

    # Check if coordinates are missing
    if not has_conformer(molecule, check_two_dimension=True):
        raise RuntimeError("molecule {} does not have a 3D conformation".format(oechem.OEMolToSmiles(molecule)))

    for mapping in range(1, len(atom_map)+1):
        idx = atom_map[mapping]
        atom = molecule.GetAtom(oechem.OEHasAtomIdx(idx))
        syb = oechem.OEGetAtomicSymbol(atom.GetAtomicNum())
        symbols.append(syb)
        for i in range(3):
            geometry.append(molecule.GetCoords()[idx][i])

    return {'symbols': symbols, 'geometry': geometry, 'molecular_charge': charge,
            'molecular_multiplicity': multiplicity}


def bond_order_from_psi4_raw_output(psi_output):
    """
    Extract Wiberg and Mayer bond order from raw psi4 output

    Parameters
    ----------
    psi_output: str
        psi4 raw output. This can be extracted from JSON_data['raw_output'] or by reading entire psi4 output
        file

    Returns
    -------
    bond_order_arrays: dict of numpy arrays
        {Wiberg_psi4: np.array, Mayer_psi4: np.array}
        N x N array. N is the number of atoms in the molecule. Indices in this array corresponds to tag in
        the molecule given by `tagged_smiles` in QC_JSON spec

    """
    size = None
    Mayer = []
    Wiberg = []
    FLAG = None
    for line in psi_output.split('\n'):
        if not line:
            continue
        if 'Wiberg' in line:
            FLAG = 'Wiberg'
        if 'Mayer' in line:
            FLAG = 'Mayer'
        if 'Size' in line:
            size = line.split()
            size = int(size[3]), int(size[-1])
        if FLAG is 'Mayer':
            Mayer.append(line.split())
        if FLAG is 'Wiberg':
            Wiberg.append(line.split())
    if not size:
        raise Warning("Wiberg and Mayer bond orders were not found")
    Wiberg_array = np.zeros(size)
    Mayer_array = np.zeros(size)

    for i, lines in enumerate(zip(Wiberg[2:], Mayer[2:])):
        line_w = lines[0]
        line_m = lines[1]
        if i == 0:
            elements = line_w
            continue
        if not i%float(size[0]+1) and i<((float(size[0]+1))*float((size[0]/5))):
            if len(line_w) !=5:
                if str(size[0]) in line_w:
                    elements = line_w
                continue
            elements = line_w
            continue
        j = line_w[0]
        for k, bo in enumerate(zip(line_w[1:], line_m[1:])):
            bo_w = bo[0]
            bo_m = bo[1]
            try:
                Wiberg_array[int(elements[k])-1][int(j)-1] = bo_w
                Mayer_array[int(elements[k])-1][int(j)-1] = bo_m

            except (ValueError, IndexError):
                pass

    return {'Wiberg_psi4': Wiberg_array, 'Mayer_psi4': Mayer_array}


def bond_order_to_bond_graph(bond_order, threshold=0.8, hydrogen_bond=True, molecule=None, atom_map=None):
    """
    Get bond graph from bond orders. This function returns a set of bonds where the bond order is above a threshold
    Parameters
    ----------
    bond_order: np array
    threshold: int

    Returns
    -------
    bonds: set

    """
    bonds = set()
    for i in range(bond_order.shape[0]):
        for j in range(bond_order.shape[1]):
            if bond_order[i, j] >= threshold:
                if not hydrogen_bond:
                    idx_1 = atom_map[i+1]
                    idx_2 = atom_map[j+1]
                    atom_1 = molecule.GetAtom(oechem.OEHasMapIdx(idx_1))
                    atom_2 = molecule.GetAtom(oechem.OEHasAtomIdx(idx_2))
                    if atom_1.IsHydrogen() or atom_2.IsHydrogen():
                        continue
                if (j+1, i+1) in bonds:
                    continue
                bonds.add((i+1, j+1))
    return bonds


def boltzman_average_bond_order(bond_orders):
    """
    Calculate the Boltzmann weighted bond order average.

    Parameters
    ----------
    bond_orders: Dictionary of bond orders. The key is the energy of the molecule.

    Returns
    -------
    bond_order_arrays: Dictionary of Boltzmann weighted bond orders.

    """
    energies = np.asarray(list(bond_orders.keys()))
    weights = np.exp(-energies/298.15)
    denominator = weights.sum()
    weights = weights/denominator

    Wiberg = np.zeros((tuple([len(energies)]) + bond_orders[energies[0]]['Wiberg_psi4'].shape))
    Mayer = np.zeros((tuple([len(energies)]) + bond_orders[energies[0]]['Wiberg_psi4'].shape))
    for i, energy in enumerate(energies):
        Wiberg[i] = bond_orders[energy]['Wiberg_psi4']
        Mayer[i] = bond_orders[energy]['Mayer_psi4']

    average_Wiberg =( weights[:, np.newaxis, np.newaxis] * Wiberg).sum(axis=0)
    average_Mayer = (weights[:, np.newaxis, np.newaxis] * Mayer).sum(axis=0)
    bond_order_arrays = {'Wiberg_psi4': average_Wiberg, 'Mayer_psi4': average_Mayer}
    return bond_order_arrays


def tag_conjugated_bond(mol, tautomers=None, tag=None, threshold=1.05):
    """
    Add conjugated bond data tag. If the bond order is above the threshold, this tag will be True, otherwise it will be
    False
    Parameters
    ----------
    mol: OEMol
    tautomers: list of OEMols, optional, Default None
        If a list is provided, the conjugation tag will be true if the bond in any of the set of molecules is double.
        The list should consist of resonance structures of the molecules. You can get that from oequacpac.EnumerateTautomres
    tag: str, optional, Default None
        If provided, will use that bond order. Options are WibergBondOrder, Wiberg_psi4, Mayer_psi4
    threshold: int, optional, Default is 1.05
        The fractional bond order threshold above which the bond will be considered conjugated.

    Returns
    -------
    atom_indices: list of atom indices in conjugated bonds.

    """
    atom_indices = []
    for bond in mol.GetBonds():
        resonance = False
        if tautomers is not None:
            for tmol in tautomers:
                t_bond = tmol.GetBond(oechem.OEHasBondIdx(bond.GetIdx()))
                if t_bond.GetOrder() > 1:
                    resonance = True
                    break
        elif bond.GetData()[tag] >= threshold:
            resonance = True
        # Add tag to bond
        conj_tag = oechem.OEGetTag("conjugated")
        bond.SetData(conj_tag, resonance)
        if resonance:
            a1 = bond.GetBgnIdx()
            a2 = bond.GetEndIdx()
            atom_indices.extend([a1, a2])
    return atom_indices


def highlight_bonds(mol_copy, fname, conjugation=True, rotor=False, width=600, height=400, label=None):
    """
    Generate image of molecule with highlighted bonds. The bonds can either be highlighted with a conjugation tag
    or if it is rotatable.

    Parameters
    ----------
    mol_copy: OEMol
    fname: str
        Name of image file
    conjugation: Bool, optional, Default is True
        If True, the bonds with conjugation tag set to True will be highlighted
    rotor: Bool, optional, Default is False
        If True, the rotatable bonds will be highlighted.
    width: int
    height: int
    label: string. Optional, Default is None
        The bond order label. The options are WibergBondOrder, Wiberg_psi4, Mayer_psi4.

    """
    mol = oechem.OEMol(mol_copy)
    bond_index_list = []
    for bond in mol.GetBonds():
        if conjugation:
            try:
                if bond.GetData('conjugated'):
                    bond_index_list.append(bond.GetIdx())
            except ValueError:
                pass
        if rotor:
            if bond.IsRotor():
                bond_index_list.append(bond.GetIdx())

    atomBondSet = oechem.OEAtomBondSet()
    for bond in mol.GetBonds():
        if bond.GetIdx() in bond_index_list:
            atomBondSet.AddBond(bond)
            atomBondSet.AddAtom(bond.GetBgn())
            atomBondSet.AddAtom(bond.GetEnd())

    dopt = oedepict.OEPrepareDepictionOptions()
    dopt.SetSuppressHydrogens(True)
    oedepict.OEPrepareDepiction(mol, dopt)

    opts = oedepict.OE2DMolDisplayOptions(width, height, oedepict.OEScale_AutoScale)
    opts.SetTitleLocation(oedepict.OETitleLocation_Hidden)
    if label is not None:
        bond_label = {'WibergBondOrder': LabelWibergBondOrder, 'Wiberg_psi4': LabelWibergPsiBondOrder,
                      'Mayer_psi4': LabelMayerPsiBondOrder}

        bondlabel = bond_label[label]
        opts.SetBondPropertyFunctor(bondlabel())

    disp = oedepict.OE2DMolDisplay(mol, opts)

    aroStyle = oedepict.OEHighlightStyle_Color
    aroColor = oechem.OEColor(oechem.OEBlack)
    oedepict.OEAddHighlighting(disp, aroColor, aroStyle,
                               oechem.OEIsAromaticAtom(), oechem.OEIsAromaticBond() )
    hstyle = oedepict.OEHighlightStyle_BallAndStick
    hcolor = oechem.OEColor(oechem.OELightBlue)
    oedepict.OEAddHighlighting(disp, hcolor, hstyle, atomBondSet)

    return oedepict.OERenderMolecule(fname, disp)


def bond_order_tag(molecule, atom_map, bond_order_array):
    """
    Add psi bond order to bond in molecule. This function adds a tag to the GetData dictionary
    in bond.GetData()

    Parameters
    ----------
    molecule: OEMol
        This molecule must have tags that corresponds to the atom_map
    atom_map: dict
        dictionary that maps atom tag to atom index
    bond_order_array: dict
        maps Wiberg and Meyer bond indices to N x N numpy arrays.
        N - atoms in molecule. This array contains the bond order for bond(i,j) where i,j correspond to
        tag on atom and index in bond_order_array
    """
    wiberg_bond_order = bond_order_array['Wiberg_psi4']
    mayer_bond_order = bond_order_array['Mayer_psi4']
    # Sanity check, both arrays are same shape
    for i, j in itertools.combinations(range(wiberg_bond_order.shape[0]), 2):
        idx_1 = atom_map[i+1]
        idx_2 = atom_map[j+1]
        atom_1 = molecule.GetAtom(oechem.OEHasAtomIdx(idx_1))
        atom_2 = molecule.GetAtom(oechem.OEHasAtomIdx(idx_2))
        bond = molecule.GetBond(atom_1, atom_2)
        if bond:
            wbo = wiberg_bond_order[i][j]
            mbo = mayer_bond_order[i][j]
            tag = oechem.OEGetTag('Wiberg_psi4')
            bond.SetData(tag, wbo)
            tag = oechem.OEGetTag('Mayer_psi4')
            bond.SetData(tag, mbo)


def mol_to_graph(molecule):
    """
    Generate networkx graph from oe molecule
    """
    import networkx as nx
    G = nx.Graph()
    for atom in molecule.GetAtoms():
        G.add_node(atom.GetIdx(), element=atom.GetElement())
    for bond in molecule.GetBonds():
        G.add_edge(bond.GetBgnIdx(), bond.GetEndIdx(), index=bond.GetIdx())
    return G


def get_charge(molecule):

    charge=0
    for atom in molecule.GetAtoms():
        charge += atom.GetFormalCharge()
    return charge


def png_atoms_labeled(smiles, fname, width=600, height=400, label_scale=2.0, scale_bondwidth=True):
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

    opts = oedepict.OE2DMolDisplayOptions(width, height, oedepict.OEScale_AutoScale)
    opts.SetAtomPropertyFunctor(oedepict.OEDisplayAtomMapIdx())
    opts.SetAtomPropLabelFont(oedepict.OEFont(oechem.OEDarkGreen))
    opts.SetAtomPropLabelFontScale(label_scale)
    opts.SetBondWidthScaling(scale_bondwidth)

    disp = oedepict.OE2DMolDisplay(mol, opts)
    return oedepict.OERenderMolecule(fname, disp)


def png_bond_labels(mol, fname, width=600, height=400, label='WibergBondOrder'):
    """
    Generate png figure of molecule. Bonds should include bond order defined in label

    Parameters
    ----------
    mol: OpenEye OEMol
    fname: str
        filename for png
    width: int
    height: int
    label: str
        Which label to print. Options are WibergBondOrder, Wiberg_psi4 and Mayer_psi4

    Returns
    -------
    bool:
    """

    oedepict.OEPrepareDepiction(mol)


    opts = oedepict.OE2DMolDisplayOptions(width, height, oedepict.OEScale_AutoScale)
    # opts.SetAtomPropertyFunctor(oedepict.OEDisplayAtomIdx())
    # opts.SetAtomPropLabelFont(oedepict.OEFont(oechem.OEDarkGreen))
    bond_label = {'WibergBondOrder': LabelWibergBondOrder, 'Wiberg_psi4': LabelWibergPsiBondOrder,
                  'Mayer_psi4': LabelMayerPsiBondOrder}
    bondlabel = bond_label[label]
    opts.SetBondPropertyFunctor(bondlabel())

    disp = oedepict.OE2DMolDisplay(mol, opts)
    return oedepict.OERenderMolecule(fname, disp)


class LabelWibergBondOrder(oedepict.OEDisplayBondPropBase):
    def __init__(self):
        oedepict.OEDisplayBondPropBase.__init__(self)

    def __call__(self, bond):
        bondOrder = bond.GetData('WibergBondOrder')
        label = "{:.2f}".format(bondOrder)
        return label

    def CreateCopy(self):
        copy = LabelWibergBondOrder()
        return copy.__disown__()


class LabelWibergPsiBondOrder(oedepict.OEDisplayBondPropBase):
    def __init__(self):
        oedepict.OEDisplayBondPropBase.__init__(self)

    def __call__(self, bond):
        bondOrder = bond.GetData('Wiberg_psi4')
        label = "{:.2f}".format(bondOrder)
        return label

    def CreateCopy(self):
        copy = LabelWibergPsiBondOrder()
        return copy.__disown__()


class LabelMayerPsiBondOrder(oedepict.OEDisplayBondPropBase):
    def __init__(self):
        oedepict.OEDisplayBondPropBase.__init__(self)

    def __call__(self, bond):
        bondOrder = bond.GetData('Mayer_psi4')
        label = "{:.2f}".format(bondOrder)
        return label

    def CreateCopy(self):
        copy = LabelMayerPsiBondOrder()
        return copy.__disown__()


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