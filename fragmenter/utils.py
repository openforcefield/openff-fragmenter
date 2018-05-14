import os
import time
import logging
import sys
import json
from uuid import UUID
import numpy as np
import itertools
import re
import codecs

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
    path: str
        path to output file
    base: str
        base name for output file

    """

    outf = open(filename, 'w')
    smiles_list = map(lambda x: x+"\n", list(smiles))
    outf.writelines(smiles_list)
    outf.close()
    if return_fname:
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


def create_mapped_smiles(molecule, tagged=True, explicit_h=True, isomeric=True):
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
    if not explicit_h and not tagged:
        return oechem.OEMolToSmiles(molecule)
    oechem.OEAddExplicitHydrogens(molecule)
    if not tagged and explicit_h and isomeric:
        return oechem.OECreateSmiString(molecule, oechem.OESMILESFlag_Hydrogens | oechem.OESMILESFlag_Isotopes | oechem.OESMILESFlag_AtomStereo
                                        | oechem.OESMILESFlag_BondStereo | oechem.OESMILESFlag_Canonical | oechem.OESMILESFlag_RGroups)
    if not tagged and explicit_h and not isomeric:
        return oechem.OECreateSmiString(molecule, oechem.OESMILESFlag_Hydrogens | oechem.OESMILESFlag_Canonical |
                                        oechem.OESMILESFlag_RGroups)

    for atom in molecule.GetAtoms():
        atom.SetMapIdx(atom.GetIdx() + 1)

    # add tag to data
    tag = oechem.OEGetTag("has_map")
    molecule.SetData(tag, bool(True))

    return oechem.OEMolToSmiles(molecule)


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


def get_atom_map(tagged_smiles, molecule=None, is_mapped=False, StrictStereo=True):
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
            molecule = openeye.generate_conformers(molecule, max_confs=1, strictStereo=StrictStereo)
            # Omega can change the ordering so whatever map existed is not there anymore

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
        print("I have no idea what's happening")
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
            tag = oechem.OEGetTag('Wiberg_psi')
            bond.SetData(tag, wbo)
            tag = oechem.OEGetTag('Mayer_psi')
            bond.SetData(tag, mbo)


# def to_bond_order_dict(molecule, atom_map, bond_order_arrays):
#     """
#
#     Parameters
#     ----------
#     molecule
#     atom_map
#     bond_oder_arrays
#
#     Returns
#     -------
#
#     """
#     # for bond in molecule.GetBonds():
#     #     try:
#     #         bond.GetData('Wiberg_ensa')
#     #     except ValueError:
#     #
#     #     tag = oechem.OEGetTag('Wiberg_ensamble')
#     #     bond.SetData(tag, [])
#     #     tag = oechem.OEGetTag('Mayer_ensamble')
#     #     bond.SetData(tag, [])
#
#     wiberg_bond_order = bond_order_arrays['wiberg_psi']
#     mayer_bond_order = bond_order_arrays['Mayer_psi']
#     for i, j in itertools.combinations(range(wiberg_bond_order.shape[0], 2)):
#         idx_1 = atom_map[i+1]
#         idx_2 = atom_map[j+1]
#         atom_1 = molecule.GetAtom(oechem.OEHasAtomIdx(idx_1))
#         atom_2 = molecule.GetAtom(oechem.OEHasAtomIdx(idx_2))
#         bond = molecule.GetBond(atom_1, atom_2)
#         if bond:
#             wbo = wiberg_bond_order[i][j]
#             mbo = mayer_bond_order[i][j]





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
        bondOrder = bond.GetData('Wiberg_psi')
        label = "{:.2f}".format(bondOrder)
        return label

    def CreateCopy(self):
        copy = LabelWibergPsiBondOrder()
        return copy.__disown__()


class LabelMayerPsiBondOrder(oedepict.OEDisplayBondPropBase):
    def __init__(self):
        oedepict.OEDisplayBondPropBase.__init__(self)

    def __call__(self, bond):
        bondOrder = bond.GetData('Mayer_psi')
        label = "{:.2f}".format(bondOrder)
        return label

    def CreateCopy(self):
        copy = LabelMayerPsiBondOrder()
        return copy.__disown__()

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


def make_python_identifier(string, namespace=None, reserved_words=None,
                           convert='drop', handle='force'):
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

    # create a working copy (and make it lowercase, while we're at it)
    s = string.lower()

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