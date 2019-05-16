"""functions to manipulate, read and write OpenEye and Psi4 molecules"""

try:
    from openeye import oechem, oeomega, oeiupac, oedepict, oequacpac, oeszybki
except ImportError:
    raise Warning("Need license for OpenEye!")
from rdkit import Chem

import cmiles
from .utils import logger, ANGSROM_2_BOHR, BOHR_2_ANGSTROM

import os
import numpy as np
import time
import itertools
import copy
from math import radians

"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Functions to charge, generate conformers and
manipulate OpenEye molecules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""


def get_charges(molecule, max_confs=800, strict_stereo=True,
                normalize=True, keep_confs=None, legacy=True):
    """Generate charges for an OpenEye OEMol molecule.
    Parameters
    ----------
    molecule : OEMol
        Molecule for which to generate conformers.
        Omega will be used to generate max_confs conformations.
    max_confs : int, optional, default=800
        Max number of conformers to generate
    strictStereo : bool, optional, default=True
        If False, permits smiles strings with unspecified stereochemistry.
        See https://docs.eyesopen.com/omega/usage.html
    normalize : bool, optional, default=True
        If True, normalize the molecule by checking aromaticity, adding
        explicit hydrogens, and renaming by IUPAC name.
    keep_confs : int, optional, default=None
        If None, apply the charges to the provided conformation and return
        this conformation, unless no conformation is present.
        Otherwise, return some or all of the generated
        conformations. If -1, all generated conformations are returned.
        Otherwise, keep_confs = N will return an OEMol with up to N
        generated conformations.  Multiple conformations are still used to
        *determine* the charges.
    legacy : bool, default=True
        If False, uses the new OpenEye charging engine.
        See https://docs.eyesopen.com/toolkits/python/quacpactk/OEProtonFunctions/OEAssignCharges.html#
    Returns
    -------
    charged_copy : OEMol
        A molecule with OpenEye's recommended AM1BCC charge selection scheme.
    Notes
    -----
    Roughly follows
    http://docs.eyesopen.com/toolkits/cookbook/python/modeling/am1-bcc.html
    """

    # If there is no geometry, return at least one conformation.
    if molecule.GetConfs() == 0:
        keep_confs = 1

    if not oechem.OEChemIsLicensed(): raise(ImportError("Need License for OEChem!"))
    if not oequacpac.OEQuacPacIsLicensed(): raise(ImportError("Need License for oequacpac!"))

    if normalize:
        molecule = normalize_molecule(molecule)
    else:
        molecule = oechem.OEMol(molecule)

    charged_copy = generate_conformers(molecule, max_confs=max_confs, strict_stereo=strict_stereo)  # Generate up to max_confs conformers

    if not legacy:
        # 2017.2.1 OEToolkits new charging function
        status = oequacpac.OEAssignCharges(charged_copy, oequacpac.OEAM1BCCCharges())
        if not status: raise(RuntimeError("OEAssignCharges failed."))
    else:
        # AM1BCCSym recommended by Chris Bayly to KAB+JDC, Oct. 20 2014.
        status = oequacpac.OEAssignPartialCharges(charged_copy, oequacpac.OECharges_AM1BCCSym)
        if not status: raise(RuntimeError("OEAssignPartialCharges returned error code %d" % status))



    #Determine conformations to return
    if keep_confs == None:
        #If returning original conformation
        original = molecule.GetCoords()
        #Delete conformers over 1
        for k, conf in enumerate( charged_copy.GetConfs() ):
            if k > 0:
                charged_copy.DeleteConf(conf)
        #Copy coordinates to single conformer
        charged_copy.SetCoords( original )
    elif keep_confs > 0:
        logger().debug("keep_confs was set to %s. Molecule positions will be reset." % keep_confs)

        #Otherwise if a number is provided, return this many confs if available
        for k, conf in enumerate( charged_copy.GetConfs() ):
            if k > keep_confs - 1:
                charged_copy.DeleteConf(conf)
    elif keep_confs == -1:
        #If we want all conformations, continue
        pass
    else:
        #Not a valid option to keep_confs
        raise(ValueError('Not a valid option to keep_confs in get_charges.'))

    return charged_copy


def generate_conformers(molecule, max_confs=800, dense=False, strict_stereo=True, ewindow=15.0, rms_threshold=1.0, strict_types=True,
                        can_order=True, copy=True):
    """Generate conformations for the supplied molecule
    Parameters
    ----------
    molecule : OEMol
        Molecule for which to generate conformers
    max_confs : int, optional, default=800
        Max number of conformers to generate.  If None, use default OE Value.
    strict_stereo : bool, optional, default=True
        If False, permits smiles strings with unspecified stereochemistry.
    strict_types : bool, optional, default=True
        If True, requires that Omega have exact MMFF types for atoms in molecule; otherwise, allows the closest atom
        type of the same element to be used.
    Returns
    -------
    molcopy : OEMol
        A multi-conformer molecule with up to max_confs conformers.
    Notes
    -----
    Roughly follows
    http://docs.eyesopen.com/toolkits/cookbook/python/modeling/am1-bcc.html
    """
    if copy:
        molcopy = oechem.OEMol(molecule)
    else:
        molcopy = molecule

    if dense:
        omega_opts = oeomega.OEOmegaOptions(oeomega.OEOmegaSampling_Dense)
        omega = oeomega.OEOmega(omega_opts)
    else:
        omega = oeomega.OEOmega()

    if cmiles.utils.has_atom_map(molcopy):
        remove_map(molcopy)

    # These parameters were chosen to match http://docs.eyesopen.com/toolkits/cookbook/python/modeling/am1-bcc.html
    omega.SetMaxConfs(max_confs)
    omega.SetIncludeInput(True)
    omega.SetCanonOrder(can_order)

    omega.SetSampleHydrogens(True)  # Word to the wise: skipping this step can lead to significantly different charges!
    omega.SetEnergyWindow(ewindow)
    omega.SetRMSThreshold(rms_threshold)  # Word to the wise: skipping this step can lead to significantly different charges!

    omega.SetStrictStereo(strict_stereo)
    omega.SetStrictAtomTypes(strict_types)

    omega.SetIncludeInput(False)  # don't include input
    if max_confs is not None:
        omega.SetMaxConfs(max_confs)

    status = omega(molcopy)  # generate conformation
    if not status:
        raise(RuntimeError("omega returned error code %d" % status))

    restore_map(molcopy)

    return molcopy


def generate_grid_conformers(molecule, dihedrals, intervals, max_rotation=360, copy_mol=True):
    """
    Generate conformers using torsion angle grids.

    Parameters
    ----------
    molecule: OEMol
    dihedrals: list of
    intervals

    Returns
    -------

    """
    # molecule must be mapped
    if copy_mol:
        molecule = copy.deepcopy(molecule)
    if cmiles.utils.has_atom_map(molecule):
        remove_map(molecule)
    else:
        raise ValueError("Molecule must have map indices")

    # Check length of dihedrals match length of intervals

    conf_mol = generate_conformers(molecule, max_confs=1)
    conf = conf_mol.GetConfs().next()
    coords = oechem.OEFloatArray(conf.GetMaxAtomIdx()*3)
    conf.GetCoords(coords)

    torsions = [[conf_mol.GetAtom(oechem.OEHasMapIdx(i+1)) for i in dih] for dih in dihedrals]

    for i, tor in enumerate(torsions):
        copy_conf_mol = copy.deepcopy(conf_mol)
        conf_mol.DeleteConfs()
        for conf in copy_conf_mol.GetConfs():
            coords = oechem.OEFloatArray(conf.GetMaxAtomIdx()*3)
            conf.GetCoords(coords)
            for angle in range(5, max_rotation+5, intervals[i]):
                newconf = conf_mol.NewConf(coords)
                oechem.OESetTorsion(newconf, tor[0], tor[1], tor[2], tor[3], radians(angle))

    restore_map(conf_mol)
    return conf_mol


def resolve_clashes(mol):
    """
    Taken from quanformer (https://github.com/MobleyLab/quanformer/blob/master/quanformer/initialize_confs.py#L54
    Minimize conformers with severe steric interaction.
    Parameters
    ----------
    mol : single OEChem molecule (single conformer)
    clashfile : string
        name of file to write output
    Returns
    -------
    boolean
        True if completed successfully, False otherwise.
    """

    # set general energy options along with the single-point specification
    spSzybki = oeszybki.OESzybkiOptions()
    spSzybki.SetForceFieldType(oeszybki.OEForceFieldType_MMFF94S)
    spSzybki.SetSolventModel(oeszybki.OESolventModel_Sheffield)
    spSzybki.SetRunType(oeszybki.OERunType_SinglePoint)

    # generate the szybki MMFF94 engine for single points
    szSP = oeszybki.OESzybki(spSzybki)

    # construct minimiz options from single-points options to get general optns
    optSzybki = oeszybki.OESzybkiOptions(spSzybki)

    # now reset the option for minimization
    optSzybki.SetRunType(oeszybki.OERunType_CartesiansOpt)

    # generate szybki MMFF94 engine for minimization
    szOpt = oeszybki.OESzybki(optSzybki)
    # add strong harmonic restraints to nonHs

    szOpt.SetHarmonicConstraints(10.0)
    # construct a results object to contain the results of a szybki calculation

    szResults = oeszybki.OESzybkiResults()
    # work on a copy of the molecule
    tmpmol = oechem.OEMol(mol)
    if not szSP(tmpmol, szResults):
        print('szybki run failed for %s' % tmpmol.GetTitle())
        return False
    Etotsp = szResults.GetTotalEnergy()
    Evdwsp = szResults.GetEnergyTerm(oeszybki.OEPotentialTerms_MMFFVdW)
    if Evdwsp > 35:
        if not szOpt(tmpmol, szResults):
            print('szybki run failed for %s' % tmpmol.GetTitle())
            return False
        Etot = szResults.GetTotalEnergy()
        Evdw = szResults.GetEnergyTerm(oeszybki.OEPotentialTerms_MMFFVdW)
        # wfile = open(clashfile, 'a')
        # wfile.write(
        #     '%s resolved bad clash: initial vdW: %.4f ; '
        #     'resolved EvdW: %.4f\n' % (tmpmol.GetTitle(), Evdwsp, Evdw))
        # wfile.close()
        mol.SetCoords(tmpmol.GetCoords())
    oechem.OESetSDData(mol, oechem.OESDDataPair('MM Szybki Single Point Energy', "%.12f" % szResults.GetTotalEnergy()))
    return True


def remove_clash(multi_conformer, in_place=True):
    # resolve clashes
    if not in_place:
        multi_conformer = copy.deepcopy(multi_conformer)
    confs_to_remove = []
    for conf in multi_conformer.GetConfs():
        if not resolve_clashes(conf):
            confs_to_remove.append(conf)
    for i in confs_to_remove:
        multi_conformer.DeleteConf(i)

    if not in_place:
        return multi_conformer


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


# def mol_to_graph(molecule):
#     """
#     Generate networkx graph from oe molecule
#     """
#     import networkx as nx
#     G = nx.Graph()
#     for atom in molecule.GetAtoms():
#         G.add_node(atom.GetIdx(), element=atom.GetElement())
#     for bond in molecule.GetBonds():
#         G.add_edge(bond.GetBgnIdx(), bond.GetEndIdx(), index=bond.GetIdx())
#     return G


def get_charge(molecule):

    charge = 0
    for atom in molecule.GetAtoms():
        charge += atom.GetFormalCharge()
    return charge


"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Functions for bond orders (Psi4 and OpenEye)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""


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


"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Functions to read and write molecules and SMILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""


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
        filename = os.path.join(os.getcwd(), filename)
        return filename


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


def file_to_oemols(filename, title=True, verbose=False):
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

    molecule = oechem.OEMol()
    while oechem.OEReadMolecule(ifs, molecule):
        molecule_copy = oechem.OEMol(molecule)
        oechem.OEPerceiveChiral(molecule)
        oechem.OE3DToAtomStereo(molecule)
        oechem.OE3DToBondStereo(molecule)
        if title:
            title = molecule_copy.GetTitle()
            if verbose:
                logger().info("Reading molecule {}".format(title))

        mollist.append(normalize_molecule(molecule_copy, title))

    # if len(mollist) <= 1:
    #     mollist = mollist[0]

    ifs.close()

    return mollist


def smifile_to_rdmols(filename):
    """
    Read SMILES file and return list of RDmols

    Parameters
    ----------
    filename: str. Path to file

    Returns
    -------
    rd_mols: list
        list of RDKit molecules

    """
    smiles_txt = open(filename, 'r').read()
    # Check first line
    first_line = smiles_txt.split('\n')[0]
    if first_line != 'SMILES':
        smiles_txt = 'SMILES\n' + smiles_txt

    rd_mol_supp = Chem.SmilesMolSupplierFromText(smiles_txt)
    rd_mols = [x for x in rd_mol_supp]

    # Check for failure to parse
    nones = []
    for i, mol in enumerate(rd_mols):
        if mol is None:
            nones.append(i)

    if len(nones) > 0:
        # Find SMILES that did not parse
        smiles_list = smiles_txt.split('\n')[1:]
        print(nones)
        missing_mols = [smiles_list[none] for none in nones]
        lines = [int(none) + 1 for none in nones]
        error = RuntimeError("Not all SMILES were parsed properly. {} indices are None in the rd_mols list. The corresponding"
                           "SMILES are {}. They are on lines {} in the file ".format(nones, missing_mols, lines))
        error.results = rd_mols
        raise error

    return rd_mols


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
        SMILES.append(cmiles.utils.mol_to_smiles(mol, mapped=False, explicit_hydrogen=False, isomeric=isomeric))

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


def standardize_molecule(molecule, title=''):

    if isinstance(molecule, oechem.OEMol):
        mol = molecule
        return molecule

    if isinstance(molecule, str):
        mol = oechem.OEMol()
        # First try reading as smiles
        if not oechem.OESmilesToMol(mol, molecule):
            # Try reading as input file
            ifs = oechem.oemolistream()
            if not ifs.open(molecule):
                raise Warning('Could not parse molecule.')

        # normalize molecule
        mol = normalize_molecule(mol, title=title)

    else:
        raise TypeError("Wrong type of input for molecule. Can be SMILES, filename or OEMol")
    return mol


# def multiconf_mol_to_qcschema(mapped_mol):
#     """
#
#     """
#     if not cmiles.utils.has_atom_map(mapped_mol):


"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Functions to work with mapped SMILES
These functions are used to keep the orders atoms consistent across different molecule graphs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""


def remove_map(molecule, keep_map_data=True):
    """
    Remove atom map but store it in atom data.
    Parameters
    ----------
    molecule

    Returns
    -------

    """
    for atom in molecule.GetAtoms():
        if atom.GetMapIdx() !=0:
            if keep_map_data:
                atom.SetData('MapIdx', atom.GetMapIdx())
            atom.SetMapIdx(0)


def restore_map(molecule):
    """
    Restore atom map from atom data
    """
    for atom in molecule.GetAtoms():
        if atom.HasData('MapIdx'):
            atom.SetMapIdx(atom.GetData('MapIdx'))


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
        smiles = cmiles.utils.mol_to_smiles(mol, mapped=True, explicit_hydrogen=True, isomeric=True)
        #ToDo:
        #  make sure this still works (probably doesn't because of copy of molecule. better to use list of molecule
        # with name if molecule has title
        oechem.OEWriteMolecule(ofs, mol)


def to_mapped_xyz(molecule, atom_map=None, conformer=None, xyz_format=True, filename=None):
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
    str: elements and xyz coordinates (in angstroms) in order of tagged SMILES

    """
    if not atom_map and not cmiles.utils.has_atom_map(molecule):
        raise ValueError("If molecule does not have atom map, you must provide an atom map")
    if not has_conformer(molecule, check_two_dimension=True):
        raise ValueError("Molecule must have conformers")
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

            for mapping in range(1, molecule.NumAtoms()+1):
                if not atom_map:
                    atom = molecule.GetAtom(oechem.OEHasMapIdx(mapping))
                    idx = atom.GetIdx()
                else:
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
    else:
        return xyz


def get_mapped_connectivity_table(molecule, atom_map=None):
    """
    generate a connectivity table with map indices

    Parameters
    ----------
    mapped_molecule: oemol or string
        A mapped molecule or a mapped SMILES
    Returns
    -------
    connectivity_table: list
        list of list of map indices of bond and order [[map_idx_1, map_idx_2, bond_order] ...]
    """
    # Should I allow mapped SMILES too?
    if isinstance(molecule, str):
        # Input is a SMILES
        molecule = smiles_to_oemol(molecule)
    if isinstance(molecule, oechem.OEMol):
        if not cmiles.utils.has_atom_map(molecule) and atom_map is None:
            raise TypeError("Molecule must contain map indices. You can get this by generating a molecule from a mapped SMILES")

    if atom_map is None:
        connectivity_table = [[bond.GetBgn().GetMapIdx()-1, bond.GetEnd().GetMapIdx()-1, bond.GetOrder()]
                              for bond in molecule.GetBonds()]
    else:
        # First convert mapping from map:idx to idx:map
        inv_map = dict(zip(atom_map.values(), atom_map.keys()))
        connectivity_table = [[inv_map[bond.GetBgnIdx()]-1, inv_map[bond.GetEndIdx()]-1, bond.GetOrder()]
                              for bond in molecule.GetBonds()]
    return connectivity_table


def from_mapped_xyz_to_mol_idx_order(mapped_coords, atom_map):
    """
    """
    # reshape
    mapped_coords = np.array(mapped_coords, dtype=float).reshape(int(len(mapped_coords)/3), 3)
    coords = np.zeros((mapped_coords.shape))
    for m in atom_map:
        coords[atom_map[m]] = mapped_coords[m-1]

    # flatten
    coords = coords.flatten()
    return coords


def qcschema_to_xyz_format(qcschema, name=None, filename=None):
    """
    Write qcschema molecule to xyz format
    Parameters
    ----------
    qcschema: dict
        qcschema molecule. Must have symbols and geometry
    name: str, optional, default None
        name of molecule
    filename: str, optional, default None
        If filename given, write out file to disk. If not, will return xyz string

    Returns
    -------
    xyz: str
        qcschema molecule in xyz format

    """
    if not isinstance(qcschema, list):
        qcschema = [qcschema]
    xyz = ""
    for qcmol in qcschema:
        symbols = qcmol['symbols']
        coords = qcmol['geometry']
        coords = np.asarray(coords)*BOHR_2_ANGSTROM
        xyz += "{}\n".format(len(symbols))
        xyz += "{}\n".format(name)
        for i, s in enumerate(symbols):
            xyz += "  {}      {:05.3f}   {:05.3f}   {:05.3f}\n".format(s,
                                                                      coords[i * 3],
                                                                      coords[i * 3 + 1],
                                                                      coords[i * 3 + 2])
    if filename:
        with open(filename, 'w') as f:
            f.write(xyz)
    else:
        return xyz


def qcschema_to_xyz_traj(final_molecule_grid, filename=None):
    """
    Generate an xyz trajectory from QCArchive final molecule output from torsion drive.
    The input should be the grid for one torsiondrive job. Remember to deserialize the output from QCArchive
    This function assumes a 1D torsion scan
    Parameters
    ----------
    final_molecule_grid: dict
        maps grid id to qcschema molecule
    filename: str, optional, default None
        If a name is given, an xyz trajectory will be written to file. If not, xyz string will be returned
    """
    xyz = ""
    angles = sorted(list(final_molecule_grid.keys()))
    for angle in angles:
        molecule = final_molecule_grid[angle]
        name = str(angle)
        xyz += qcschema_to_xyz_format(molecule, name)
    if filename:
        with open(filename, 'w') as f:
            f.write(xyz)
    else:
        return xyz

"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Functions for molecule visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

_KELLYS_COLORS = ['#ebce2b', '#702c8c', '#db6917', '#96cde6', '#ba1c30', '#c0bd7f', '#7f7e80',
                  '#5fa641', '#d485b2', '#4277b6', '#df8461', '#463397', '#e1a11a', '#91218c', '#e8e948', '#7e1510',
                  '#92ae31', '#6f340d', '#d32b1e', '#2b3514']


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


def highlight_torsion(mapped_molecule, dihedrals, fname, width=600, height=400, combine_central_bond=True, color=None):

    mol = oechem.OEMol(mapped_molecule)

    atom_bond_sets = []

    if combine_central_bond:
        central_bonds = [(tor[1], tor[2]) for tor in dihedrals]
        eq_torsions = {cb : [tor for tor in dihedrals if cb == (tor[1], tor[2]) or cb ==(tor[2], tor[1])] for cb in central_bonds}

        for cb in eq_torsions:
            atom_bond_set = oechem.OEAtomBondSet()
            for dihedral in eq_torsions[cb]:
                a = mol.GetAtom(oechem.OEHasMapIdx(dihedral[0]+1))
                atom_bond_set.AddAtom(a)

                for idx in dihedral[1:]:
                    a2 = mol.GetAtom(oechem.OEHasMapIdx(idx+1))
                    atom_bond_set.AddAtom((a2))
                    bond = mol.GetBond(a, a2)
                    atom_bond_set.AddBond(bond)
                    a = a2
            atom_bond_sets.append(atom_bond_set)

    if not combine_central_bond:
        for dihedral in dihedrals:
            atom_bond_set = oechem.OEAtomBondSet()
            a = mol.GetAtom(oechem.OEHasMapIdx(dihedral[0]+1))
            atom_bond_set.AddAtom(a)

            for idx in dihedral[1:]:
                a2 = mol.GetAtom(oechem.OEHasMapIdx(idx+1))
                atom_bond_set.AddAtom((a2))
                bond = mol.GetBond(a, a2)
                atom_bond_set.AddBond(bond)
                a = a2
            atom_bond_sets.append(atom_bond_set)

    dopt = oedepict.OEPrepareDepictionOptions()
    dopt.SetSuppressHydrogens(False)
    oedepict.OEPrepareDepiction(mol, dopt)

    opts = oedepict.OE2DMolDisplayOptions(width, height, oedepict.OEScale_AutoScale)
    opts.SetTitleLocation(oedepict.OETitleLocation_Hidden)

    disp = oedepict.OE2DMolDisplay(mol, opts)

    aroStyle = oedepict.OEHighlightStyle_Color
    aroColor = oechem.OEColor(oechem.OEBlack)
    oedepict.OEAddHighlighting(disp, aroColor, aroStyle,
                               oechem.OEIsAromaticAtom(), oechem.OEIsAromaticBond() )
    hstyle = oedepict.OEHighlightStyle_BallAndStick

    if color:
        highlight = oechem.OEColor(color)
        # combine all atom_bond_sets
        atom_bond_set = oechem.OEAtomBondSet()
        for ab_set in atom_bond_sets:
            for a in ab_set.GetAtoms():
                atom_bond_set.AddAtom(a)
            for b in ab_set.GetBonds():
                atom_bond_set.AddBond(b)
        oedepict.OEAddHighlighting(disp, highlight, hstyle, atom_bond_set)
    else:
        highlight = oedepict.OEHighlightOverlayByBallAndStick(oechem.OEGetContrastColors())
        oedepict.OEAddHighlightOverlay(disp, highlight, atom_bond_sets)
    #hcolor = oechem.OEColor(oechem.OELightBlue)
    #oedepict.OEAddHighlighting(disp, hcolor, hstyle, atom_bond_sets)

    return oedepict.OERenderMolecule(fname, disp)


def highltigh_torsion_by_cluster(mapped_molecule, clustered_dihedrals, fname, width=600, height=400):
    """
    Highlight torsion by cluster. This is used to visualize clustering output.

    Parameters
    ----------
    mapped_molecule: oemol with map indices
    clustered_dihedrals
    fname
    width
    height

    Returns
    -------

    """
    mol = oechem.OEMol(mapped_molecule)
    atom_bond_sets = []

    for cluster in clustered_dihedrals:
        atom_bond_set = oechem.OEAtomBondSet()
        for dihedral in clustered_dihedrals[cluster]:
            a = mol.GetAtom(oechem.OEHasMapIdx(dihedral[0]+1))
            atom_bond_set.AddAtom(a)
            for idx in dihedral[1:]:
                a2 = mol.GetAtom(oechem.OEHasMapIdx(idx+1))
                atom_bond_set.AddAtom(a2)
                bond = mol.GetBond(a, a2)
                atom_bond_set.AddBond(bond)
                a=a2
        atom_bond_sets.append(atom_bond_set)

    dopt = oedepict.OEPrepareDepictionOptions()
    dopt.SetSuppressHydrogens(False)
    oedepict.OEPrepareDepiction(mol, dopt)

    opts = oedepict.OE2DMolDisplayOptions(width, height, oedepict.OEScale_AutoScale)
    opts.SetTitleLocation(oedepict.OETitleLocation_Hidden)

    disp = oedepict.OE2DMolDisplay(mol, opts)

    aroStyle = oedepict.OEHighlightStyle_Color
    aroColor = oechem.OEColor(oechem.OEBlack)
    oedepict.OEAddHighlighting(disp, aroColor, aroStyle,
                               oechem.OEIsAromaticAtom(), oechem.OEIsAromaticBond() )
    hstyle = oedepict.OEHighlightStyle_BallAndStick

    # if color:
    #     highlight = oechem.OEColor(color)
    #     # combine all atom_bond_sets
    #     atom_bond_set = oechem.OEAtomBondSet()
    #     for ab_set in atom_bond_sets:
    #         for a in ab_set.GetAtoms():
    #             atom_bond_set.AddAtom(a)
    #         for b in ab_set.GetBonds():
    #             atom_bond_set.AddBond(b)
    #     oedepict.OEAddHighlighting(disp, highlight, hstyle, atom_bond_set)
    # else:
    highlight = oedepict.OEHighlightOverlayByBallAndStick(oechem.OEGetContrastColors())
    oedepict.OEAddHighlightOverlay(disp, highlight, atom_bond_sets)
    #hcolor = oechem.OEColor(oechem.OELightBlue)
    #oedepict.OEAddHighlighting(disp, hcolor, hstyle, atom_bond_sets)

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


def png_atoms_labeled(mol, fname, map_idx=True, width=600, height=400, label_scale=2.0, scale_bondwidth=True):
    """Write out png file of molecule with atoms labeled with their map index.

    Parameters
    ----------
    smiles: str
        SMILES
    fname: str
        absolute path and filename for png
    map_idx: bool
        If True, lable atoms with map index instead of atom index. If set to True, input SMILES must have map indices.

    """

    if isinstance(mol, str):
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smiles)
    oedepict.OEPrepareDepiction(mol)
    opts = oedepict.OE2DMolDisplayOptions(width, height, oedepict.OEScale_AutoScale)

    if map_idx:
        # check if molecule has map
        if not cmiles.utils.has_atom_map(mol):
            raise ValueError("Input SMILES must have atom maps to display map indices in image")
        opts.SetAtomPropertyFunctor(oedepict.OEDisplayAtomMapIdx())
        opts.SetAtomPropertyFunctor(oedepict.OEDisplayAtomMapIdx())
    if not map_idx:
        opts.SetAtomPropertyFunctor(oedepict.OEDisplayAtomIdx())

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

def to_pdf(molecules, oname, rows=5, cols=3):
    itf = oechem.OEInterface()
    PageByPage = True

    ropts = oedepict.OEReportOptions(rows, cols)
    ropts.SetHeaderHeight(25)
    ropts.SetFooterHeight(25)
    ropts.SetCellGap(2)
    ropts.SetPageMargins(10)
    report = oedepict.OEReport(ropts)

    cellwidth, cellheight = report.GetCellWidth(), report.GetCellHeight()
    opts = oedepict.OE2DMolDisplayOptions(cellwidth, cellheight, oedepict.OEScale_AutoScale)
    oedepict.OESetup2DMolDisplayOptions(opts, itf)

    for mol in molecules:
        cell = report.NewCell()
        oedepict.OEPrepareDepiction(mol)
        disp = oedepict.OE2DMolDisplay(mol, opts)
        oedepict.OERenderMolecule(cell, disp)
        oedepict.OEDrawCurvedBorder(cell, oedepict.OELightGreyPen, 10.0)

    oedepict.OEWriteReport(oname, report)
