__author__ = 'Chaya D. Stern'

try:
    import openeye.oechem as oechem
except ImportError:
    pass
import warnings
import numpy as np
import os
import time
import json
from openmoltools import openeye

from . import utils

warnings.simplefilter('always')


def fragment_to_torsion_scan(fragments, json_filename=None):
    """
    Wrapper function for finding torsions to drive in fragments and generating crank jobs

    Parameters
    ----------
    fragments: dict
        This dictionary has 2 dictionaries. 1) provenance: Information on how the fragments were generated
                                            2) fragments: Maps molecule SMILES to fragment SMILES

    json_filename: str
        If None, will not save molecules to file.
    Returns
    -------
    molecules: dict
        dictionary that maps fragment SMILES to crank job specs (includes torsions to drive, crank specs and provenance)
    """
    provenance = fragments['provenance']
    molecules = {}
    for parent in fragments['fragments']:
        for frag in fragments['fragments'][parent]:
            json_specs = dict()
            json_specs['provenance'] = provenance
            json_specs['provenance']['parent_molecule'] = parent
            json_specs['canonical_isomeric_SMILES'] = frag
            molecule = openeye.smiles_to_oemol(frag)
            json_specs['canonical_SMILES'] = oechem.OECreateCanSmiString(molecule)
            explicit_h_isomeric = utils.create_mapped_smiles(molecule, tagged=False)
            json_specs['explicit_hydrogen_canonical_isomeric_SMILES'] = explicit_h_isomeric
            explicit_h = utils.create_mapped_smiles(molecule, tagged=False, isomeric=False)
            json_specs['explicit_hydrogen_canonical_SMILES'] = explicit_h
            tagged_SMARTS = utils.create_mapped_smiles(molecule)
            json_specs['tagged_SMARTS'] = tagged_SMARTS
            molecule, atom_map = utils.get_atom_map(tagged_SMARTS, is_mapped=True)
            # Find formal charge
            charge = 0
            for atom in molecule.GetAtoms():
                charge += atom.GetFormalCharge()
            QC_JSON_molecule = utils.to_mapped_QC_JSON_geometry(molecule, atom_map, charge=charge)
            json_specs['molecule'] = QC_JSON_molecule
            needed_torsion_drives = find_torsions(molecule)
            json_specs['needed_torsion_drives'] = needed_torsion_drives
            define_crank_job(json_specs)
            molecules[frag] = json_specs

    if json_filename:
        f = open(json_filename, 'w')
        j = json.dump(molecules, f, indent=4, sort_keys=True, cls=utils.UUIDEncoder)
        f.close()

    return molecules


def find_torsions(molecule):
    """
    This function takes an OEMol (atoms must be tagged with index map) and finds the map indices for torsion that need
    to be driven.
    Parameters
    ----------
    molecule : OEMol
        The atoms in the molecule need to be tagged with map indices

    """
    # Check if molecule has map
    is_mapped = utils.is_mapped(molecule)
    if not is_mapped:
        warnings.warn('Molecule does not have atom map. A new map will be generated. You might need a new tagged SMARTS if the ordering was changed')
        tagged_smiles = utils.create_mapped_smiles(molecule)
        utils.logger().info('If you already have a tagged SMARTS, compare it with the new one to ensure the ordering did not change')
        utils.logger().info('The new tagged SMARTS is: {}'.format(tagged_smiles))


    mid_tors = [[tor.a, tor.b, tor.c, tor.d ] for tor in oechem.OEGetTorsions(molecule)]

    # This smarts should match terminal torsions such as -CH3, -NH2, -NH3+, -OH, and -SH
    smarts = '[*]~[*]-[X2H1,X3H2,X4H3]-[#1]'
    qmol=oechem.OEQMol()
    if not oechem.OEParseSmarts(qmol, smarts):
        warnings.warn('OEParseSmarts failed')
    ss = oechem.OESubSearch(qmol)
    mol = oechem.OEMol(molecule)
    h_tors = []
    oechem.OEPrepareSearch(mol, ss)
    unique = True
    for match in ss.Match(mol, unique):
        tor = []
        for ma in match.GetAtoms():
            tor.append(ma.target)
        h_tors.append(tor)

    # Combine middle and terminal torsions
    all_tors = mid_tors + h_tors
    # Sort all_tors so that it's grouped by central bond
    central_bonds = np.zeros((len(all_tors), 3), dtype=int)
    for i, tor in enumerate(all_tors):
        central_bonds[i][0] = i
        central_bonds[i][1] = tor[1].GetIdx()
        central_bonds[i][2] = tor[2].GetIdx()

    grouped = central_bonds[central_bonds[:, 2].argsort()]
    sorted_tors = [all_tors[i] for i in grouped[:, 0]]

    # Keep only one torsion per rotatable bond
    tors = []
    best_tor = [sorted_tors[0][0], sorted_tors[0][0], sorted_tors[0][0], sorted_tors[0][0]]
    first_pass = True
    for tor in sorted_tors:
        utils.logger().info("Map Idxs: {} {} {} {}".format(tor[0].GetMapIdx(), tor[1].GetMapIdx(), tor[2].GetMapIdx(), tor[3].GetMapIdx()))
        utils.logger().info("Atom Numbers: {} {} {} {}".format(tor[0].GetAtomicNum(), tor[1].GetAtomicNum(), tor[2].GetAtomicNum(), tor[3].GetAtomicNum()))
        if tor[1].GetMapIdx() != best_tor[1].GetMapIdx() or tor[2].GetMapIdx() != best_tor[2].GetMapIdx():
            new_tor = True
            if not first_pass:
                utils.logger().info("Adding to list: {} {} {} {}".format(best_tor[0].GetMapIdx(), best_tor[1].GetMapIdx(), best_tor[2].GetMapIdx(), best_tor[3].GetMapIdx()))
                tors.append(best_tor)
            first_pass = False
            best_tor = tor
            best_tor_order = tor[0].GetAtomicNum() + tor[3].GetAtomicNum()
            utils.logger().info("new_tor with central bond across atoms: {} {}".format(tor[1].GetMapIdx(), tor[2].GetMapIdx()))
        else:
            utils.logger().info("Not a new_tor but now with end atoms: {} {}".format(tor[0].GetMapIdx(), tor[3].GetMapIdx()))
            tor_order = tor[0].GetAtomicNum() + tor[3].GetAtomicNum()
            if tor_order > best_tor_order:
                best_tor = tor
                best_tor_order = tor_order
    utils.logger().info("Adding to list: {} {} {} {}".format(best_tor[0].GetMapIdx(), best_tor[1].GetMapIdx(), best_tor[2].GetMapIdx(), best_tor[3].GetMapIdx()))
    tors.append(best_tor)

    utils.logger().info("List of torsion to drive:")
    for tor in tors:
        utils.logger().info("Idx: {} {} {} {}".format(tor[0].GetMapIdx(), tor[1].GetMapIdx(), tor[2].GetMapIdx(), tor[3].GetMapIdx()))
        utils.logger().info("Atom numbers: {} {} {} {}".format(tor[0].GetAtomicNum(), tor[1].GetAtomicNum(), tor[2].GetAtomicNum(), tor[3].GetAtomicNum()))

    needed_torsion_scans = dict()
    for i, tor in enumerate(tors):
        tor_name = ((tor[0].GetMapIdx()), (tor[1].GetMapIdx()), (tor[2].GetMapIdx()), (tor[3].GetMapIdx()))
        needed_torsion_scans['torsion_{}'.format(str(i))] = tor_name
    return needed_torsion_scans


def define_crank_job(fragment_data, grid=None, combinations=None, qc_program='Psi4', method='B3LYP', basis='aug-cc-pVDZ', **kwargs):
    """

    Parameters
    ----------
    fragment_data
    grid: list of int
        spacing for dihedral scan grid points in degree. Must be divisible by 360. If only one value is given, all
        dimensions will have the same grid spacing.
        Default None. When none, all dimensions will have a grid spacing of 30 degrees.
        Can also be a list of lists for combinations.
    combinations: list of ints
        can be list of list. Each list defines which torsions in needed_torsion_drives to run on crank.
        Default is None. Only one crank job will be defined for all torsions in needed_torsion_drives
    qc_program
    method
    kwargs

    Returns
    -------

    """
    if not combinations:
        #Only one crank job including all needed torsions
        needed_torsion_drives = fragment_data['needed_torsion_drives']
        scan_dimension = len(needed_torsion_drives)
        if not grid:
            grid = [30]*scan_dimension
        if type(grid) is not list:
            # Evenly spaced grid for all dimensions
            grid = [grid]*scan_dimension
        grid_dimension = len(grid)
        if grid_dimension != scan_dimension:
                raise Exception("scan dimension {} must be equal to grid dimension {}".format(scan_dimension, grid_dimension))
        # Check that grid is divisible by 360
        for spacing in grid:
            if 360 % spacing:
                raise ValueError("grid spacing must be a factor of 360")

        fragment_data['crank_torsion_drives'] = dict()
        fragment_data['crank_torsion_drives']['crank_job_1'] = dict()
        fragment_data['crank_torsion_drives']['crank_job_1']['crank_specs'] = dict()

        model = {'method': method, 'basis': basis}
        options = {'scf_type': 'df'}
        if kwargs:
            options = kwargs['options']
        fragment_data['crank_torsion_drives']['crank_job_1']['crank_specs']['model'] = model

        fragment_data['crank_torsion_drives']['crank_job_1']['crank_specs']['options'] = options

        for d, spacing in enumerate(grid):
            fragment_data['crank_torsion_drives']['crank_job_1']['torsion_{}'.format(d)] = spacing

    # ToDo define jobs combining different torsions

    return fragment_data


def get_initial_crank_state(fragment, fragment_name=None):
    """
    Generate initial crank state JSON for each crank job in fragment
    Parameters
    ----------
    fragment: dict
        A fragment from JSON crank jobs
    fragment_name: str
        Name for path and file for crank jsonstatefile. Default is None. If None the file does not get written

    Returns
    -------
    crank_initial_states: dict
        dictionary containing JSON specs for initial states for all crank jobs in a fragment.
    """
    crank_initial_states = {}
    init_geometry = fragment['molecule']['geometry']
    needed_torsions = fragment['needed_torsion_drives']
    crank_jobs = fragment['crank_torsion_drives']
    for job in crank_jobs:
        dihedrals = []
        grid_spacing = []
        for torsion in needed_torsions:
            if torsion in crank_jobs[job]:
                # Convert to 0 based numbering which crank uses. Tags start counting at 1
                dihedrals.append([i-1 for i in needed_torsions[torsion]])
                grid_spacing.append(crank_jobs[job][torsion])
        crank_state = {}
        crank_state['dihedrals'] = dihedrals
        crank_state['grid_spacing'] = grid_spacing
        crank_state['elements'] = fragment['molecule']['symbols']

        #ToDo add ability to start with many geomotries
        crank_state['init_coords'] = [init_geometry]
        crank_state['grid_status'] = {}
        if fragment_name:
            # Make directory for job
            current_path = os.getcwd()
            path = os.path.join(current_path, fragment_name + '_{}'.format(job))
            try:
                os.mkdir(path)
            except FileExistsError:
                utils.logger().info('Warning: overwriting {}'.format(path))
            # extend with job number because several jobs can exist in a fragment
            jsonfilename = os.path.join(path, fragment_name + '_{}.json'.format(job))
            outfile = open(jsonfilename, 'w')
            json.dump(crank_state, outfile, indent=2, sort_keys=True, cls=utils.UUIDEncoder)
            outfile.close()

        crank_initial_states[job] = crank_state
    return crank_initial_states


def launch_crank(fragment, inputfile, dihedralfile, init_coords=None, grid=30, engine='psi4', native_opt=False,
                 wq_port=None, verbose=True):
    """
    Launch crank-launch for a specific crank job in fragment.

    Parameters
    ----------
    fragment: dict
        A fragment from JSON crank jobs
    inputfile: str
        Path to psi4 inputfile
    dihedralfile: str
        path to dihedralfile
    init_coords: str
        path to coordinates trajectory if starting the scan from several starting configurations. Default is None and
        the geometry in fragment will be used.
    grid: int or list of ints
        grid spacing for crank torsion scan in degrees. Default is 30. Must be divisor of 360.
    engine: str
        QM engine crank should use. Default psi4. Allowed options:['qchem', 'psi4', 'terachem']
    native_opt: bool
        If True, use QM program native optimizer. Default is False - crank will use geometric for optimizaiton
    wq_port: int
        Specify port number ot use Work Queue to distribute optimization jobs. Default is None - crank will run sequentially.
    verbose: bool
        If True, output will be verbose. Default is True

    """

    #launch crank
    command = 'crank-launch {} {} -g {} -e {}'.format(
            inputfile, dihedralfile, grid, engine)
    if init_coords:
        command = command + ' --init_coords {}'.format(init_coords)
    if native_opt:
        command += ' --native_opt '
    if wq_port:
        command = command + ' --wq_port {} '.format(wq_port)
    if verbose:
        command += ' -v'

    outfile = inputfile.replace('dat', 'out')
    os.system(command + ' > {}'.format(outfile))


def to_crank_input(fragment, mol_name=None, path=None, crank_job='crank_job_1', launch=False, **kwargs):
    """
    Generate crank input files for a fragment (psi4 input file and dihedral file containing torsions that should
    be restrained

    Parameters
    ----------
    fragment: dict
    mol_name: str
        name for molecule. Will be used for filename and molecule name in psi4 input file so it should be a valid Python
        identifier. If None, a name will be generated from SMILES with invalid characters converted to hex. Default is
        None
    path: str
        path to write files to. If None, the directory will be created in current directory. Default None
    crank_job: str
        key to crank job in fragment. Default is crank_job_1
    **kwargs: key arguments for generating psi4 input file

    Returns
    -------
    path: str
        path to where the files were written
    inputfile: str
        absolute path to psi4 input file
    dihedralfile: str
        absolute path to dihedral file

    """
    if not mol_name:
        # Generate Python valid identifier string from smiles by converting forbidden characters to _hex_
        smiles = fragment['canonical_isomeric_SMILES']
        mol_name, namespace = utils.make_python_identifier(smiles, convert='hex')
    if not path:
        cwd = os.getcwd()
        path = os.path.join(cwd, mol_name + '_{}'.format(crank_job))

    # create folder to launch crank in
    try:
        os.mkdir(path)
    except FileExistsError:
        warnings.warn("Overwriting {}".format(path))

    # Write out input files
    inputfile = os.path.join(path, mol_name + '_{}.dat'.format(crank_job))
    utils.to_psi4_input(fragment, molecule_name=mol_name, crank_job=crank_job, filename=inputfile, **kwargs)
    dihedralfile = os.path.join(path, mol_name + '_{}.txt'.format(crank_job))
    utils.to_dihedraltxt(fragment, crank_job=crank_job, filename=dihedralfile)

    return path, inputfile, dihedralfile
