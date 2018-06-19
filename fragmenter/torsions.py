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


def fragment_to_torsion_scan(fragments, json_filename=None, grid=30, terminal_torsion_spacing=30):
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
        j = json.dump(molecules, f, indent=4, sort_keys=True)
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
        # ToDo: save the new tagged SMILES somewhere. Maybe return it?

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

    needed_torsion_scans = {'mid': {}, 'terminal': {}}
    if h_tors:
        h_tors_min = one_torsion_per_rotatable_bond(h_tors)
        for i, tor in enumerate(h_tors_min):
            tor_name = ((tor[0].GetMapIdx()), (tor[1].GetMapIdx()), (tor[2].GetMapIdx()), (tor[3].GetMapIdx()))
            needed_torsion_scans['terminal']['torsion_{}'.format(str(i))] = tor_name
    if mid_tors:
        mid_tors_min = one_torsion_per_rotatable_bond(mid_tors)
        for i, tor in enumerate(mid_tors_min):
            tor_name = ((tor[0].GetMapIdx()), (tor[1].GetMapIdx()), (tor[2].GetMapIdx()), (tor[3].GetMapIdx()))
            needed_torsion_scans['mid']['torsion_{}'.format(str(i))] = tor_name

    # Check that there are no duplicate torsions in mid and h_torsions
    list_tor = list(needed_torsion_scans['mid'].values()) + list(needed_torsion_scans['terminal'].values())
    set_tor = set(list_tor)

    if not len(set_tor) == len(list_tor):
        raise Warning("There is a torsion defined in both mid and terminal torsions. This should not happen. Check "
                      "your molecule and the atom mapping")

    for torsion_type in needed_torsion_scans:
        needed_torsion_scans[torsion_type]['grid_spacing'] = 30

    return needed_torsion_scans


def one_torsion_per_rotatable_bond(torsion_list):
    """
    Keep only one torsion per rotatable bond
    Parameters
    ----------
    torsion_list: list
        list of torsion in molecule

    Returns
    -------
    list of only one torsion per rotatable bonds

    """

    central_bonds = np.zeros((len(torsion_list), 3), dtype=int)
    for i, tor in enumerate(torsion_list):
        central_bonds[i][0] = i
        central_bonds[i][1] = tor[1].GetIdx()
        central_bonds[i][2] = tor[2].GetIdx()

    grouped = central_bonds[central_bonds[:, 2].argsort()]
    sorted_tors = [torsion_list[i] for i in grouped[:, 0]]

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

    return tors


def customize_grid_spacing(fragment_data, mid_grid=15, terminal_grid=None):
    """
    Optional function to create more complicated grids for crank jobs.

    Parameters
    ----------
    fragment_data
    mid_grid: int
        Can be either None, int, list of ints or list of lists of ints. If None, no mid torsion will be driven.
        If int, all mid torsions will be driven at that grid resolution
        If list of ints, the list length must be the same as the amount of mid torsions. Each element in the list corresponds
        to the torsion with that index. Values can also be None and those torsions will not be driven.
        If list of ints, several crank jobs will be generated - each job with grid defined by that list.
    terminal_grid: int, optional, default=None
        The options are the same as mid_grid

    Returns
    -------

    """
    # Check grid specification
    if mid_grid is None and terminal_grid is None:
        warnings.warn("No grid points are specified. Are you sure you do not want to drive any torsions?", Warning)
    mid_torsion_dimension = len(fragment_data['needed_torsion_drives']['mid']) - 1
    terminal_torsion_dimension = len(fragment_data['needed_torsion_drives']['terminal']) - 1

    # Check mid torsion grid specification
    if isinstance(mid_grid, int):
        mid_grid_type = 'int'
    elif isinstance(mid_grid, list):
        # Check for list of list
        mid_grid_type = 'list'
        if any(isinstance(el, list) for el in mid_grid):
            mid_grid_type = 'list_of_lists'
    else:
        mid_grid_type = None
    if mid_grid_type == 'list':
        if len(mid_grid) != mid_torsion_dimension:
            raise Warning("length of grid spacing must be equal to torsion dimension")
    if mid_grid_type == 'list_of_lists':
        # Check length of every list
        if not all(len(el) == mid_torsion_dimension for el in mid_grid):
            raise Warning("dimension of the grid must be equal to dimension of torsion scan")

    # Check terminal torsion grid specification
    if isinstance(terminal_grid, int):
        term_grid_type = 'int'
    if isinstance(terminal_grid, list):
        # Check for list of list
        term_grid_type = 'list'
        if any(isinstance(el, list) for el in terminal_grid):
            term_grid_type = 'list_of_lists'
    else:
        term_grid_type = None
    if term_grid_type == 'list':
        if len(terminal_grid) != terminal_torsion_dimension:
            raise Warning("length of grid spacing must be equal to torsion dimension")
    if term_grid_type == 'list_of_lists':
        # Check length of every list
        if not all(len(el) == terminal_torsion_dimension for el in terminal_grid):
            raise Warning("dimension of the grid must be equal to dimension of torsion scan")

    fragment_data['needed_torsion_drives']['mid']['grid_spacing'] = mid_grid
    fragment_data['needed_torsion_drives']['terminal']['grid_spacing'] = terminal_grid


def define_crank_job(fragment_data, qc_program='Psi4', method='B3LYP', basis='aug-cc-pVDZ', **kwargs):
    """

    Parameters
    ----------
    fragment_data
    grid: int, optional, default 30
        spacing for mid dihedral scan grid points in degree. Must be divisible by 360.
    terminal_torsion_spacing: int, optional, defualt 30
        spacing for terminal (usually trivial) torsions. If None, will not specify crank to drive terminal torsions. Only
        mid torsion will be cranked.
    combinations: list of ints
        can be list of list. Each list defines which torsions in needed_torsion_drives to run on crank.
        Default is None. Only one crank job will be defined for all torsions in needed_torsion_drives
    qc_program: str, optional, default Psi4
    method: str, optional, default B3LYP
    basis: str, optional, default aug-cc-pVDZ
    kwargs

    Returns
    -------
    JSON crank job specs

    """

    needed_torsion_drives = fragment_data['needed_torsion_drives']

    for spacing in (needed_torsion_drives['mid']['grid_spacing'], needed_torsion_drives['terminal']['grid_spacing']):
        if not isinstance(spacing, list):
            spacing = [spacing]
        flattened = utils.flatten(spacing)
        for grid_interval in flattened:
            if grid_interval is None:
                break
            if 360 % grid_interval:
                raise ValueError("grid spacing must be a factor of 360")

    fragment_data['crank_torsion_drives'] = dict()

    model = {'method': method, 'basis': basis}
    options = {'scf_type': 'df'}
    if kwargs:
        options = kwargs['options']

    # Figure out how many crank jobs are needed
    mid_grids = needed_torsion_drives['mid']['grid_spacing']
    terminal_grids = needed_torsion_drives['terminal']['grid_spacing']
    crank_jobs = 1
    if isinstance(mid_grids, list):
        # Check if it is a list of lists
        if any(isinstance(el, list) for el in mid_grids):
            crank_jobs = len(mid_grids)
    if mid_grids is None:
        # No mid torsions are specified. Check terminal torsions
        if terminal_grids is None:
            raise Warning("Are you sure you do not want to specify any crank jobs for fragment {}".format(
                    fragment_data['canonical_isomeric_SMILES']))
        elif isinstance(terminal_grids, int):
            pass
        elif any(isinstance(el, list) for el in terminal_grids):
            crank_jobs = len(terminal_grids)

    if crank_jobs == 1:
        fragment_data['crank_torsion_drives'] = {'crank_job_0': {'crank_specs': {'model': model, 'options': options},
                                                 'mid_torsions': {}, 'terminal_torsions': {}}}
        if isinstance(mid_grids, int):
            # All torsions are driven at the same resolution
            torsions_to_drive = list(needed_torsion_drives['mid'].keys())
            torsions_to_drive.remove('grid_spacing')
            fragment_data['crank_torsion_drives']['crank_job_0']['mid_torsions'] = {torsion: mid_grids for torsion in torsions_to_drive}
        if isinstance(mid_grids, list):
            # Only constrain torsions that have a corresponding interval
            for i, interval in enumerate(mid_grids):
                if interval is not None:
                    fragment_data['crank_torsion_drives']['crank_job_0']['mid_torsions']['torsion_{}'.format(str(i))] = interval
        if isinstance(terminal_grids, int):
            torsions_to_drive = list(needed_torsion_drives['terminal'].keys())
            torsions_to_drive.remove('grid_spacing')
            fragment_data['crank_torsion_drives']['crank_job_0']['terminal_torsions'] = {torsion: terminal_grids for torsion in torsions_to_drive}
        if isinstance(terminal_grids, list):
            for i, interval in enumerate(terminal_grids):
                if interval is not None:
                    fragment_data['crank_torsion_drives']['crank_job_0']['terminal_torsions']['torsion_{}'.format(str(i))] = interval

    else:

        for i in range(crank_jobs):
            fragment_data['crank_torsion_drives']['crank_job_{}'.format(str(i))] = {'crank_specs': {
                'model': model, 'options': options}, 'mid_torsions': {}, 'terminal_torsions': {}}
            mid_grid = mid_grids[i]
            for e, interval in enumerate(mid_grid):
                if interval is not None:
                    fragment_data['crank_torsion_drives']['crank_job_{}'.format(str(i))]['mid_torsions']['torsion_{}'.format(str(e))] = interval
            # Check terminal
            if isinstance(terminal_grids, int):
                # Terminal torsions are all driven at the same resolution
                torsions_to_drive = list(needed_torsion_drives['terminal'].keys())
                torsions_to_drive.remove('grid_spacing')
                print(torsions_to_drive)
                fragment_data['crank_torsion_drives']['crank_job_{}'.format(str(i))]['terminal_torsions'] = {torsion: terminal_grids for torsion in torsions_to_drive}
            if isinstance(terminal_grids, list):
                terminal_grid = terminal_grids[i]
                for e, interval in enumerate(terminal_grid):
                    if interval is not None:
                        fragment_data['crank_torsion_drives']['crank_job_{}'.format(str(i))]['terminal_torsions']['torsion_{}'.format(str(e))] = interval

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
    for i, job in enumerate(crank_jobs):
        dihedrals = []
        grid_spacing = []
        needed_mid_torsions = needed_torsions['mid']
        for mid_torsion in crank_jobs[job]['mid_torsions']:
            dihedrals.append([j-1 for j in needed_mid_torsions[mid_torsion]])
            grid_spacing.append(crank_jobs[job]['mid_torsions'][mid_torsion])
        needed_terminal_torsions = needed_torsions['terminal']
        for terminal_torsion in crank_jobs[job]['terminal_torsions']:
            dihedrals.append([j-1 for j in needed_terminal_torsions[terminal_torsion]])
            grid_spacing.append(crank_jobs[job]['terminal_torsions'][terminal_torsion])

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
