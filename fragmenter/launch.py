import socket
import uuid
import getpass
import yaml
import json
import os
import warnings

import fragmenter
from fragmenter import utils
import openeye as oe
from openeye import oechem


_OPENEYE_VERSION = oe.__name__ + '-v' + oe.__version__
_canonicalization_details = {'package': _OPENEYE_VERSION,
                             'canonical_isomeric_SMILES': {'Flags': ['ISOMERIC', 'Isotopes', 'AtomStereo',
                                                                    'BondStereo', 'Canonical', 'AtomMaps', 'RGroups'],
                                                           'oe_function': 'openeye.oechem.OEMolToSmiles(molecule)'},
                             'canonical_SMILES': {'Flags': ['DEFAULT', 'Canonical', 'AtomMaps', 'RGroups'],
                                                  'oe_function': 'openeye.oechem.OECreateCanSmiString(molecule)'},
                             'canonical_isomeric_explicit_hydrogen_SMILES': {'Flags': ['Hydrogens', 'Isotopes', 'AtomStereo',
                                                                                      'BondStereo', 'Canonical', 'RGroups'],
                                                                            'oe_function': 'openeye.oechem.OECreateSmiString()'},
                             'canonical_explicit_hydrogen_SMILES': {'Flags': ['Hydrogens', 'Canonical', 'RGroups'],
                                                                   'oe_function': 'openeye.oechem.OECreateSmiString()'},
                             'notes': 'All other available OESMIELSFlag are set to False'}

_default_options = {}
_default_options['expand_states'] = {'protonation': True,
                                     "tautomers":False,
                                     'stereoisomers':True,
                                     'max_states': 200,
                                     'level': 0,
                                     'reasonable': True,
                                     'carbon_hybridization': True,
                                     'suppress_hydrogen': True}
_default_options['generate_fragments'] = {'strict_stereo': True,
                                          'combinatorial': True,
                                          'MAX_ROTORS': 2,
                                          'remove_map': True}
_default_options['torsion_scan'] = {'mid_grid':30,
                                    'terminal_grid':30,
                                    '1D_scans': False,
                                    'options':{
                                    'qc_program': 'psi4',
                                    'method': 'B3LYP',
                                    'basis': 'aug-cc-pVDZ'}}


def get_provenance():
    provenance = {'job_id': str(uuid.uuid4()),
                  'creator': fragmenter.__package__,
                  'version': fragmenter.__version__,
                  'canonicalization_details': _canonicalization_details,
                  'hostname': socket.gethostname(),
                  'username': getpass.getuser()}

    return provenance


def load_options(load_path):

    options = {}
    options['config_path'] = load_path

    with open(load_path) as stream:
        user_config = yaml.load(stream)

        # override default
        options = {}
        for option in _default_options:
            try:
                if user_config[option] is False:
                    options[option] = False
                    continue
            except KeyError:
                options[option] = _default_options[option]
            try:
                options[option] = {**_default_options[option], **user_config[option]}
            except KeyError:
                options[option] = _default_options[option]

    return options

def specify_1D_grid(molecule_data):
    mid = molecule_data['needed_torsion_drives']['mid']
    mid_grid = None

    terminal = molecule_data['needed_torsion_drives']['terminal']
    terminal_grid = None

    if not mid['grid_spacing'] and not terminal['grid_spacing']:
        warnings.warn("No grid points are specified. Are you sure you do not want to drive any torsions?", Warning)

    mid_dimension = len(mid) -1
    terminal_dimesion = len(terminal) -1
    dimension = mid_dimension + terminal_dimesion

    if mid['grid_spacing'] is not None:
        mid_grid = [[None]*dimension for i in range(mid_dimension)]
        for i, j in enumerate(mid_grid):
            j[i] = mid['grid_spacing']
    if terminal['grid_spacing'] is not None:
        terminal_grid = [[None]*dimension for i in range(terminal_dimesion)]
        for i, j in enumerate(terminal_grid):
            j[i+mid_dimension] = terminal['grid_spacing']

    return mid_grid, terminal_grid

def check_molecule(molecule):

    mol = oechem.OEMol()
    # First try reading as smiles
    if not oechem.OESmilesToMol(mol, molecule):
        # Try reading as input file
        ifs = oechem.oemolistream()
        if not ifs.open(molecule):
            raise Warning('Could not parse molecule.')

    # normalize molecule
    title = mol.GetTitle()
    molecule = utils.normalize_molecule(mol, title=title)
    return molecule


def launch_fragmenter(molecule, options=None, json_filename=None):
    """
    Launch fragmenter

    Parameters
    ----------
    molecule: molecule to fragment.
        Can be in any format that OpenEye accepts (file or string)
    options: str, optional. Default is None
        path to yaml file with options.
    json_filename: str, optional, Default is None
        name of json output file. If given will write crank JSON specs to file

    Returns
    -------
    molecules: dict
        JSON specs for crank jobs. Keys are canonical SMILES of fragment.

    """

    provenance = get_provenance()

    if options is not None:
        user_options = load_options(options)
        expand_state_options = user_options['expand_states']
        fragment_options = user_options['generate_fragments']
        scan_options = user_options['torsion_scan']
        if not expand_state_options:
            routine_options = {}
            routine_options['expand_states'] = False
            user_options.pop('expand_states')
        # remove options not needed for reproducibility from provenance
        routine_options = {routine: {option: user_options[routine][option] for option in _default_options[routine]} for routine in user_options}
    else:
        # Use default options
        expand_state_options = _default_options['expand_states']
        fragment_options = _default_options['generate_fragments']
        scan_options = _default_options['torsion_scan']

        routine_options = _default_options

    provenance['routine_options'] = routine_options

    molecule = check_molecule(molecule)

    utils.logger().info('fragmenting {}...'.format(molecule.GetTitle()))

    if expand_state_options:
        states = fragmenter.expand_states(molecule, **expand_state_options)
        fragments = fragmenter.generate_fragments(states, **fragment_options)
    else:
        fragments = fragmenter.generate_fragments(molecule, **fragment_options)

    # Find torsions
    molecules = {}
    for parent in fragments:
        for frag in fragments[parent]:
            json_specs = dict()
            json_specs['provenance'] = provenance
            json_specs['provenance']['parent_molecule'] = parent
            json_specs['canonical_isomeric_SMILES'] = frag
            molecule = utils.smiles_to_oemol(frag)
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
            needed_torsion_drives = fragmenter.find_torsions(molecule)
            json_specs['needed_torsion_drives'] = needed_torsion_drives

            if scan_options['mid_grid'] != 30 or scan_options['terminal_grid'] != 30:
                # Customize scan grid
                fragmenter.customize_grid_spacing(json_specs, scan_options['mid_grid'], scan_options['terminal_grid'])
            if scan_options['1D_scans']:
                print('1D scan')
                # create 1D torsion scans
                mid_grid, term_grid = specify_1D_grid(json_specs)
                json_specs['needed_torsion_drives']['mid']['grid_spacing'] = mid_grid
                json_specs['needed_torsion_drives']['terminal']['grid_spacing'] = term_grid

            fragmenter.define_crank_job(json_specs, **scan_options['options'])
            molecules[frag] = json_specs

    if json_filename:
        f = open(json_filename, 'w')
        j = json.dump(molecules, f, indent=2, sort_keys=True)
        f.close()

    return molecules

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
            json.dump(crank_state, outfile, indent=2, sort_keys=True)
            outfile.close()

        crank_initial_states[job] = crank_state
    return crank_initial_states


