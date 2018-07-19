import socket
import uuid
import getpass
import yaml
import json
import os
import warnings

import fragmenter
from fragmenter import fragment, torsions, utils
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
_default_options['enumerate_states'] = {'protonation': True,
                                     "tautomers":False,
                                     'stereoisomers':True,
                                     'max_states': 200,
                                     'level': 0,
                                     'reasonable': True,
                                     'carbon_hybridization': True,
                                     'suppress_hydrogen': True}
_default_options['enumerate_fragments'] = {'strict_stereo': True,
                                          'combinatorial': True,
                                          'MAX_ROTORS': 2,
                                          'remove_map': True}
_default_options['generate_crank_jobs'] = {'max_conf': 1,
                                    'terminal_torsion_resolution': 30,
                                    'internal_torsion_resolution': 30,
                                    'scan_internal_external_combination': 0,
                                    'scan_dimension': 2,
                                    'options':{
                                    'qc_program': 'psi4',
                                    'method': 'B3LYP',
                                    'basis': 'aug-cc-pVDZ'}}


def enumerate_states(molecule, options=None, json_filename=None):
    """
    enumerate protonation, tautomers and stereoisomers for molecule. Default does not enumerate tautomers, only
    protonation/ionization states and stereoisomers

    Parameters
    ----------
    molecule: any format that OpenEye pareses. Can be path to file containing molecule or SMILES/Inchi string
    options: str, optional. Default None
        path to user options yaml file. If None will use default options
    json_filename: str, optional. default None
        path to json output file. If None, no file will be written

    Returns
    -------
    json_dict: dict
        dictoinary containing canonicla isomeric SMILES for states and provenane.

    """
    # Load options for enumerate states
    provenance = _get_provenance('enumerate_states', options=options)
    options = provenance['routine']['enumerate_states']['keywords']

    molecule = utils.check_molecule(molecule)
    can_iso_smiles = oechem.OEMolToSmiles(molecule)
    states = fragment.expand_states(molecule, **options)

    # Remove extraneous options
    routine_options = _remove_extraneous_options(user_options=options, routine='enumerate_states')

    provenance['routine']['enumerate_states']['keywords'] = routine_options
    provenance['routine']['enumerate_states']['parent_molecule'] = can_iso_smiles
    json_dict = {'provenance': provenance, 'states': states}

    if json_filename:
        json_dict['states'] = list(json_dict['states'])
        with open(json_filename, 'w') as f:
            json.dump(json_dict, f, indent=2, sort_keys=True)

    return json_dict


def enumerate_fragments(molecule, mol_provenance=None, options=None, json_filename=None):
    """
    Fragment molecule

    Parameters
    ----------
    molecule: str
        SMILES string of molecule to fragment
    mol_provenance: dict, optional. Default is None
        provenance for molecule. If the molecule is a state from enumerate_states, the provenance from enumerate_states
        should be used
    options: str, optional, Default None
        path to yaml file with user options. If None will use default options
    json_filename: str, optional. Default None
        If a filename is provided, will write output to json file.

    Returns
    -------
    json_dict: dict
        dictionary containing provenance and fragments.

    """
    provenance = _get_provenance(routine='enumerate_fragments', options=options)
    provenance['routine']['enumerate_fragments']['parent_molecule'] = molecule
    options = provenance['routine']['enumerate_fragments']['keywords']

    parent_molecule = utils.check_molecule(molecule)
    fragments = fragment.generate_fragments(parent_molecule, **options)

    routine_options = _remove_extraneous_options(user_options=options, routine='enumerate_fragments')

    if mol_provenance:
        provenance['routine']['enumerate_states'] = mol_provenance['routine']['enumerate_states']
    provenance['routine']['enumerate_fragments']['keywords'] = routine_options

    # Generate SMILES
    fragments_json_dict = {}
    for fragm in fragments:
        for frag in fragments[fragm]:
            SMILES = {}
            fragment_mol = utils.smiles_to_oemol(frag)
            SMILES['canonical_SMILES'] = utils.create_mapped_smiles(fragment_mol, tagged=False, isomeric=False,
                                                                    explicit_hydrogen=False)
            SMILES['canonical_isomeric_SMILES'] = utils.create_mapped_smiles(fragment_mol, tagged=False, explicit_hydrogen=False)
            SMILES['canonical_explicit_hydrogen_SMILES'] = utils.create_mapped_smiles(fragment_mol, tagged=False, isomeric=False)
            SMILES['canonical_isomeric_explicit_hydrogen_SMILES'] = utils.create_mapped_smiles(fragment_mol, tagged=False)
            SMILES['canonical_isomeric_explicit_hydrogen_mapped_SMILES'] = utils.create_mapped_smiles(fragment_mol)

            fragments_json_dict[frag] = {'SMILES': SMILES}

            # Generate QM molecule
            mol, atom_map = utils.get_atom_map(tagged_smiles=SMILES['canonical_isomeric_explicit_hydrogen_mapped_SMILES'],
                                               molecule=fragment_mol, is_mapped=True)

            charge = utils.get_charge(mol)
            qm_mol = utils.to_mapped_QC_JSON_geometry(mol, atom_map, charge=charge)
            fragments_json_dict[frag]['molecule'] = qm_mol
            fragments_json_dict[frag]['provenance'] = provenance


    if json_filename:
        with open(json_filename, 'w') as f:
            json.dump(fragments_json_dict, f, indent=2, sort_keys=True)

    return fragments_json_dict

def generate_crank_jobs(fragment_dict, options=None, fragment_name=None):
    """
    Generate crank initial state from fragment json specs

    Parameters
    ----------
    fragment_dict: dict
        JSON spec for fragment
    options: str, Optional. Default None
        path to yaml file with options. If None will use default options
    fragment_name: str, optional. Default None
        The base name for crank job names. If name is given, each crank job will be written to its own directory with the
        initial state json file
    Returns
    -------
    crank_initial_states: dict
        dictionary of crank jobs

    """

    provenance = _get_provenance(routine='generate_crank_jobs', options=options)
    options = provenance['routine']['generate_crank_jobs']['keywords']

    mapped_smiles = fragment_dict['SMILES']['canonical_isomeric_explicit_hydrogen_mapped_SMILES']
    OEMol = utils.smiles_to_oemol(mapped_smiles)

    # Check if molecule is mapped
    if not utils.is_mapped(OEMol):
        warnings.warn("OEMol is not mapped. Creating a new mapped SMILES")
        fragment_dict['SMILES']['canonical_isomeric_explicit_hydrogen_mapped_SMILES'] = utils.create_mapped_smiles(OEMol)

    needed_torsion_drives = torsions.find_torsions(OEMol)
    fragment_dict['needed_torsion_drives'] = needed_torsion_drives

    crank_job_specs = torsions.define_crank_job(fragment_dict, **options)
    crank_initial_states = torsions.get_initial_crank_state(crank_job_specs)

    routine_options = _remove_extraneous_options(user_options=options, routine='generate_crank_jobs')

    for crank_job in crank_initial_states:
        crank_initial_states[crank_job]['provenance'] = fragment_dict['provenance']
        crank_initial_states[crank_job]['provenance']['SMILES'] = fragment_dict['SMILES']
        crank_initial_states[crank_job]['provenance']['routine']['generate_crank_jobs'] = provenance['routine']['generate_crank_jobs']
        crank_initial_states[crank_job]['provenance']['routine']['generate_crank_jobs']['keywords'] = routine_options

    return crank_initial_states


def workflow(molecules_smiles, options=None, write_json_intermediate=False, write_json_crank_job=True):
    """
    Convenience function to run Fragmenter workflow.

    Parameters
    ----------
    molecules_smiles: list of str
        list of SMILES
    options: str, optional. Default is None
        path to yaml file with options.
    write_json_intermediate: bool, optional. Default False
        If True will write JSON files for intermediate steps (enumerating states and fragments)
    write_json_crank_job: bool, optional. Default Tre
        If True, will create a directory for crank job with crank initial state. The name of the directory and file will
        be the canonical SMILES of the fragment.

    Returns
    -------
    molecules: dict
        JSON specs for crank jobs. Keys are canonical SMILES of fragment.

    """

    # Check input
    if not isinstance(molecules_smiles, list):
        molecules_smiles = [molecules_smiles]

    all_frags = {}
    all_crank_jobs = {}
    for i, molecule in enumerate(molecules_smiles):
        json_filename = None
        if write_json_intermediate:
            json_filename = molecule  + '_states.json'
        states = enumerate_states(molecule, options=options, json_filename=json_filename)
        for state in states['states']:
            if write_json_intermediate:
                json_filename = state  + '_fragments.json'
            fragments = enumerate_fragments(state, mol_provenance=states['provenance'], options=options,
                                            json_filename=json_filename)
            all_frags.update(**fragments)
    for frag in all_frags:
        crank_jobs = generate_crank_jobs(all_frags[frag], options=options)
        all_crank_jobs[frag] = crank_jobs
        if write_json_crank_job:
            for job in crank_jobs:
                # Make directory for job
                current_path = os.getcwd()
                path = os.path.join(current_path, frag + '_{}'.format(job))
                try:
                    os.mkdir(path)
                except FileExistsError:
                    utils.logger().info('Warning: overwriting {}'.format(path))
                # extend with job number because several jobs can exist in a fragment
                jsonfilename = os.path.join(path, frag + '_{}.json'.format(job))
                outfile = open(jsonfilename, 'w')
                json.dump(crank_jobs[job], outfile, indent=2, sort_keys=True)
                outfile.close()

    return all_crank_jobs


def _get_provenance(routine, options=None):
    """
    Get provenance with keywords for routine

    Parameters
    ----------
    routine: str
        routine to get provenance for. Options are 'enumerate_states', 'enumerate_fragments', and 'generate_crank_jobs'
    options: str, optional. Default is None
        path to yaml file containing user specified options.

    Returns
    -------
    provenance: dict
        dictionary with provenance and routine keywords.

    """
    provenance = {'creator': fragmenter.__package__,
                  'routine': {
                      routine: {
                          'version': fragmenter.__version__,
                          'keywords': {},
                          'job_id': str(uuid.uuid4()),
                          'hostname': socket.gethostname()}
                  },
                  'username': getpass.getuser()}

    if routine == 'enumerate_states':
        package = _canonicalization_details['package']
        can_is_smiles = _canonicalization_details['canonical_isomeric_SMILES']
        canonicalization_details = {'package': package,
                                    'canonical_isomeric_SMILES': can_is_smiles}
        provenance['canonicalization_details'] = canonicalization_details
    if routine == 'enumerate_fragments':
        canonicalization_details = _canonicalization_details
        provenance['canonicalization_details'] = canonicalization_details

    provenance['routine'][routine]['keywords'] = _load_options(routine, options)

    return provenance


def _load_options(routine, load_path=None):
    """
    load user options

    Parameters
    ----------
    routine: str
        routine to load options for. Options are 'enumerate_states', 'enumerate_fragments', and 'generate_crank_jobs'
    load_path: str, optional. Default None
        path to user option yaml file. If none will use default options

    Returns
    -------
    options: dict
        keywords for routine

    """

    if load_path:
        with open(load_path) as stream:
            user_config = yaml.load(stream)

        options = {**_default_options[routine], **user_config[routine]}
    else:
        options = _default_options[routine]
    return options

def _remove_extraneous_options(user_options, routine):

    options = {option: user_options[option] for option in _default_options[routine]}

    return options

def combine_json_fragments(json_inputs, json_output=None):
    """
    This function takes a list of json input fragment files and returns a dictionary with redundant fragments removed.

    Parameters
    ----------
    json_inputs: list of json fragments input files
    json_output: str, optional, Default None
        If not none, the new json fragment dictionary will be written to json file

    Returns
    -------
    fragments: dict
        A dictionary containing all fragments without redundant fragments removed.
    """
    if not isinstance(json_inputs, list):
        json_inputs = [json_inputs]

    fragments = {}
    for json_file in json_inputs:
        with open(json_file, 'r') as f:
            fragments.update(json.load(f))

    if json_output:
        with open(json_output, 'r') as f:
            json.dump(fragments, f, indent=2, sort_keys=True)

    return fragments


