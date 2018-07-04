import socket
import uuid
import getpass
import yaml
import json
import os
import warnings

import fragmenter
from fragmenter import fragment, torsions, chemi, utils
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


def get_provenance(routine, options=None):
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

    provenance['routine'][routine]['keywords'] = load_options(routine, options)

    return provenance


def load_options(routine, load_path=None):
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

def remove_extraneous_options(user_options, routine):

    options = {option: user_options[option] for option in _default_options[routine]}

    return options

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
    provenance = get_provenance('enumerate_states', options=options)
    options = provenance['routine']['enumerate_states']['keywords']

    molecule = utils.check_molecule(molecule)
    can_iso_smiles = oechem.OEMolToSmiles(molecule)
    states = fragment.expand_states(molecule, **options)

    # Remove extraneous options
    routine_options = remove_extraneous_options(user_options=options, routine='enumerate_states')

    provenance['routine']['enumerate_states']['keywords'] = routine_options
    json_dict = {'provenance': provenance,
                  can_iso_smiles: states}
    if json_filename:
        with open(json_filename, 'w') as f:
            json.dump(json_dict, f, indent=2, sort_keys=True)

    return json_dict


def enumerate_fragments(molecule, options=None, json_filename=None):
    pass

def generate_crank_jobs(molecule, options=None, json_filename=None):
    pass


def workflow(molecule, options=None, json_filename=None):
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




