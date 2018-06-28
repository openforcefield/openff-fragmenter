import socket
import uuid
import getpass
import yaml

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
_default_options['torsion_scan'] = {'grid':30,
                                     'terminal_torsion_spacing':30}


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


def launch_fragmenter(molecule, options=None):
    """
    Launch fragmenter

    Parameters
    ----------
    molecule
    options

    Returns
    -------

    """

    provenance = get_provenance()

    if options is not None:
        user_options = load_options(options)
        expand_state_options = user_options['expand_states']
        fragment_options = user_options['generate_fragments']
        torsions_options = user_options['torsion_scan']
    else:
        # Use default options
        expand_state_options = _default_options['expand_states']
        fragment_options = _default_options['generate_fragments']
        torsions_options = _default_options['torsion_scan']

    provenance['routine_options'] = {'expand_states': expand_state_options,
                                     'generate_fragments': fragment_options,
                                     'torsion_scan': torsions_options}

    molecule = check_molecule(molecule)

    utils.logger().info('fragmenting {}...'.format(molecule.GetTitle()))

    if expand_state_options:
        states = fragmenter.expand_states(molecule, **expand_state_options)
        fragments = fragmenter.generate_fragments(states, **fragment_options)
    else:
        fragments = fragmenter.generate_fragments(molecule, **fragment_options)

    return fragments
    #
    # # Find torsions
    # for parent in fragments:
    #     for frag in fragments[parent]:
    #         json_specs = dict()
    #         json_specs['provenance'] = provenance
    #         json_specs['provenance']['parent_molecule'] = parent
    #         json_specs['canonical_isomeric_SMILES'] = frag
    #         molecule = utils.smiles_to_oemol(frag)
    #         json_specs['canonical_SMILES'] = oechem.OECreateCanSmiString(molecule)
    #         explicit_h_isomeric = utils.create_mapped_smiles(molecule, tagged=False)
    #         json_specs['explicit_hydrogen_canonical_isomeric_SMILES'] = explicit_h_isomeric
    #         explicit_h = utils.create_mapped_smiles(molecule, tagged=False, isomeric=False)
    #         json_specs['explicit_hydrogen_canonical_SMILES'] = explicit_h
    #         tagged_SMARTS = utils.create_mapped_smiles(molecule)
    #         json_specs['tagged_SMARTS'] = tagged_SMARTS
    #         molecule, atom_map = utils.get_atom_map(tagged_SMARTS, is_mapped=True)
    #         # Find formal charge
    #         charge = 0
    #         for atom in molecule.GetAtoms():
    #             charge += atom.GetFormalCharge()
    #         QC_JSON_molecule = utils.to_mapped_QC_JSON_geometry(molecule, atom_map, charge=charge)
    #         json_specs['molecule'] = QC_JSON_molecule
    #         needed_torsion_drives = fragmenter.find_torsions(molecule)
    #         json_specs['needed_torsion_drives'] = needed_torsion_drives
    #         define_crank_job(json_specs)
    #         molecules[frag] = json_specs
    #
    # if json_filename:
    #     f = open(json_filename, 'w')
    #     j = json.dump(molecules, f, indent=4, sort_keys=True)
    #     f.close()

    #return molecules






