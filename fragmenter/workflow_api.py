import socket
import uuid
import getpass
import yaml
import json
import os
from pkg_resources import resource_filename
import fragmenter
from fragmenter import fragment, torsions, utils, chemi
from openeye import oechem
from cmiles import to_molecule_id


def enumerate_states(molecule, workflow_id, options=None, title='', json_filename=None):
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
        dictoinary containing canonical isomeric SMILES for states and provenane.

    """
    # Load options for enumerate states
    routine = 'enumerate_states'
    provenance = _get_provenance(workflow_id=workflow_id, routine=routine)
    if not options:
        options = _get_options(workflow_id, routine)
    #options = provenance['routine']['enumerate_states']['keywords']

    molecule = chemi.standardize_molecule(molecule, title=title)
    can_iso_smiles = oechem.OEMolToSmiles(molecule)
    states = fragment.expand_states(molecule, **options)

    # Remove extraneous options
    #routine_options = _remove_extraneous_options(user_options=options, routine='enumerate_states')

    #provenance['routine']['enumerate_states']['keywords'] = routine_options
    provenance['routine']['enumerate_states']['parent_molecule'] = can_iso_smiles
    provenance['routine']['enumerate_states']['parent_molecule_name'] = molecule.GetTitle()
    json_dict = {'provenance': provenance, 'states': states}

    if json_filename:
        json_dict['states'] = list(json_dict['states'])
        with open(json_filename, 'w') as f:
            json.dump(json_dict, f, indent=2, sort_keys=True)

    return json_dict


def enumerate_fragments(molecule, workflow_id, options=None, mol_provenance=None, title='', json_filename=None,
                        generate_vis=False):
    """
    Fragment molecule

    Parameters
    ----------
    molecule: Input molecule. Very permissive. Can be anything that OpenEye can parse
        SMILES string of molecule to fragment
    title: str, optional. Default empty str
        The title or name of the molecule. If empty stirng will use the IUPAC name for molecule title.
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
    routine = 'enumerate_fragments'
    provenance = _get_provenance(workflow_id=workflow_id, routine=routine)
    if options is None:
        options = _get_options(workflow_id, routine)

    parent_molecule = chemi.standardize_molecule(molecule, title)
    provenance['routine']['enumerate_fragments']['parent_molecule_name'] = parent_molecule.GetTitle()
    provenance['routine']['enumerate_fragments']['parent_molecule'] = oechem.OEMolToSmiles(parent_molecule)

    fragments = fragment.generate_fragments(parent_molecule, generate_vis, **options)

    if mol_provenance:
        provenance['routine']['enumerate_states'] = mol_provenance['routine']['enumerate_states']

    # Generate identifiers for fragments
    fragments_json_dict = {}
    for fragm in fragments:
        for i, frag in enumerate(fragments[fragm]):
            fragment_mol = chemi.smiles_to_oemol(frag)
            identifiers = to_molecule_id(fragment_mol, canonicalization='openeye')

            frag = identifiers['canonical_isomeric_smiles']
            fragments_json_dict[frag] = {'identifiers': identifiers}
            fragments_json_dict[frag]['provenance'] = provenance
            fragments_json_dict[frag]['provenance']['canonicalization'] = identifiers.pop('provenance')

    if json_filename:
        with open(json_filename, 'w') as f:
            json.dump(fragments_json_dict, f, indent=2, sort_keys=True)

    return fragments_json_dict


def generate_torsiondrive_input(fragment_dict, workflow_id, options=None, json_filename=None):
    """

    Parameters
    ----------
    fragment_dict
    workflow

    Returns
    -------

    """

    options = _get_options(workflow_id=workflow_id, routine='torsiondriver_input')
    provenance = _get_provenance(workflow_id=workflow_id, routine='torsiondriver_input')
    fragment_dict['provenance']['routine']['torsiondriver_input'] = provenance['routine']['torsiondriver_input']
    provenance = fragment_dict['provenance']

    mol_id = fragment_dict['identifiers']

    if not options:
        options = _get_options(workflow_id, 'torsiondrive_input')

    mapped_smiles = mol_id['canonical_isomeric_explicit_hydrogen_mapped_smiles']
    mapped_mol = chemi.smiles_to_oemol(mapped_smiles)

    if options.pop('max_conf') == 1:
        # Generate 1 conformation for all jobs
        try:
            qm_mol = chemi.to_mapped_QC_JSON_geometry(mapped_mol)
        except RuntimeError:
            utils.logger().warning("{} does not have coordinates. This can happen for several reasons related to Omega. "
                          "{} will not be included in fragments dictionary".format(
                    mol_id['canonical_isomeric_smiles'], mol_id['canonical_isomeric_smiles']))

    torsiondriver_inputs = []
    needed_torsion_drives = torsions.define_torsiondrive_jobs(torsions.find_torsions(mapped_mol), **options)
    for job in needed_torsion_drives:
        print(needed_torsion_drives[job])
        torsiondriver_input = dict()
        torsiondriver_input['molecule'] = qm_mol
        torsiondriver_input['molecule']['identifiers'] = mol_id
        torsiondriver_input['dihedrals'] = needed_torsion_drives[job]['dihedrals']
        torsiondriver_input['grid_spacing'] = needed_torsion_drives[job]['grid_spacing']
        torsiondriver_input['provenance'] = provenance
        torsiondriver_inputs.append(torsiondriver_input)

    if json_filename:
        with open(json_filename, 'w') as f:
            json.dump(torsiondriver_inputs, f, indent=2, sort_keys=True)

    return torsiondriver_inputs


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

    mapped_smiles = fragment_dict['SMILES']['canonical_isomeric_explicit_hydrogen_mapped_smiles']
    OEMol = chemi.smiles_to_oemol(mapped_smiles)

    # Check if molecule is mapped
    if not chemi.is_mapped(OEMol):
        utils.logger().warning("OEMol is not mapped. Creating a new mapped SMILES")
        fragment_dict['SMILES']['canonical_isomeric_explicit_hydrogen_mapped_smiles'] = chemi.create_mapped_smiles(OEMol)

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


def workflow(molecules_smiles, molecule_titles=None, options=None, write_json_intermediate=False, write_json_crank_job=True):
    """
    Convenience function to run Fragmenter workflow.

    Parameters
    ----------
    molecules_smiles: list of str
        list of SMILES
    molecule_titles: list of molecule names, optional. Default None
        If list of molecule names is provided, the name will be added to provenance and used for filenames
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
    for i, molecule_smile in enumerate(molecules_smiles):
        json_filename = None
        title = ''
        if molecule_titles:
            title = molecule_titles[i]
        if write_json_intermediate and title:
            json_filename = 'states_{}.json'.format(title)
        if write_json_intermediate and not title:
            json_filename = 'states_{}.json'.format(utils.make_python_identifier(molecule_smile)[0])
        states = enumerate_states(molecule_smile, title=title, options=options, json_filename=json_filename)
        for j, state in enumerate(states['states']):
            if write_json_intermediate and title:
                json_filename =  'fragments_{}_{}.json'.format(title, j)
            if write_json_intermediate and not title:
                json_filename = 'fragment_{}_{}.json'.format(utils.make_python_identifier(state)[0], j)
            fragments = enumerate_fragments(state, title=title, mol_provenance=states['provenance'], options=options,
                                            json_filename=json_filename)
            all_frags.update(**fragments)
    namespaces = {}
    for frag in all_frags:
        crank_jobs = generate_crank_jobs(all_frags[frag], options=options)
        all_crank_jobs[frag] = crank_jobs

        if write_json_crank_job:
            name, namespace = utils.make_python_identifier(frag, namespace=namespaces, convert='hex', handle='force')
            namespaces.update(namespace)
            for job in crank_jobs:
                # Make directory for job
                current_path = os.getcwd()
                path = os.path.join(current_path, name + '_{}'.format(job))
                try:
                    os.mkdir(path)
                except FileExistsError:
                    utils.logger().warning('Warning: overwriting {}'.format(path))
                # extend with job number because several jobs can exist in a fragment
                jsonfilename = os.path.join(path, name + '_{}.json'.format(job))
                outfile = open(jsonfilename, 'w')
                json.dump(crank_jobs[job], outfile, indent=2, sort_keys=True)
                outfile.close()

    return all_crank_jobs


def _get_provenance(workflow_id, routine):
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
    # ToDo Figure out what to do about routine specific provenances (parent molecules)
    provenance = {'creator': fragmenter.__package__,
                  'job_id': str(uuid.uuid4()),
                  'hostname': socket.gethostname(),
                  'username': getpass.getuser(),
                  'workflow_id': workflow_id,
                  'routine': {routine: {
                      'version': fragmenter.__version__
                  }}}

    # if routine == 'enumerate_fragments':
    #     #canonicalization_details = _canonicalization_details
    #     provenance['canonicalization_details'] = canonicalization_details

    #provenance['routine'][routine]['keywords'] = _load_options(routine, options)

    return provenance


def _get_options(workflow_id, routine):

    fn = resource_filename('fragmenter', os.path.join('data', 'workflows.json'))
    f = open(fn)
    options = json.load(f)[workflow_id]['fragmenter'][routine]['options']
    f.close()

    return options
# def _load_options(routine, load_path=None):
#     """
#     load user options
# 
#     Parameters
#     ----------
#     routine: str
#         routine to load options for. Options are 'enumerate_states', 'enumerate_fragments', and 'generate_crank_jobs'
#     load_path: str, optional. Default None
#         path to user option yaml file. If none will use default options
# 
#     Returns
#     -------
#     options: dict
#         keywords for routine
# 
#     """
# 
#     if load_path:
#         with open(load_path) as stream:
#             user_config = yaml.load(stream)
#             # Check options
#             functions = list(user_config.keys())
#             allowed_functions = list(_default_options.keys())
#             for f in functions:
#                 if f not in allowed_functions:
#                     raise KeyError("{} is not a function. Only function names allowed are {}".format(f, allowed_functions))
#                 user_options = user_config[f].keys()
#                 allowed_options = _allowed_options[f]
#                 for option in user_options:
#                     if option not in allowed_options:
#                         if f == 'generate_crank_jobs':
#                             utils.logger().warning("Is {} an allowed keyword for {}? Please double check".format(option, f))
#                             continue
#                         raise KeyError("{} is not an allowed option for {}. Allowed options are {}".format(option, f, allowed_options))
# 
#         options = {**_default_options[routine], **user_config[routine]}
#     else:
#         options = _default_options[routine]
#     return options
# 
# def _remove_extraneous_options(user_options, routine):
# 
#     options = {option: user_options[option] for option in _default_options[routine]}
# 
#     return options

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


