import socket
import uuid
import getpass
import json
import os
from pkg_resources import resource_filename
import fragmenter
from fragmenter import fragment, torsions, utils, chemi
from openeye import oechem
from cmiles import to_molecule_id
import qcportal as portal


class WorkFlow(object):

    def __init__(self, workflow_id, client, workflow_json=None):
        """

        Parameters
        ----------
        id
        client

        Returns
        -------

        """
        self.workflow_id = workflow_id

        if workflow_json is not None:
            with open(workflow_json) as file:
                workflow_json = json.load(file)[workflow_id]
        # Check that all fields exist

        # Check if id already in database
        try:
            off_workflow = portal.collections.OpenFFWorkflow.from_server(client, workflow_id)
            if workflow_json is not None:
                # Check if workflows are the same
                _check_workflow(workflow_json, off_workflow)
                utils.logger().warning("The workflow ID provided already exits in the database. The options you "
                                       "provided are the same as in the database. The database options will be used.")
        except KeyError:
            # Get workflow from json file and register
            off_workflow = portal.collections.OpenFFWorkflow(workflow_id, options=workflow_json, client=client)

        self.off_workflow = off_workflow


def enumerate_states(molecule, workflow_id, options=None, title='', json_filename=None):
    """
    enumerate protonation, tautomers and stereoisomers for molecule.

    Parameters
    ----------
    molecule: any format that OpenEye pareses. Can be path to file containing molecule or SMILES/Inchi string
    workflow_id: str
        Which workflow to use as defined in data/workflows.json
    options: dict, optional, default None
        dictionary of keyword options. Default is None. If None, will use options defined in workflow ID
    title: str, optional, default empty string
        title of molecule. If None, the title of the molecule will be the IUPAC name
    json_filename: str, optional, default None
        json filename for states generated. If None will not write json file

    Returns
    -------
    json_dict: dict
        dictionary containing canonical isomeric SMILES for states and provenance.

    """
    # Load options for enumerate states
    routine = 'enumerate_states'
    provenance = _get_provenance(workflow_id=workflow_id, routine=routine)
    if not options:
        options = _get_options(workflow_id, routine)

    molecule = chemi.standardize_molecule(molecule, title=title)
    can_iso_smiles = oechem.OEMolToSmiles(molecule)
    states = fragment.expand_states(molecule, **options)

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
    workflow_id: str
        Which workflow to use for options.
    options: dictionary, optional, default None
        Dictionary of keyword options. If None, will use optiond defined in workflows
    title: str, optional. Default empty str
        The title or name of the molecule. If empty stirng will use the IUPAC name for molecule title.
    mol_provenance: dict, optional. Default is None
        provenance for molecule. If the molecule is a state from enumerate_states, the provenance from enumerate_states
        should be used
    json_filename: str, optional. Default None
        If a filename is provided, will write output to json file.
    generate_vis: bool, optional, default False
        If True, will generate visualization of fragments from parent molecule

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
    Generate input for torsiondrive QCFractal portal

    Parameters
    ----------
    fragment_dict: dict
        dictionary with fragment identifiers and provenance
    workflow_id: str
        workflow to use for options
    options: dict, optional, default None
        Keyword options. If None will use options defined in workflow
    json_filename: str, optional, default None
        If given will write jobs to json file

    Returns
    -------
    torsiondrive_inputs: dictionary defining the molecule and torsiondrive job options.

    """

    options = _get_options(workflow_id=workflow_id, routine='torsiondrive_input')
    provenance = _get_provenance(workflow_id=workflow_id, routine='torsiondrive_input')
    fragment_dict['provenance']['routine']['torsiondrive_input'] = provenance['routine']['torsiondrive_input']
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

    identifier = mol_id['canonical_isomeric_explicit_hydrogen_mapped_smiles']
    torsiondrive_inputs = {identifier: {'torsiondrive_input': {}, 'provenance': provenance}}
    needed_torsion_drives = torsions.define_torsiondrive_jobs(torsions.find_torsions(mapped_mol), **options)
    for i, job in enumerate(needed_torsion_drives):
        torsiondrive_input = dict()
        torsiondrive_input['initial_molecule'] = qm_mol
        torsiondrive_input['initial_molecule']['identifiers'] = mol_id
        torsiondrive_input['dihedrals'] = needed_torsion_drives[job]['dihedrals']
        torsiondrive_input['grid_spacing'] = needed_torsion_drives[job]['grid_spacing']
        torsiondrive_inputs[identifier]['torsiondrive_input']['job_{}'.format(i)] = torsiondrive_input

    if json_filename:
        with open(json_filename, 'w') as f:
            json.dump(torsiondrive_inputs, f, indent=2, sort_keys=True)

    return torsiondrive_inputs


def workflow(molecules_smiles, workflow_id, molecule_titles=None, generate_vis=False, write_json_intermediate=False,
             json_filename=None):
    """
    Convenience function to run Fragmenter workflow.

    Parameters
    ----------
    molecules_smiles: list of str
        list of SMILES
    molecule_titles: list of molecule names, optional. Default None
        If list of molecule names is provided, the name will be added to provenance and used for intermediate filenames.
        If not, the name will be the IUPAC name
    generate_vis: bool, optional, default None
        If True, will generate visualization of fragments
    write_json_intermediate: bool, optional. Default False
        If True will write JSON files for intermediate steps (enumerating states and fragments)
    json_filename: str, optional, default None
        filename to write jobs out to.

    Returns
    -------
    molecules: dict
        JSON specs for torsiondrive jobs. Keys are canonical isomeric explicit hydrogen mapped SMILES of fragment.

    """

    # Check input
    if not isinstance(molecules_smiles, list):
        molecules_smiles = [molecules_smiles]

    all_frags = {}
    options = _get_options(workflow_id=workflow_id, all=True)
    for i, molecule_smile in enumerate(molecules_smiles):
        filename = None
        title = ''
        if molecule_titles:
            title = molecule_titles[i]
        if write_json_intermediate and title:
            filename = 'states_{}.json'.format(title)
        if write_json_intermediate and not title:
            filename = 'states_{}.json'.format(utils.make_python_identifier(molecule_smile)[0])
        states = enumerate_states(molecule_smile, workflow_id=workflow_id, title=title, options=options['enumerate_states']['options'], json_filename=filename)
        for j, state in enumerate(states['states']):
            if write_json_intermediate and title:
                filename = 'fragments_{}_{}.json'.format(title, j)
            if write_json_intermediate and not title:
                filename = 'fragment_{}_{}.json'.format(utils.make_python_identifier(state)[0], j)
            fragments = enumerate_fragments(state, title=title, workflow_id=workflow_id,
                                            mol_provenance=states['provenance'], options=options['enumerate_fragments']['options'],
                                            json_filename=filename, generate_vis=generate_vis)
            all_frags.update(**fragments)

    all_jobs = {}
    for frag in all_frags:
        crank_jobs = generate_torsiondrive_input(all_frags[frag], workflow_id=workflow_id, options=options['torsiondrive_input']['options'])
        all_jobs.update(crank_jobs)

    if json_filename:
        with open(json_filename, 'w') as f:
            json.dump(all_jobs, f, indent=2, sort_keys=True)
    return all_jobs


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

    return provenance


def _get_options(workflow_id, routine=None, all=False):

    fn = resource_filename('fragmenter', os.path.join('data', 'workflows.json'))
    f = open(fn)
    options = json.load(f)[workflow_id]
    f.close()
    if routine:
        options = options['fragmenter'][routine]['options']
    elif all:
        options = options['fragmenter']
    else:
        raise RuntimeError("You must define a routine or set all to True for options for all fragmenter routines. ")

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


def register_workflow():
    pass

def pull_workflow():
    pass


def _check_workflow(workflow_json, off_workflow):

    for key in workflow_json:
        print(key)
        print(workflow_json[key])
        print(off_workflow.get_options(key))
        if workflow_json[key] != off_workflow.get_options(key):
            raise ValueError("The workflow ID provided already exists in the database. The options for {} are different"
                             "in the registered workflow and provided workflow. The options provided are {} and "
                             "the options in the database are {}".format(key, workflow_json[key], off_workflow.get_options(key), key))