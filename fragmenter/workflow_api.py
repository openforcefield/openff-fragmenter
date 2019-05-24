import socket
import uuid
import getpass
import json
import numpy as np
import fragmenter
from fragmenter import fragment, torsions, utils, chemi
from cmiles import get_molecule_ids
from cmiles.utils import mol_to_smiles, mol_to_map_ordered_qcschema, has_atom_map
from openeye import oechem
import copy
import warnings
#import qcportal as portal
# For Travis CI
try:
    import qcfractal.interface as portal
except ImportError:
    pass

class WorkFlow(object):

    def __init__(self, workflow_id, client, workflow_json=None, verbose=False):
        """

        Parameters
        ----------
        id
        client

        Returns
        -------

        """
        self.workflow_id = workflow_id
        self.verbose = verbose

        if workflow_json is not None:
            with open(workflow_json) as file:
                workflow_json = json.load(file)[workflow_id]['fragmenter']
        # Check that all fields exist

        # Check if id already in database
        try:
            off_workflow = client.get_collection('OpenFFWorkflow', workflow_id)
            if workflow_json is not None:
                # Check if workflows are the same
                _check_workflow(workflow_json, off_workflow)
                utils.logger().warning("The workflow ID provided already exits in the database. The options you "
                                       "provided are the same as in the database. The database options will be used.")
        except KeyError:
            # Get workflow from json file and register
            off_workflow = portal.collections.OpenFFWorkflow(workflow_id, client, **workflow_json)

        self.off_workflow = off_workflow
        self.states = {}
        self.fragments = {}
        self.qcfractal_jobs = {}
        self.final_energies = {}
        self.final_geometries = {}
        self.failed_jobs = {}

    def enumerate_states(self, molecule, title='', json_filename=None):
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
        options = self.off_workflow.get_options('enumerate_states')['options']
        provenance = _get_provenance(workflow_id=self.workflow_id, routine=routine)
        # if not options:
        #     options = _get_options(workflow_id, routine)

        molecule = chemi.standardize_molecule(molecule, title=title)
        if not options['stereoisomers']:
            can_smiles = mol_to_smiles(molecule, isomeric=True, mapped=False, explicit_hydrogen=False)
        else:
            can_smiles = mol_to_smiles(molecule, isomeric=False, mapped=False, explicit_hydrogen=False)
        states = fragment.expand_states(molecule, **options)

        provenance['routine']['enumerate_states']['parent_molecule'] = can_smiles
        provenance['routine']['enumerate_states']['parent_molecule_name'] = molecule.GetTitle()
        json_dict = {'provenance': provenance, 'states': states}

        if json_filename:
            json_dict['states'] = list(json_dict['states'])
            with open(json_filename, 'w') as f:
                json.dump(json_dict, f, indent=2, sort_keys=True)

        return json_dict

    def enumerate_fragments(self, molecule, title='', mol_provenance=None, json_filename=None,
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
        provenance = _get_provenance(workflow_id=self.workflow_id, routine=routine)
        options = self.off_workflow.get_options('enumerate_fragments')
        scheme = options['scheme']
        options = options['options']

        parent_molecule = chemi.standardize_molecule(molecule, title)
        parent_molecule_smiles = mol_to_smiles(parent_molecule, isomeric=True, explicit_hydrogen=False,
                                                        mapped=False)
        provenance['routine']['enumerate_fragments']['parent_molecule_name'] = parent_molecule.GetTitle()
        provenance['routine']['enumerate_fragments']['parent_molecule'] = parent_molecule_smiles

        functional_groups = options.pop('functional_groups')
        if scheme == 'combinatorial':
            fragment_engine = fragment.CombinatorialFragmenter(parent_molecule, functional_groups=functional_groups)
        elif scheme == 'wiberg_bond_order':
            fragment_engine = fragment.WBOFragmenter(parent_molecule, functional_groups=functional_groups)
        else:
            raise ValueError("Only combinatorial and wiberg_bond_order are supported fragmenting schemes")

        fragment_engine.fragment(**options)

        if generate_vis:
            fname = '{}.pdf'.format(parent_molecule.GetTitle())
            fragment_engine.depict_fragments(fname=fname)

        if not fragment_engine.fragments:
            warnings.warn("No fragments were generated for {}".format(parent_molecule_smiles))
            # ToDo: if combinatorial does not have fragments it means that there was no point to cut and if
            # wbo does not have fragments it means no internal rotatable bond was found but it might still have torsions
            # to drive. Not yet sure how to handle these cases
            return

        if self.states:
            # Check if current state exists
            if parent_molecule_smiles in self.states['states']:
                provenance['routine']['enumerate_states'] = self.states['provenance']['routine']['enumerate_states']
            elif mol_provenance:
                provenance['routine']['enumerate_states'] = mol_provenance['routine']['enumerate_states']

        # Generate identifiers for fragments
        fragments_json_dict = {}
        for frag in fragment_engine.fragments:
            mols = fragment_engine.fragments[frag]
            if not isinstance(mols, list):
                mols = [mols]
            for mol in mols:
                parent_map_smiles = oechem.OEMolToSmiles(mol)
                smiles = mol_to_smiles(mol, mapped=False)
                identifiers = get_molecule_ids(smiles, canonicalization='openeye')
                frag_smiles = identifiers['canonical_isomeric_smiles']
                fragments_json_dict[frag_smiles] = {'identifiers': identifiers}
                fragments_json_dict[frag_smiles]['provenance'] = provenance
                fragments_json_dict[frag_smiles]['provenance']['canonicalization'] = identifiers.pop('provenance')
                if has_atom_map(mol):
                    fragments_json_dict[frag_smiles]['provenance']['routine']['enumerate_fragments']['map_to_parent'] \
                        = parent_map_smiles
                if scheme == 'wiberg_bond_order':
                    # Add bond to provenance
                    fragments_json_dict[frag_smiles]['provenance']['routine']['enumerate_fragments']['rot_bond'] = frag

        if json_filename:
            with open(json_filename, 'w') as f:
                json.dump(fragments_json_dict, f, indent=2, sort_keys=True)

        return fragments_json_dict

    def generate_torsiondrive_input(self, frag, json_filename=None):
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

        options = self.off_workflow.get_options('torsiondrive_input')['options']
        provenance = _get_provenance(workflow_id=self.workflow_id, routine='torsiondrive_input')
        frag['provenance']['routine']['torsiondrive_input'] = provenance['routine']['torsiondrive_input']
        provenance = frag['provenance']

        mol_id = frag['identifiers']

        mapped_smiles = mol_id['canonical_isomeric_explicit_hydrogen_mapped_smiles']
        mapped_mol = chemi.smiles_to_oemol(mapped_smiles)
        print(options)
        needed_torsions = torsions.find_torsions(mapped_mol, options['restricted'])

        if 'conf_grid' in options:
            # Generate grid of multiple conformers
            dihedrals = []
            for torsion_type in needed_torsions:
                for tor in needed_torsions[torsion_type]:
                    dihedrals.append(needed_torsions[torsion_type][tor])
            intervals = options['initial_conf_grid_resolution']
            if not isinstance(intervals, list):
                intervals = [intervals]*len(dihedrals)
            try:
                conformers = chemi.generate_grid_conformers(mapped_mol, dihedrals=dihedrals, intervals=intervals)
            except RuntimeError:
                utils.logger().warning("{} does not have coordinates. This can happen for several reasons related to Omega. "
                              "{} will not be included in fragments dictionary".format(
                        mol_id['canonical_isomeric_smiles'], mol_id['canonical_isomeric_smiles']))
                return False

            chemi.resolve_clashes(conformers)
        else:
            try:
                max_confs = options['torsiondrive_options'].pop('max_confs')
                conformers = chemi.generate_conformers(mapped_mol, max_confs=max_confs)
            except RuntimeError:
                utils.logger().warning("{} does not have coordinates. This can happen for several reasons related to Omega. "
                                  "{} will not be included in fragments dictionary".format(
                            mol_id['canonical_isomeric_smiles'], mol_id['canonical_isomeric_smiles']))
                return False

        qcschema_molecules = [mol_to_map_ordered_qcschema(conf, mol_id) for conf in conformers.GetConfs()]
        identifier = mol_id['canonical_isomeric_smiles']
        torsiondrive_inputs = {identifier: {'torsiondrive_input': {}, 'provenance': provenance}}
        restricted_torsions = needed_torsions.pop('restricted')
        if restricted_torsions:
            optimization_jobs = torsions.generate_constraint_opt_input(qcschema_molecules, restricted_torsions,
                                                             **options['restricted_optimization_options'])
            torsiondrive_inputs[identifier]['optimization_input'] = optimization_jobs
        torsiondrive_jobs = torsions.define_torsiondrive_jobs(needed_torsions, **options['torsiondrive_options'])

        #if options['multiple_confs']:
        #    qcschema_molecule = qcschema_molecules

        # Currently, all jobs are started from same initial conformation
        # ToDo Start later job from optimized conformers from last job
        for i, job in enumerate(torsiondrive_jobs):
            torsiondrive_input = {'type': 'torsiondrive_input'}
            torsiondrive_input['initial_molecule'] = qcschema_molecules
            #torsiondrive_input['initial_molecule']['identifiers'] = mol_id
            torsiondrive_input['dihedrals'] = torsiondrive_jobs[job]['dihedrals']
            torsiondrive_input['grid_spacing'] = torsiondrive_jobs[job]['grid_spacing']
            job_name = ''
            mapped_smiles = qcschema_molecules[0]['identifiers']['canonical_isomeric_explicit_hydrogen_mapped_smiles']
            for i, torsion in enumerate(torsiondrive_input['dihedrals']):
                label = to_canonical_label(mapped_smiles, torsion)
                if i > 0:
                    job_name += ' {}'.format(label)
                else:
                    job_name += '{}'.format(label)

            torsiondrive_inputs[identifier]['torsiondrive_input'][job_name] = torsiondrive_input

        if json_filename:
            with open(json_filename, 'w') as f:
                json.dump(torsiondrive_inputs, f, indent=2, sort_keys=True)

        return torsiondrive_inputs

    def workflow(self, molecules_smiles, molecule_titles=None, generate_vis=False, write_json_intermediate=False,
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
        if molecule_titles:
            if not isinstance(molecule_titles, list):
                molecule_titles = [molecule_titles]

        all_frags = {}

        for i, molecule_smile in enumerate(molecules_smiles):
            filename = None
            title = ''
            if molecule_titles:
                title = molecule_titles[i]
            if write_json_intermediate and title:
                filename = 'states_{}.json'.format(title)
            if write_json_intermediate and not title:
                filename = 'states_{}.json'.format(utils.make_python_identifier(molecule_smile)[0])
            self.states = self.enumerate_states(molecule_smile, title=title, json_filename=filename)
            for j, state in enumerate(self.states['states']):
                if write_json_intermediate and title:
                    filename = 'fragments_{}_{}.json'.format(title, j)
                if write_json_intermediate and not title:
                    filename = 'fragment_{}_{}.json'.format(utils.make_python_identifier(state)[0], j)
                fragments = self.enumerate_fragments(state, title=title, mol_provenance=self.states['provenance'],
                                                json_filename=filename, generate_vis=generate_vis)
                if not fragments:
                    continue
                all_frags.update(**fragments)
        self.fragments = all_frags

        all_jobs = {}
        for frag in all_frags:
            crank_jobs = self.generate_torsiondrive_input(all_frags[frag])
            if not crank_jobs:
                continue
            all_jobs.update(crank_jobs)
        self.qcfractal_jobs = all_jobs

        if json_filename:
            with open(json_filename, 'w') as f:
                json.dump(all_jobs, f, indent=2, sort_keys=True)
        #return all_jobs

    def add_fragments_to_db(self):
        for frag in self.qcfractal_jobs:
            torsiondrive_input = self.qcfractal_jobs[frag]['torsiondrive_input']
            optimization_input = self.qcfractal_jobs[frag]['optimization_input']
            # combine both dictionaries
            input_data = {**torsiondrive_input, **optimization_input}
            if input_data:
                self.off_workflow.add_fragment(frag, input_data, self.qcfractal_jobs[frag]['provenance'])

    def add_fragments_from_json(self, json_filenam):
        with open(json_filenam) as f:
            self.qcfractal_jobs = json.load(f)
        self.add_fragments_to_db()

    def get_final_molecules(self, json_filename=None):
        """
        Get final molecule geometries and energies from db and serialize keys for JSON
        Parameters
        ----------
        json_filename: str, optional. Default None
            If name is given, final energies and geometries will be written out to a JSON file

        """
        final_energies = copy.deepcopy(self.off_workflow.list_final_energies())
        self.final_energies = self._to_json_format(final_energies)
        if json_filename:
            filename = '{}_energies.json'.format(json_filename)
            with open(filename, 'w') as f:
                json.dump(self.final_energies, f, indent=2, sort_keys=True)

        if self.verbose:
            utils.logger().info("Pulling final geometries from database. This takes some time...")
        final_geometries = copy.deepcopy(self.off_workflow.list_final_molecules())
        self.final_geometries = self._to_json_format(final_geometries)
        if json_filename:
            filename = '{}_geometries.json'.format(json_filename)
            with open(filename, 'w') as f:
                json.dump(self.final_geometries, f, indent=2, sort_keys=True)

    def _to_json_format(self, final_dict):

        serialized_dict = {}
        for frag in final_dict:
            for job in final_dict[frag]:
                if not final_dict[frag][job]:
                    # job failed.
                    if not frag in self.failed_jobs:
                        self.failed_jobs[frag] = []
                        if not job in self.failed_jobs[frag]:
                            self.failed_jobs[frag].append(job)
                    continue
                if not isinstance(final_dict[frag][job], dict):
                    # This is an optimization job. No serialization needed
                    if frag not in serialized_dict:
                        serialized_dict[frag] = {}
                    serialized_dict[frag][job] = final_dict[frag][job]
                    continue

                for key in final_dict[frag][job]:
                    value = final_dict[frag][job][key]
                    new_key = serialize_key(key)
                    if frag not in serialized_dict:
                        serialized_dict[frag] = {}
                    if job not in serialized_dict[frag]:
                        serialized_dict[frag][job] = {}
                    serialized_dict[frag][job][new_key] = value
        return serialized_dict



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

    provenance = {'creator': fragmenter.__package__,
                  'job_id': str(uuid.uuid4()),
                  'hostname': socket.gethostname(),
                  'username': getpass.getuser(),
                  'workflow_id': workflow_id,
                  'routine': {routine: {
                      'version': fragmenter.__version__
                  }}}

    return provenance

def _check_workflow(workflow_json, off_workflow):
    #ToDo this does not always raise an error (it raises one for enumerate states but not enumerate fragments)
    for key in workflow_json:
        if workflow_json[key] != off_workflow.get_options(key):
            raise ValueError("The workflow ID provided already exists in the database. The options for {} are different "
                             "in the registered workflow and provided workflow. The options provided are {} and "
                             "the options in the database are {}".format(key, workflow_json[key], off_workflow.get_options(key), key))
        else:
            return True


def serialize_key(key):
        if isinstance(key, (int, float)):
            key = (int(key), )

        return json.dumps(key)

def grid_id_from_str(grid_id_str):
    """
    Only works for 1D grids.
    Deserialize grid ID key

    Parameters
    ----------
    grid_id_str

    Returns
    -------

    """
    return int(grid_id_str.split('[')[-1].split(']')[0])

def serialize_fractal_output(final_json):
    """

    Parameters
    ----------
    final_json

    Returns
    -------

    """
    import qcelemental

    serialized_dict = {}
    for frag in final_json:
        serialized_dict[frag] = {}
        for job in final_json[frag]:
            serialized_dict[frag][job] = {}
            for key in final_json[frag][job]:
                value = final_json[frag][job][key]
                new_key = serialize_key(key)
                if isinstance(value, qcelemental.models.molecule.Molecule):
                    value = json.loads(value.json())
                serialized_dict[frag][job][new_key] = value
    return serialized_dict


def sort_energies(final_energies):
    """
    Sort energies by angle in place

    Parameters
    ----------
    final_energies: dictionary
        output from workflow.list_final_energies()


    """
    for frag in final_energies:
        for job in final_energies[frag]:
            angles = []
            energies = []
            for angle in final_energies[frag][job]:
                energy = final_energies[frag][job][angle]
                angle = grid_id_from_str(angle)
                angles.append(angle)
                energies.append(energy)
            energies = np.asarray(energies)
            energies = energies * utils.HARTREE_2_KJMOL
            rel_energies = energies - energies.min()
            sorted_energies = [x for _, x in sorted(zip(angles, rel_energies))]
            sorted_angles = sorted(angles)
            final_energies[frag][job] = (sorted_angles, sorted_energies)


def deserialze_molecules(final_molecules):
    deserialized = {}
    for frag in final_molecules:
        deserialized[frag] = {}
        for job in final_molecules[frag]:
            deserialized[frag][job] = {}
            for angle in final_molecules[frag][job]:
                molecule = final_molecules[frag][job][angle]
                angle_int = grid_id_from_str(angle)
                deserialized[frag][job][angle_int] = molecule
    return deserialized


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

