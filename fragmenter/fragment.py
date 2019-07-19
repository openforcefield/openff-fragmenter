from itertools import combinations
import openeye as oe
from openeye import oechem, oedepict, oegrapheme, oequacpac, oeomega
import cmiles.utils
from cmiles.utils import mol_to_smiles, has_stereo_defined, has_atom_map, is_missing_atom_map, remove_atom_map, restore_atom_map

import yaml
import os
from pkg_resources import resource_filename

import warnings
import networkx as nx
import time

from .utils import logger
from .chemi import get_charges, generate_conformers, LabelWibergBondOrder, LabelFragBondOrder, ColorAtomByFragmentIndex
from .torsions import find_torsion_around_bond


OPENEYE_VERSION = oe.__name__ + '-v' + oe.__version__



def enumerate_states(molecule, tautomers=True, stereoisomers=True, verbose=False, return_mols=False,
                     explicit_h=True, return_names=False, max_stereo_returns=1, filter_nitro=True, **kwargs):
    """
    Expand tautomeric state and stereoisomers for molecule.

    Parameters
    ----------
    molecule : OEMol
    tautomers : bool, optional, default True
        If False, will not generate tautomers
    stereoisomers : bool, optional, default True
        If False, will not generate all stereoisomers.
    verbose : bool, optional, default False
    return_mols : bool, optional, default False
        If True, will return oemols instead of SMILES. Some molecules might be duplicate states
    explicit_h : bool, optional, default True
        If True, SMILES of states will have explicit hydrogen
    return_names : bool, optional, default True
        If True, will return names of molecules with SMILES
    max_stereo_returns : int, optional, default 1
        If stereoisomers is set to False, and the incoming molecule is missing stereo information, OEFlipper will
        generate stereoisomers for missing stereo center. max_stereo_returns controls how many of those will be returned
    ** max_states: int, optional, default 200
        This gets passed to `_enumerate_tautomers` and `_enumerate_stereoisomers`
        max number of states `_enumerate_tautomers` and `_enumerate_stereoisomers` generate
    ** pka_norm: bool, optional, default True
        This gets passed to `_enumerate_tautomers`. If True, ionization state of each tautomer will be assigned to a predominate
        state at pH ~7.4
    ** warts: bool, optional, default True
        This gets passed to `_enumerate_tautomers` and _enumerate_stereoisomers`
        If True, adds a wart to each new state. A 'wart' is a systematic
    ** force_flip: bool, optional, default True
        This gets passed to `_enumerate_stereoisomers`
        Force flipping all stereocenters. If False, will only generate stereoisomers for stereocenters that are undefined
    ** enum_nitorgen: bool, optional, default True
        This gets passed to `_enumerate_stereoisomers`
        If true, invert non-planer nitrogens

    Returns
    -------
    states: list
        list of oemols or SMILES of states generated for molecule

    """

    # If incoming molecule has nitro in form ([NX3](=O)=O), do not filter out later
    if _check_nitro(molecule):
        filter_nitro = False
    title = molecule.GetTitle()
    states = []
    if return_names:
        names = []

    if verbose:
        logger().info("Enumerating states for {}".format(title))

    if stereoisomers:
        if verbose:
            logger().info("Enumerating stereoisomers for {}".format(title))
        stereo_mols = (_enumerate_stereoisomers(molecule, **kwargs))
        if verbose:
            logger().info('Enumerated {} stereoisomers'.format(len(stereo_mols)))

    if tautomers:
        if not stereoisomers:
            stereo_mols = [molecule]
        tau_mols = []
        if verbose:
            logger().info("Enumerating tautomers states for {}".format(title))
        for mol in stereo_mols:
            tau_mols.extend(_enumerate_tautomers(mol, **kwargs))
        if verbose:
            logger().info('Enumerated {} tautomers'.format(len(tau_mols)))

        # check for nitro in ([NX3](=O)=O) form
        if filter_nitro:
            tau_mols[:] = [mol for mol in tau_mols if not _check_nitro(mol)]

    if stereoisomers and tautomers:
        all_mols = stereo_mols + tau_mols
    elif stereoisomers and not tautomers:
        all_mols = stereo_mols
    elif not stereoisomers and tautomers:
        all_mols = tau_mols
        all_mols.append(molecule)
    else:
        all_mols = [molecule]

    if return_mols:
        return all_mols

    for mol in all_mols:
        try:
            smiles = mol_to_smiles(mol, isomeric=True, mapped=False, explicit_hydrogen=explicit_h)
            if smiles not in states:
                states.append(smiles)
                if return_names:
                    names.append(mol.GetTitle())

        except ValueError:
            # Stereo is not fully defined. Use flipper with force_flip set to False
            stereo_states = _enumerate_stereoisomers(mol, force_flip=False, enum_nitrogen=True, warts=True)
            if len(stereo_states) > max_stereo_returns:
                stereo_states = stereo_states[:max_stereo_returns]

            for state in stereo_states:
                try:
                    smiles = mol_to_smiles(state, isomeric=True, mapped=False, explicit_hydrogen=explicit_h)
                except ValueError:
                    stereo_states_forced = _enumerate_stereoisomers(mol, force_flip=True, enum_nitrogen=True, warts=True)
                    if len(stereo_states_forced) > max_stereo_returns:
                        stereo_states_forced = stereo_states_forced[:max_stereo_returns]
                    for state_forced in stereo_states_forced:
                        smiles = mol_to_smiles(state_forced, isomeric=True, mapped=False, explicit_hydrogen=explicit_h)
                        if smiles not in states:
                            states.append(smiles)
                            if return_names:
                                names.append(state.GetTitle())
                if smiles not in states:
                    states.append(smiles)
                    if return_names:
                        names.append(state.GetTitle())

    if verbose:
        logger().info("{} states were generated for {}".format(len(states), oechem.OEMolToSmiles(molecule)))

    if return_names:
        return states, names

    return states

def _enumerate_tautomers(molecule, max_states=200, pka_norm=True, warts=True):
    """
    Expand reasonable tautomer states. This function generates tautomers (which might be different ionization states
    than parent) that are normalized to the predominant state at pH ~7.4
    Parameters
    ----------
    molecule : OEMol to expand states
    max_states : int
        max number of states
    pka_norm: bool, optional, default True
    warts: bool, optional default True

    Returns
    -------
    tautomers: list of oemols

    """
    tautomers = []
    tautomer_options = oequacpac.OETautomerOptions()
    tautomer_options.SetApplyWarts(warts)
    tautomer_options.SetMaxTautomersGenerated(max_states)
    i = 0
    for tautomer in oequacpac.OEGetReasonableTautomers(molecule, tautomer_options, pka_norm):
        i += 1
        tautomers.append(tautomer)
    return tautomers

def _enumerate_stereoisomers(molecule, max_states=200, force_flip=True, enum_nitrogen=True, warts=True, verbose=True):
    """
    Enumerate stereoisomers
    Parameters
    ----------
    molecule : OEMol
    max_states : int, optional, default 200
        max number of states to enumerate
    force_flip : bool, optional, default True
        If True, will flip all steocenters. If False, will only flip centers that are undefined
    enum_nitrogen : bool, optional, default True
        Invert non-planar nitrogen
    warts : bool, optional, default True
        If True, add int to molecule name
    verbose : bool, optional, default True

    Returns
    -------
    stereoisomers: list of oemols

    """
    stereoisomers = []
    if verbose:
        logger().debug("Enumerating stereoisomers...")
    i = 0
    for enantiomer in oeomega.OEFlipper(molecule, max_states, force_flip, enum_nitrogen, warts):
        i += 1
        enantiomer = oechem.OEMol(enantiomer)
        stereoisomers.append(enantiomer)
    return stereoisomers

def _check_nitro(molecule):
    """
    Filter out nitro that is in ([NX3](=O)=O) form. OEGetReasonableTautomers generates this form.
    Parameters
    ----------
    molecule :

    Returns
    -------

    """
    qmol = oechem.OEQMol()
    if not oechem.OEParseSmarts(qmol, '([NX3](=O)=O)'):
        print('OEParseSmarts failed')
    ss = oechem.OESubSearch(qmol)
    oechem.OEPrepareSearch(molecule, ss)
    matches = [m for m in ss.Match(molecule)]
    return bool(matches)


class Fragmenter(object):
    """
    Base fragmenter class.
    This class is inherited by CombinatorialFragmenter and WBOFragmenter
    """

    def __init__(self, molecule):
        self.molecule = molecule
        self.functional_groups = {}

         # Add explicit hydrogen
        if not cmiles.utils.has_explicit_hydrogen(self.molecule):
            oechem.OEAddExplicitHydrogens(self.molecule)
        # Add canonical atom map
        cmiles.utils.add_atom_map(self.molecule)

        if is_missing_atom_map(self.molecule):
            warnings.warn("Some atoms are missing atom maps. This might cause a problem at several points during "
                          "the fragmentation process. Make sure you know what you are doing. ")

        # Keep track of stereo to make sure it does not flip
        self.stereo = {}
        self._atom_stereo_map = {oechem.OECIPAtomStereo_S: 'S', oechem.OECIPAtomStereo_R: 'R'}
        self._bond_stereo_map = {oechem.OECIPBondStereo_E: 'E', oechem.OECIPBondStereo_Z: 'Z'}
        self._find_stereo()

        # keep track of fragments that form new stereocenters
        self.new_stereo = []

        # For provenance
        self._options = {}
    @property
    def n_rotors(self):
        return sum([bond.IsRotor() for bond in self.molecule.GetBonds()])

    @staticmethod
    def _count_rotors_in_fragment(fragment):
        return sum([bond.IsRotor() for bond in fragment.GetBonds()])

    @staticmethod
    def _count_heavy_atoms_in_fragment(fragment):
        return sum([not atom.IsHydrogen() for atom in fragment.GetAtoms()])

    def _find_stereo(self):
        """
        Find chiral atoms and bonds, store the chirality.
        This is needed to check if fragments flipped chirality. Currently this can happen and it
        is a bug

        """
        for a in self.molecule.GetAtoms():
            if a.IsChiral():
                self.stereo[a.GetMapIdx()] = self._atom_stereo_map[oechem.OEPerceiveCIPStereo(self.molecule, a)]
        for bond in self.molecule.GetBonds():
            if bond.IsChiral():
                m1 = bond.GetBgn().GetMapIdx()
                m2 = bond.GetEnd().GetMapIdx()
                self.stereo[(m1, m2)] = self._bond_stereo_map[oechem.OEPerceiveCIPStereo(self.molecule, bond)]

    def _check_stereo(self, fragment, use_parent=False):
        """
        Check if stereo in fragment is different than stereo in parent. If stereo flips, it raises
        an exception.

        Parameters
        ----------
        fragment : oemol or AtomBondSet
        use_parent : bool, optional, default False
            When checking if an AtomBondSet has the same stereo, the parent molecule needs to be used.

        """
        for a in fragment.GetAtoms():
            if a.IsChiral():
                if use_parent:
                    s = self._atom_stereo_map[oechem.OEPerceiveCIPStereo(self.molecule, a)]
                else:
                    s = self._atom_stereo_map[oechem.OEPerceiveCIPStereo(fragment, a)]
                if not a.GetMapIdx() in self.stereo:
                    logger().warning('A new stereocenter formed at atom {} {}'.format(a.GetMapIdx(),
                                                                                      oechem.OEGetAtomicSymbol(a.GetAtomicNum())))
                    return False
                if not s == self.stereo[a.GetMapIdx()]:
                    logger().warning('Stereochemistry for atom {} flipped from {} to {}'.format(a.GetMapIdx(), self.stereo[a.GetMapIdx()],
                                                                                            s))
                    return False
        for b in fragment.GetBonds():
            if b.IsChiral():
                b_tuple = (b.GetBgn().GetMapIdx(), b.GetEnd().GetMapIdx())
                if use_parent:
                    s = self._bond_stereo_map[oechem.OEPerceiveCIPStereo(self.molecule, b)]
                else:
                    s = self._bond_stereo_map[oechem.OEPerceiveCIPStereo(fragment, b)]
                if not s == self.stereo[b_tuple]:
                    logger().warning('Stereochemistry fro bond {} flipped from {} to {}'.format(b_tuple, self.stereo[b_tuple],
                                     s))
                    return False
        return True

    def _fix_stereo(self, fragment):
        """
        Flip all stereocenters and find the stereoisomer that matches the parent

        Parameters
        ----------
        fragment : oemol

        Returns
        -------
        mol: oemol with same stereochemistry as parent molecule

        """
        all_stereo = _enumerate_stereoisomers(fragment, force_flip=True)
        for mol in all_stereo:
            if self._check_stereo(mol):
                return mol
        # Could not fix stereo because it is a new chiral center. Returning all stereoisomers around new stereocenter
        return _enumerate_stereoisomers(fragment, force_flip=False)


    def _tag_functional_groups(self, functional_groups):
        """
        This function tags atoms and bonds of functional groups defined in fgroup_smarts. fgroup_smarts is a dictionary
        that maps functional groups to their smarts pattern. It can be user generated or from yaml file.

        Parameters
        ----------
        functional_groups: dict, optional, default None.
            fgroup_smarts is a dictionary of SMARTS of functional groups that should not be fragmented.
            It can be user generated.
            If it is None, fragmenter/fragmenter/data/fgroup_smarts.yaml will be used.
            If False, no functional groups will be tagged and they will all be fragmented.

        """
        if functional_groups:
            fgroups_smarts = functional_groups
        if functional_groups is None:
            # Load yaml file
            fn = resource_filename('fragmenter', os.path.join('data', 'fgroup_smarts.yml'))
            with open(fn, 'r') as f:
                fgroups_smarts = yaml.safe_load(f)
        elif functional_groups is False:
            # Don't tag fgroups
            return

        for f_group in fgroups_smarts:
            qmol = oechem.OEQMol()
            if not oechem.OEParseSmarts(qmol, fgroups_smarts[f_group]):
                print('OEParseSmarts failed')
            ss = oechem.OESubSearch(qmol)
            oechem.OEPrepareSearch(self.molecule, ss)

            for i, match in enumerate(ss.Match(self.molecule, True)):
                fgroup_atoms = set()
                for ma in match.GetAtoms():
                    fgroup_atoms.add(ma.target.GetMapIdx())
                    tag = oechem.OEGetTag('fgroup')
                    ma.target.SetData(tag, '{}_{}'.format(f_group, str(i)))
                fgroup_bonds = set()
                for ma in match.GetBonds():
                    #if not ma.target.IsInRing():
                    m1 = ma.target.GetBgn().GetMapIdx()
                    m2 = ma.target.GetEnd().GetMapIdx()
                    fgroup_bonds.add((m1, m2))
                    #fgroup_bonds.add(ma.target.GetIdx())
                    tag =oechem.OEGetTag('fgroup')
                    ma.target.SetData(tag, '{}_{}'.format(f_group, str(i)))

                self.functional_groups['{}_{}'.format(f_group, str(i))] = (fgroup_atoms, fgroup_bonds)

    def _to_atom_bond_set(self, atoms, bonds):
        """
        Convert sets of atom map indices and bond map indices tuples to OEAtomBondSet
        Parameters
        ----------
        atoms : set of ints
            set of map indices
        bonds : set of tuples of ints
            set of bond tuples (m1, m2)

        Returns
        -------
        atom_bond_set: OEAtomBondSet with atoms and bonds specified in atoms bonds set

        """
        atom_bond_set = oechem.OEAtomBondSet()
        for a_mdx in atoms:
            atom = self.molecule.GetAtom(oechem.OEHasMapIdx(a_mdx))
            atom_bond_set.AddAtom(atom)
        for b_tuple in bonds:
            a1 = self.molecule.GetAtom(oechem.OEHasMapIdx(b_tuple[0]))
            a2 = self.molecule.GetAtom(oechem.OEHasMapIdx(b_tuple[-1]))
            bond = self.molecule.GetBond(a1, a2)
            if not bond:
                raise ValueError("{} is a disconnected bond".format(b_tuple))
            atom_bond_set.AddBond(bond)
        self._check_stereo(atom_bond_set, use_parent=True)
        return atom_bond_set

    def _atom_bond_set_to_mol(self, frag, adjust_hcount=True):
        """
        Convert fragments (AtomBondSet) to OEMol
        Parameters
        ----------
        frag: OEAtomBondSet
        adjust_hcount: bool, optional, default True
            If False, hydrogen counts will not be adjusted. Not recommended.

        Returns
        -------
        fragment: OEMol
        """
        fragatompred = oechem.OEIsAtomMember(frag.GetAtoms())
        fragbondpred = oechem.OEIsBondMember(frag.GetBonds())

        fragment = oechem.OEMol()
        adjustHCount = adjust_hcount
        oechem.OESubsetMol(fragment, self.molecule, fragatompred, fragbondpred, adjustHCount)

        oechem.OEAddExplicitHydrogens(fragment)
        oechem.OEPerceiveChiral(fragment)
        # sanity check that all atoms are bonded
        for atom in fragment.GetAtoms():
            if not list(atom.GetBonds()):
                warnings.warn("Yikes!!! An atom that is not bonded to any other atom in the fragment. "
                              "You probably ran into a bug. Please report the input molecule to the issue tracker")
        #Always restore map?
        #if restore_maps:
        restore_atom_map(fragment)

        # Perceive stereo and check that defined stereo did not change
        oechem.OEPerceiveChiral(fragment)
        oechem.OE3DToAtomStereo(fragment)
        oechem.OE3DToBondStereo(fragment)
        if isinstance(self, CombinatorialFragmenter):
            # I only saw chirality flip when calculating WBOs - skip the fix and check for combinatorial
            # fragmentation because it can change the relative stereo
            return fragment
        if has_stereo_defined(fragment):
            if not self._check_stereo(fragment):
                fragment = self._fix_stereo(fragment)

        # check for stereo defined
        # if not has_stereo_defined(fragment) and expand_stereoisomers:
        #     # Try to convert to smiles and back. A molecule might look like it's missing stereo because of submol
        #     # first restore atom map so map on mol is not lost
        #     #restore_atom_map(fragment)
        #     new_smiles = oechem.OEMolToSmiles(fragment)
        #     fragment = oechem.OEMol()
        #     oechem.OESmilesToMol(fragment, new_smiles)
        #     # add explicit H
        #     oechem.OEAddExplicitHydrogens(fragment)
        #     oechem.OEPerceiveChiral(fragment)
        #     #remove_atom_map(fragment, keep_map_data=True)
        #     # If it's still missing stereo, expand states
        #     if not has_stereo_defined(fragment):
        #         if expand_stereoisomers:
        #             enantiomers = _expand_states(fragment, enumerate='stereoisomers')
        #             return enantiomers

        return fragment

    def get_provenance(self):
        """
        Get version of fragmenter and options used
        Returns
        -------

        """
        import fragmenter
        import uuid
        import socket
        import getpass
        fragmenter_version = fragmenter.__version__
        # Restore map to parent
        restore_atom_map(self.molecule)
        provenance = {'creator': fragmenter.__package__,
                      'job_id': str(uuid.uuid4()),
                      'hostname': socket.gethostname(),
                      'username': getpass.getuser(),
                      'routine': {'fragment_molecule': {
                          'version': fragmenter_version,
                          'options': self._options,
                          'parent_molecule': cmiles.utils.mol_to_smiles(self.molecule, mapped=False,
                                                                        explicit_hydrogen=False),
                          'parent_name': self.molecule.GetTitle(),
                          'mapped_parent_smiles': oechem.OEMolToSmiles(self.molecule)
                      }}}
        return provenance


class CombinatorialFragmenter(Fragmenter):
    """
    This fragmenter will fragment all bonds besides the ones in rings and sepcified functional
    groups. Then it will generate all possible connected fragment.
    This class should only be used to generate validation sets. It is not recommended for general
    fragmentation because it generates a lot more fragments than is needed.
    """
    def __init__(self, molecule, functional_groups=None):
        super().__init__(molecule)

        if functional_groups is None:
            fn = resource_filename('fragmenter', os.path.join('data', 'fgroup_smarts_comb.yml'))
            with open(fn, 'r') as f:
                functional_groups = yaml.safe_load(f)
        self._options['functional_groups'] = functional_groups
        self._tag_functional_groups(functional_groups)
        self._nx_graph = self._mol_to_graph()
        self._fragments = list()  # all possible fragments without breaking rings
        self._fragment_combinations = list() # AtomBondSets of combined fragments. Used internally to generate PDF
        self.fragments = {} # Dict that maps SMILES to all equal combined fragments

    def fragment(self, min_rotors=1, max_rotors=None, min_heavy_atoms=4):
        """
        combinatorial fragmentation. Fragment along every bond that is not in a ring or specified
        functional group.
        Parameters:
        ----------
        min_rotor: int, optional, default 1
            The minimum number of rotatable bond in resutling fragments
        max_rotor: int, optional, default None
            The maximum number of rotatable bond in resulting fragment
            If None, the maximum number of rotors will be the amount in the parent molecule.
        min_heavy_atoms: int, optional, default 4
            minimum number of heavy atoms in the resulting fragments

        """
        # Store options
        self._options['min_rotors'] = min_rotors
        self._options['min_heavy_atoms'] = min_heavy_atoms
        # Remove atom map and store it in data. This is needed to keep track of the fragments
        if has_atom_map(self.molecule):
            remove_atom_map(self.molecule, keep_map_data=True)
        self.fragment_all_bonds_not_in_ring_systems()
        if max_rotors is None:
            max_rotors = self.n_rotors + 1
        self._options['max_rotors'] = max_rotors
        self.combine_fragments(min_rotors=min_rotors, max_rotors=max_rotors, min_heavy_atoms=min_heavy_atoms)

    def _mol_to_graph(self):
        """
        Convert molecule to networkx graph.
        Nodes are the atoms, edges are the bonds.

        """

        G = nx.Graph()
        for atom in self.molecule.GetAtoms():
           # print(atom, atom.GetMapIdx())
            G.add_node(atom.GetIdx(), name=oechem.OEGetAtomicSymbol(atom.GetAtomicNum()), halogen=atom.IsHalogen())
        for bond in self.molecule.GetBonds():
            # Check for functional group tags
            try:
                fgroup = bond.GetData('fgroup')
            except ValueError:
                fgroup = False
            G.add_edge(bond.GetBgnIdx(), bond.GetEndIdx(),  index=bond.GetIdx(),
                       aromatic=bond.IsAromatic(), in_ring=bond.IsInRing(), bond_order=bond.GetOrder(),
                       rotor=bond.IsRotor(), fgroup=fgroup)
        return G

    def _fragment_graph(self):
        """
        Fragment all bonds that are not in rings or not tagged as a functional group

        """
        import networkx as nx

        ebunch = []
        for node in self._nx_graph.edges:
            if not self._nx_graph.edges[node]['aromatic'] and not self._nx_graph.edges[node]['in_ring'] \
                    and not (self._nx_graph.node[node[0]]['name'] == 'H' or self._nx_graph.node[node[-1]]['name'] == 'H')\
                    and not (self._nx_graph.edges[node]['fgroup']):
                ebunch.append(node)
        # Cut molecule
        self._nx_graph.remove_edges_from(ebunch)
        # Generate fragments
        subgraphs = list(nx.connected_component_subgraphs(self._nx_graph))
        return subgraphs

    def _subgraphs_to_atom_bond_sets(self, subgraphs):
        """
        Build Openeye AtomBondSet from subrgaphs for enumerating fragments recipe

        Parameters
        ----------
        subgraph: NetworkX subgraph

        """
        # Build openeye atombondset from subgraphs
        for subgraph in subgraphs:
            atom_bond_set = oechem.OEAtomBondSet()
            for node in subgraph.node:
                atom_bond_set.AddAtom(self.molecule.GetAtom(oechem.OEHasAtomIdx(node)))
            for node1, node2 in subgraph.edges():
                a1 = self.molecule.GetAtom(oechem.OEHasAtomIdx(node1))
                a2 = self.molecule.GetAtom(oechem.OEHasAtomIdx(node2))
                bond = self.molecule.GetBond(a1, a2)
                atom_bond_set.AddBond(bond)
            self._fragments.append(atom_bond_set)

    def fragment_all_bonds_not_in_ring_systems(self):
        """

        Returns
        -------

        """
        subgraphs = self._fragment_graph()
        self._subgraphs_to_atom_bond_sets(subgraphs)

    def combine_fragments(self, max_rotors=3, min_rotors=1, min_heavy_atoms=5, **kwargs):
        """

        Parameters
        ----------
        max_rotors :
        min_rotors :
        min_heavy_atoms :
        kwargs :

        Returns
        -------

        """
        #Only keep unique molecules
        unique_mols = set()
        nrfrags = len(self._fragments)
        for n in range(1, nrfrags):
            for fragcomb in combinations(self._fragments, n):
                if self._is_combination_adjacent(fragcomb):
                    frag = self._combine_and_connect_atom_bond_sets(fragcomb)
                    if (self._count_rotors_in_fragment(frag) <= max_rotors) and \
                            (self._count_rotors_in_fragment(frag) >= min_rotors):
                        if self._count_heavy_atoms_in_fragment(frag) >= min_heavy_atoms:
                            self._fragment_combinations.append(frag)
                            # convert to mol
                            mol = self._atom_bond_set_to_mol(frag, **kwargs)
                            smiles = []
                            if not isinstance(mol, list):
                                mol = [mol]
                            new_mols = []
                            for m in mol:
                                mapped_smiles = oechem.OEMolToSmiles(m)
                                if not mapped_smiles in unique_mols:
                                    # ToDo what to do when fragment loses stereo info?
                                    unique_mols.add(mapped_smiles)
                                    new_mols.append(m)
                                    try:
                                        smiles.append(mol_to_smiles(m, isomeric=True, mapped=False, explicit_hydrogen=False))
                                    except ValueError:
                                        logger().warning("Fragment lost stereo information or now has a new stereocenter")
                                        fixed_m = self._fix_stereo(m)
                                        if isinstance(fixed_m, list):
                                            logger().warning('Could not fix missing stereo. Adding all steroisomers')
                                            for new_stereo in fixed_m:
                                                smiles.append(mol_to_smiles(new_stereo, isomeric=True, mapped=False,
                                                                            explicit_hydrogen=False))
                                            self.new_stereo.append(oechem.OEMolToSmiles(m))
                                        else:
                                            smiles.append(mol_to_smiles(fixed_m, isomeric=True, mapped=False, explicit_hydrogen=False))
                            for sm in smiles:
                                if sm not in self.fragments:
                                    self.fragments[sm] = []
                                self.fragments[sm].extend(new_mols)

    def _is_combination_adjacent(self, frag_combination):
        """
        This function was taken from Openeye cookbook
        https://docs.eyesopen.com/toolkits/cookbook/python/cheminfo/enumfrags.html
        """

        parts = [0] * len(frag_combination)
        nrparts = 0

        for idx, frag in enumerate(frag_combination):
            if parts[idx] != 0:
                continue

            nrparts += 1
            parts[idx] = nrparts
            self._traverse_fragments(frag, frag_combination, parts, nrparts)

        return nrparts == 1

    def _traverse_fragments(self, acting_fragment, frag_combination, parts, nrparts):
        """
        This function was taken from openeye cookbook
        https://docs.eyesopen.com/toolkits/cookbook/python/cheminfo/enumfrags.html
        """
        for idx, frag in enumerate(frag_combination):
            if parts[idx] != 0:
                continue

            if not self._are_atom_bond_sets_adjacent(acting_fragment, frag):
                continue

            parts[idx] = nrparts
            self._traverse_fragments(frag, frag_combination, parts, nrparts)

    @staticmethod
    def _are_atom_bond_sets_adjacent(frag_a, frag_b):
        """

        Parameters
        ----------
        frag_a :
        frag_b :

        Returns
        -------

        """

        for atom_a in frag_a.GetAtoms():
            for atom_b in frag_b.GetAtoms():
                if atom_a.GetBond(atom_b) is not None:
                    return True
        return False

    @staticmethod
    def _combine_and_connect_atom_bond_sets(adjacent_frag_list):
        """
        This function was taken from OpeneEye Cookbook
        https://docs.eyesopen.com/toolkits/cookbook/python/cheminfo/enumfrags.html
        """

        # combine atom and bond sets

        combined = oechem.OEAtomBondSet()
        for frag in adjacent_frag_list:
            for atom in frag.GetAtoms():
                combined.AddAtom(atom)
            for bond in frag.GetBonds():
                combined.AddBond(bond)

        # add connecting bonds

        for atom_a in combined.GetAtoms():
            for atom_b in combined.GetAtoms():
                if atom_a.GetIdx() < atom_b.GetIdx():
                    continue

                bond = atom_a.GetBond(atom_b)
                if bond is None:
                    continue
                if combined.HasBond(bond):
                    continue

                combined.AddBond(bond)

        return combined

    def to_json(self):
        """
        Write out fragments to JSON with provenance
        Returns
        -------

        """
        json_dict = {}
        for frag_smiles in self.fragments:
            json_dict[frag_smiles] = {}
            identifiers = cmiles.get_molecule_ids(frag_smiles, strict=False)
            json_dict[frag_smiles]['cmiles_identifiers'] = identifiers
            json_dict[frag_smiles]['provenance'] = self.get_provenance()
            json_dict[frag_smiles]['provenance']['routine']['fragment_molecule']['map_to_parent'] = []
            for mol in self.fragments[frag_smiles]:
                json_dict[frag_smiles]['provenance']['routine']['fragment_molecule']['map_to_parent'].append(oechem.OEMolToSmiles(mol))

        return json_dict


    def depict_fragments(self, fname, line_width=0.75):
        """
        Generate PDF of all combinatorial fragments with individual fragments color coded
        Parameters
        ----------
        fname : str
            filename for output PDF
        line_width : float
            width of drawn molecule lines
        """

        itf = oechem.OEInterface()
        oedepict.OEConfigure2DMolDisplayOptions(itf)
        oedepict.OEConfigureReportOptions(itf)

        oedepict.OEPrepareDepiction(self.molecule)


        ropts = oedepict.OEReportOptions()
        oedepict.OESetupReportOptions(ropts, itf)
        ropts.SetFooterHeight(25.0)
        ropts.SetHeaderHeight(ropts.GetPageHeight() / 4.0)
        report = oedepict.OEReport(ropts)

        opts = oedepict.OE2DMolDisplayOptions()
        oedepict.OESetup2DMolDisplayOptions(opts, itf)
        cellwidth, cellheight = report.GetCellWidth(), report.GetCellHeight()
        opts.SetDimensions(cellwidth, cellheight, oedepict.OEScale_AutoScale)
        opts.SetTitleLocation(oedepict.OETitleLocation_Hidden)
        opts.SetAtomColorStyle(oedepict.OEAtomColorStyle_WhiteMonochrome)
        opts.SetAtomLabelFontScale(1.2)

        # add index to keep track of bonds
        idx_tag = 'fragment_idx'
        tag = oechem.OEGetTag(idx_tag)
        for f_idx, frag in enumerate(self._fragments):
            for atom in frag.GetAtoms():
                atom.SetData(tag, f_idx)

        # setup depiction style
        n_frags = len(self._fragments)
        colors = [c for c in oechem.OEGetLightColors()]
        if len(colors) < n_frags:
            n = n_frags - len(colors)
            colors.extend([c for c in oechem.OEGetColors(oechem.OEBlueTint, oechem.OERed, n)])

        atom_glyph = ColorAtomByFragmentIndex(colors, tag)
        fade_hightlight = oedepict.OEHighlightByColor(oechem.OEGrey, line_width)

        for frag in self._fragment_combinations:
            cell = report.NewCell()
            display = oedepict.OE2DMolDisplay(self.molecule, opts)

            frag_atoms = oechem.OEIsAtomMember(frag.GetAtoms())
            frag_bonds = oechem.OEIsBondMember(frag.GetBonds())

            not_fragment_atoms = oechem.OENotAtom(frag_atoms)
            not_fragment_bonds = oechem.OENotBond(frag_bonds)

            oedepict.OEAddHighlighting(display, fade_hightlight, not_fragment_atoms, not_fragment_bonds)
            oegrapheme.OEAddGlyph(display, atom_glyph, frag_atoms)

            oedepict.OERenderMolecule(cell, display)
        cellwidth, cellheight = report.GetHeaderWidth(), report.GetHeaderHeight()
        opts.SetDimensions(cellwidth, cellheight, oedepict.OEScale_AutoScale)
        opts.SetAtomColorStyle(oedepict.OEAtomColorStyle_WhiteMonochrome)
        disp = oedepict.OE2DMolDisplay(self.molecule, opts)
        oegrapheme.OEAddGlyph(disp, atom_glyph, oechem.OEIsTrueAtom())

        headerpen = oedepict.OEPen(oechem.OEWhite, oechem.OELightGrey, oedepict.OEFill_Off, 2.0)
        for header in report.GetHeaders():
            oedepict.OERenderMolecule(header, disp)
            oedepict.OEDrawBorder(header, headerpen)

        return oedepict.OEWriteReport(fname, report)


class WBOFragmenter(Fragmenter):
    """
    Fragment engine for fragmenting molecules using Wiberg Bond Order
    """
    def __init__(self, molecule, functional_groups=None, verbose=False):
        super().__init__(molecule)

        self.verbose = verbose

        if functional_groups is None:
            fn = resource_filename('fragmenter', os.path.join('data', 'fgroup_smarts.yml'))
            with open(fn, 'r') as f:
                functional_groups = yaml.safe_load(f)
        self._options['functional_groups'] = functional_groups
        self._tag_functional_groups(functional_groups)

        self.rotors_wbo = {}
        self.ring_systems = {}
        self.fragments = {} # Fragments from fragmentation scheme for each rotatable bond
        self._fragments = {} # Used internally for depiction


    def fragment(self,  threshold=0.01, keep_non_rotor_ring_substituents=True, **kwargs):
        """
        Fragment molecules using the Wiberg Bond Order as a surrogate

        Parameters
        ----------
        keep_non_rotor_ring_substituents: bool
            If True, will always keep all non rotor substiuents on ring. If False, will only add them
            if they are ortho to rotatable bond or if it's needed for WBO to be within the threshold
        **heuristic : str, optional, default 'path_length'
            The path fragmenter should take when fragment needs to be grown out. The other option is
            'wbo'
        **threshold : float, optional, default 0.01
            The threshold for the central bond WBO. If the fragment WBO is below this threshold, fragmenter
            will grow out the fragment one bond at a time via the path specified by the heuristic option

        """
        # Capture options used
        self._options['keep_non_rotor_ring_substituents'] = keep_non_rotor_ring_substituents
        if 'threshold' not in self._options:
            self._options['threshold'] = threshold

        # Add threshold as attribute because it is used in more than one function
        setattr(self, 'threshold', threshold)
        self._options.update(kwargs)
        # Calculate WBO for molecule
        self.calculate_wbo()
        self._get_rotor_wbo()
        # Find ring systems
        self._find_ring_systems(keep_non_rotor_ring_substituents=keep_non_rotor_ring_substituents)

        # Build fragments
        for bond in self.rotors_wbo:
            self.build_fragment(bond,  **kwargs)


    def calculate_wbo(self, fragment=None, **kwargs):
        """
        Calculate WBO
        Parameters
        ----------
        fragment : oechem.OEMol
            fragment to recalculate WBO. When fragment is None, fragmenter assumes it's the full molecule and saves the
            calculated values in self.molecule

        Returns
        -------
        fragment with WBOs

        """
        if not fragment:
            time1 = time.time()
            self.molecule = get_charges(self.molecule)
            time2 = time.time()
            if self.verbose:
                logger().info('WBO took {} seconds to calculate'.format(time2-time1))
        else:
            time1 = time.time()
            fragment = get_charges(fragment, **kwargs)
            time2 = time.time()
            if self.verbose:
                logger().info('WBO took {} seconds to calculate'.format(time2-time1))
        return fragment


    def _find_rotatable_bonds(self):
        #ToDo: Add option to build fragments around terminal torsions (-OH, -NH2, -CH3)
        """
        Using SMARTS instead of OpenEye's built in IsRotor function so double bonds are also captured
        This does not find terminal rotatable bonds such as -OH, -NH2 -CH3.


        Returns
        -------
        rotatable_bonds: list of tuples
            list of rotatable bonds map indices [(m1, m2),...]

        """
        rotatable_bonds = []
        smarts = '[!$(*#*)&!D1]-,=;!@[!$(*#*)&!D1]'
        # Suppress H to avoid finding terminal bonds
        copy_mol = oechem.OEMol(self.molecule)
        oechem.OESuppressHydrogens(copy_mol)
        oechem.OESuppressHydrogens(copy_mol)
        qmol = oechem.OEQMol()
        if not oechem.OEParseSmarts(qmol, smarts):
            raise RuntimeError('Cannot parse SMARTS {}'.format(smarts))
        ss = oechem.OESubSearch(qmol)
        oechem.OEPrepareSearch(copy_mol, ss)
        unique = True
        for match in ss.Match(copy_mol, unique):
            b = []
            for ma in match.GetAtoms():
                b.append(ma.target.GetMapIdx())
            rotatable_bonds.append(tuple(b))
        return rotatable_bonds

    def _get_rotor_wbo(self):
        """
        Cache WBO for each rotatable bond
        """

        if 'nrotor' not in self.molecule.GetData() and 'cputime' not in self.molecule.GetData():
            logger().info("WBO was not calculated for this molecule. Calculating WBO...")
            self.calculate_wbo()
        rotatable_bonds = self._find_rotatable_bonds()
        for bond in rotatable_bonds:
            b = self.get_bond(bond)
            self.rotors_wbo[bond] = b.GetData('WibergBondOrder')


    def _find_ring_systems(self, keep_non_rotor_ring_substituents=True):

        """
        This function tags ring atom and bonds with ringsystem index

        Parameters
        ----------
        mol: OpenEye OEMolGraph

        Returns
        -------
        tagged_rings: dict
            maps ringsystem index to ring atom and bond indices

        """
        nringsystems, parts = oechem.OEDetermineRingSystems(self.molecule)
        for ringidx in range(1, nringsystems +1):
            ringidx_atoms = set()
            for atom in self.molecule.GetAtoms():
                if parts[atom.GetIdx()] == ringidx:
                    ringidx_atoms.add(atom.GetMapIdx())
                    tag = oechem.OEGetTag('ringsystem')
                    atom.SetData(tag, ringidx)
            # Find bonds in ring and tag
            ringidx_bonds = set()
            for a_idx in ringidx_atoms:
                atom = self.molecule.GetAtom(oechem.OEHasMapIdx(a_idx))
                for bond in atom.GetBonds():
                    nbrAtom = bond.GetNbr(atom)
                    nbrIdx = nbrAtom.GetMapIdx()
                    if nbrIdx in ringidx_atoms and nbrIdx != a_idx:
                        ringidx_bonds.add((a_idx, nbrIdx))
                        tag = oechem.OEGetTag('ringsystem')
                        bond.SetData(tag, ringidx)
            non_rotor_atoms = set()
            non_rotor_bond = set()
            if keep_non_rotor_ring_substituents:
                for m in ringidx_atoms:
                    ring_atom = self.molecule.GetAtom(oechem.OEHasMapIdx(m))
                    for a in ring_atom.GetAtoms():
                        if a.GetMapIdx() not in ringidx_atoms and not a.IsHydrogen():
                            # Check bond
                            bond = self.molecule.GetBond(ring_atom, a)
                            if not bond.IsRotor():
                                # Add atom and bond to ring system
                                non_rotor_atoms.add(a.GetMapIdx())
                                non_rotor_bond.add((a.GetMapIdx(), m))
                                # Cehck if bond is in a functional group
                                if 'fgroup' in bond.GetData():
                                    # Grab all atoms and bonds
                                    fgroup = a.GetData('fgroup')
                                    non_rotor_atoms.update(self.functional_groups[fgroup][0])
                                    non_rotor_bond.update(self.functional_groups[fgroup][-1])
            ringidx_atoms.update(non_rotor_atoms)
            ringidx_bonds.update(non_rotor_bond)
            self.ring_systems[ringidx] = (ringidx_atoms, ringidx_bonds)

    def _find_ortho_substituent(self, ring_idx, rot_bond):
        """
        Find ring substituents that are ortho to the rotatable bond that fragment is being built around

        Parameters
        ----------
        ring_idx : int
            index of ring
        rot_bond : tuple of ints
            map indices of bond (m1, m2)

        Returns
        -------
        ortho_atoms: set of ints
            set of map indices of atoms in ortho group
        ortho_bonds: set of tuples of ints
            set of bond tuples in ortho group

        """
        # Get the ring atom
        ring_atom = None
        for m in rot_bond:
            a = self.molecule.GetAtom(oechem.OEHasMapIdx(m))
            if a.IsInRing() and m in self.ring_systems[ring_idx][0]:
                ring_atom = a
        if not ring_atom:
            # rotatable bond is not directly bonded to this ring
            return

        # Get all atoms in the ring
        ring_atoms = [self.molecule.GetAtom(oechem.OEHasMapIdx(i)) for i in self.ring_systems[ring_idx][0]]
        ortho_atoms = set()
        ortho_bonds = set()
        for atom in ring_atoms:
            for a in atom.GetAtoms():
                if not a.IsHydrogen() and not a.GetMapIdx() in self.ring_systems[ring_idx][0]:
                    # Check if atom is bonded to ring atom of rotatable bond
                    if self.molecule.GetBond(atom, ring_atom):
                        # This is an ortho group
                        ortho_atoms.add(a.GetMapIdx())
                        ortho_bonds.add((a.GetMapIdx(), atom.GetMapIdx()))
                        # Check if substituent is part of functional group
                        if 'fgroup' in a.GetData():
                            fgroup = a.GetData('fgroup')
                            ortho_atoms.update(self.functional_groups[fgroup][0])
                            ortho_bonds.update(self.functional_groups[fgroup][-1])
                        # Check if substituent is a ring
                        if 'ringsystem' in a.GetData() and a.GetData('ringsystem') != ring_idx:
                            ring_system = a.GetData('ringsystem')
                            ortho_atoms.update(self.ring_systems[ring_system][0])
                            ortho_bonds.update(self.ring_systems[ring_system][-1])
        return ortho_atoms, ortho_bonds


    def compare_wbo(self, fragment, bond_tuple, **kwargs):
        """
        Compare Wiberg Bond order of rotatable bond in fragment to parent

        Parameters
        ----------
        fragment :
        bond_tuple :
        kwargs :

        Returns
        -------

        """

        restore_atom_map(fragment)
        try:
            charged_fragment = self.calculate_wbo(fragment=fragment, normalize=False,  **kwargs)
        except RuntimeError:
            raise RuntimeError("Cannot calculate WBO for fragment built around bond {}".format(bond_tuple))
        # Get new WBO
        a1 = charged_fragment.GetAtom(oechem.OEHasMapIdx(bond_tuple[0]))
        a2 = charged_fragment.GetAtom(oechem.OEHasMapIdx(bond_tuple[-1]))
        bond = charged_fragment.GetBond(a1, a2)
        if bond is None:
            raise RuntimeError('{} with _idx {} and {} with map_idx {} are not bonded'.format(a1,
                                                                                              bond_tuple[0],
                                                                                              a2,
                                                                                     bond_tuple[1]))
        frag_wbo = bond.GetData('WibergBondOrder')
        mol_wbo = self.rotors_wbo[bond_tuple]
        self.fragments[bond_tuple] = charged_fragment
        return abs(mol_wbo - frag_wbo)

    def get_bond(self, bond_tuple):
        """
        Get bond in molecule by atom indices of atom A and atom B

        Parameters
        ----------
        bond_tuple : tuple
            (mapidx, mapidx)

        Returns
        -------
        oechem.OEBondBase

        """
        if is_missing_atom_map(self.molecule):
            restore_atom_map(self.molecule)
        a1 = self.molecule.GetAtom(oechem.OEHasMapIdx(bond_tuple[0]))
        a2 = self.molecule.GetAtom(oechem.OEHasMapIdx(bond_tuple[-1]))
        bond = self.molecule.GetBond(a1, a2)
        if not bond:
            raise ValueError("({}) atoms are not connected".format(bond_tuple))
        return bond

    def build_fragment(self, bond_tuple, heuristic='path_length', **kwargs):
        """
        Build fragment around bond
        Parameters
        ----------
        bond_tuple : tuple
            tuple of atom maps of atoms in bond

        Returns
        -------

        """
        # Capture options

        if 'heuristic' not in self._options:
            self._options['heuristic'] = heuristic
        bond = self.get_bond(bond_tuple=bond_tuple)
        atom_map_idx = set()
        bond_tuples = set()
        #bond_idx = set()
        a1, a2 = bond.GetBgn(), bond.GetEnd()
        m1, m2 = a1.GetMapIdx(), a2.GetMapIdx()
        atom_map_idx.update((m1, m2))
        #bond = self.molecule.GetBond(a1, a2)
        #bond_idx.add(bond.GetIdx())
        bond_tuples.add((m1, m2))
        for atom in (a1, a2):
            for a in atom.GetAtoms():
                m = a.GetMapIdx()
                atom_map_idx.add(m)
                bond_tuples.add((atom.GetMapIdx(), m))
                if 'ringsystem' in a.GetData():
                    # Grab the ring
                    ring_idx = a.GetData('ringsystem')
                    atom_map_idx.update(self.ring_systems[ring_idx][0])
                    bond_tuples.update(self.ring_systems[ring_idx][-1])
                    #bond_idx.update(self.ring_systems[ring_idx][-1])
                    # Check for ortho substituents here
                    ortho = self._find_ortho_substituent(ring_idx=ring_idx, rot_bond=bond_tuple)
                    if ortho:
                        atom_map_idx.update(ortho[0])
                        bond_tuples.update(ortho[1])
                if 'fgroup' in a.GetData():
                    # Grab rest of fgroup
                    fgroup = a.GetData('fgroup')
                    atom_map_idx.update(self.functional_groups[fgroup][0])
                    bond_tuples.update(self.functional_groups[fgroup][-1])

        atom_bond_set = self._to_atom_bond_set(atom_map_idx, bond_tuples)
        fragment_mol = self._atom_bond_set_to_mol(atom_bond_set)
        # #return fragment_mol
        diff = self.compare_wbo(fragment_mol, bond_tuple, **kwargs)
        if diff <= self.threshold:
            self._fragments[bond_tuple] = atom_bond_set
        if diff > self.threshold:
            self._add_next_substituent(atom_map_idx, bond_tuples, target_bond=bond_tuple,
                                       heuristic=heuristic, **kwargs)


    def _add_next_substituent(self, atoms, bonds, target_bond, heuristic='path_length', **kwargs):
        """
        If the difference between WBO in fragment and molecules is greater than threshold, add substituents to
        fragment until difference is within threshold
        Parameters
        ----------
        atoms :
        bonds :
        target_bond :
        threshold :
        heuristic: str
            How to add substituents. Choices are path_length or wbo

        Returns
        -------

        """
        bond_atom_1 = self.molecule.GetAtom(oechem.OEHasMapIdx(target_bond[0]))
        bond_atom_2 = self.molecule.GetAtom(oechem.OEHasMapIdx(target_bond[1]))
        atoms_to_add = []
        sort_by_1 = []
        sort_by_2 = []
        sort_by = []
        for m_idx in atoms:
            a = self.molecule.GetAtom(oechem.OEHasMapIdx(m_idx))
            for nbr in a.GetAtoms():
                nbr_m_idx = nbr.GetMapIdx()
                if not nbr.IsHydrogen() and not nbr_m_idx in atoms:
                    b = self.molecule.GetBond(a, nbr)
                    atoms_to_add.append((nbr_m_idx, (m_idx, nbr_m_idx)))
                    if heuristic == 'wbo':
                        sort_by.append(b.GetData('WibergBondOrder'))
                        reverse = True
                    elif heuristic == 'path_length':
                        sort_by_1.append(oechem.OEGetPathLength(bond_atom_1, nbr))
                        sort_by_2.append(oechem.OEGetPathLength(bond_atom_2, nbr))
                        reverse = False
                    else:
                        raise ValueError('Only wbo and path_lenght are supported heuristics')
        if heuristic == 'path_length':
            min_1 = min(sort_by_1)
            min_2 = min(sort_by_2)
            if min_1 < min_2:
                sort_by = sort_by_1
            elif min_2 < min_1:
                sort_by = sort_by_2
            elif min_1 == min_2:
                indices = []
                for l in [sort_by_1, sort_by_2]:
                    indices.extend([i for i, x in enumerate(l) if x == min_1])
                atoms_to_add = [atoms_to_add[i] for i in indices]
                for a, b in atoms_to_add:
                    bond = self.get_bond(b)
                    sort_by.append(bond.GetData('WibergBondOrder'))
                    reverse = True

        sorted_atoms = [a for _, a in sorted(zip(sort_by, atoms_to_add), reverse=reverse)]
        for atom, bond in sorted_atoms:
            a = self.molecule.GetAtom(oechem.OEHasMapIdx(atom))
            if 'ringsystem' in a.GetData():
                ring_idx = a.GetData('ringsystem')
                atoms.update(self.ring_systems[ring_idx][0])
                bonds.update(self.ring_systems[ring_idx][1])
            if 'fgroup' in a.GetData():
                fgroup = a.GetData('fgroup')
                atoms.update(self.functional_groups[fgroup][0])
                bonds.update(self.functional_groups[fgroup][1])
            atoms.add(atom)
            bonds.add(bond)
            # Check new WBO
            atom_bond_set = self._to_atom_bond_set(atoms, bonds)
            fragment_mol = self._atom_bond_set_to_mol(atom_bond_set)
            if self.compare_wbo(fragment_mol, target_bond, **kwargs) < self.threshold:
                self._fragments[target_bond] = atom_bond_set
                return fragment_mol
            else:
                return self._add_next_substituent(atoms=atoms, bonds=bonds, target_bond=target_bond,
                                                  heuristic=heuristic, **kwargs)

    def to_json(self):
        """
        Write out fragments to JSON with provenance
        Returns
        -------

        """
        json_dict = {}
        for bond in self.fragments:
            can_iso_smiles = cmiles.utils.mol_to_smiles(self.fragments[bond], mapped=False)
            identifiers = cmiles.get_molecule_ids(can_iso_smiles)
            can_iso_smiles = identifiers['canonical_isomeric_smiles']
            json_dict[can_iso_smiles] = {'cmiles_identifiers': identifiers}
            json_dict[can_iso_smiles]['provenance'] = self.get_provenance()
        return json_dict

    def _to_qcschema_mol(self, molecule, **kwargs):
        """

        Parameters
        ----------
        molecule :
        kwargs :

        Returns
        -------

        """
        self._options.update(kwargs)
        mol_copy = oechem.OEMol(molecule)
        oechem.OEAddExplicitHydrogens(mol_copy)
        explicit_h_smiles = cmiles.utils.mol_to_smiles(mol_copy, mapped=False)
        cmiles_identifiers = cmiles.get_molecule_ids(explicit_h_smiles)
        can_mapped_smiles = cmiles_identifiers['canonical_isomeric_explicit_hydrogen_mapped_smiles']

        conformers = generate_conformers(mol_copy, **kwargs)
        qcschema_mols = [cmiles.utils.mol_to_map_ordered_qcschema(conf, can_mapped_smiles) for conf in conformers.GetConfs()]

        return {'initial_molecule': qcschema_mols,
                'identifiers': cmiles_identifiers,
                'provenance': self.get_provenance()}

    def to_qcschema_mols(self, **kwargs):
        """

        Parameters
        ----------
        kwargs :

        Returns
        -------

        """
        qcschema_fragments = []
        equivelant_frags = []
        for bond in self.fragments:
            molecule = self.fragments[bond]
            qcschema_mol = self._to_qcschema_mol(molecule, **kwargs)
            smiles = qcschema_mol['identifiers']['canonical_isomeric_smiles']
            if smiles not in equivelant_frags:
                equivelant_frags.append(smiles)
                qcschema_fragments.append(qcschema_mol)
        return qcschema_fragments


    def to_torsiondrive_json(self, **kwargs):
        """

        Returns
        -------

        """
        # capture options
        self._options.update(kwargs)
        torsiondrive_json_dict = {}
        for bond in self.fragments:
            # find torsion around bond in fragment
            molecule = self.fragments[bond]
            mapped_smiles = oechem.OEMolToSmiles(molecule)
            torsion_map_idx = find_torsion_around_bond(molecule, bond)
            torsiondrive_job_label = cmiles.utils.to_canonical_label(mapped_smiles, torsion_map_idx)
            torsiondrive_json_dict[torsiondrive_job_label] = {}

            # prepare torsiondrive input
            mol_copy = oechem.OEMol(molecule)

            # Map torsion to canonical ordered molecule
            # First canonical order the molecule
            cmiles._cmiles_oe.canonical_order_atoms(mol_copy)
            dih = []
            for m_idx in torsion_map_idx:
                atom = mol_copy.GetAtom(oechem.OEHasMapIdx(m_idx + 1))
                dih.append(atom.GetIdx())

            torsiondrive_json_dict[torsiondrive_job_label] = self._to_qcschema_mol(mol_copy, **kwargs)
            torsiondrive_json_dict[torsiondrive_job_label].update({'dihedral': [dih],
                                                              'grid': [15],
                                                              'provenance': self.get_provenance()})
        return torsiondrive_json_dict


    def depict_fragments(self, fname):
        """

        Parameters
        ----------
        fname :

        Returns
        -------

        """
        itf = oechem.OEInterface()
        oedepict.OEPrepareDepiction(self.molecule)

        ropts = oedepict.OEReportOptions()
        oedepict.OESetupReportOptions(ropts, itf)
        ropts.SetFooterHeight(25.0)
        ropts.SetHeaderHeight(ropts.GetPageHeight() / 4.0)
        report = oedepict.OEReport(ropts)

        # setup decpiction options
        opts = oedepict.OE2DMolDisplayOptions()
        oedepict.OESetup2DMolDisplayOptions(opts, itf)
        cellwidth, cellheight = report.GetCellWidth(), report.GetCellHeight()
        opts.SetDimensions(cellwidth, cellheight, oedepict.OEScale_AutoScale)
        opts.SetTitleLocation(oedepict.OETitleLocation_Hidden)
        opts.SetAtomColorStyle(oedepict.OEAtomColorStyle_WhiteMonochrome)
        opts.SetAtomLabelFontScale(1.2)

        lineWidthScale = 0.75
        fadehighlight = oedepict.OEHighlightByColor(oechem.OEGrey, lineWidthScale)

        # depict each fragment combinations

        for bond_tuple in self._fragments:

            cell = report.NewCell()
            bond = self.get_bond(bond_tuple)
            # Get bond in fragment
            fragment = self.fragments[bond_tuple]
            a1 = fragment.GetAtom(oechem.OEHasMapIdx(bond_tuple[0]))
            a2 = fragment.GetAtom(oechem.OEHasMapIdx(bond_tuple[1]))
            bond_in_frag = fragment.GetBond(a1, a2)
            wbo_frag = bond_in_frag.GetData('WibergBondOrder')
            bond.SetData('WibergBondOrder_frag', wbo_frag)

            bondlabel_2 = LabelFragBondOrder()
            opts.SetBondPropertyFunctor(bondlabel_2)
            disp = oedepict.OE2DMolDisplay(self.molecule, opts)

            bond.DeleteData('WibergBondOrder_frag')

            # ToDo get AtomBondSet for fragments so depiction works properly
            fragatoms = oechem.OEIsAtomMember(self._fragments[bond_tuple].GetAtoms())
            fragbonds = oechem.OEIsBondMember(self._fragments[bond_tuple].GetBonds())

            notfragatoms = oechem.OENotAtom(fragatoms)
            notfragbonds = oechem.OENotBond(fragbonds)

            oedepict.OEAddHighlighting(disp, fadehighlight, notfragatoms, notfragbonds)

            atomBondSet = oechem.OEAtomBondSet()
            atomBondSet.AddBond(bond)
            atomBondSet.AddAtom(bond.GetBgn())
            atomBondSet.AddAtom(bond.GetEnd())

            hstyle = oedepict.OEHighlightStyle_BallAndStick
            hcolor = oechem.OEColor(oechem.OELightBlue)
            oedepict.OEAddHighlighting(disp, hcolor, hstyle, atomBondSet)

            #oegrapheme.OEAddGlyph(disp, bondglyph, fragbonds)

            oedepict.OERenderMolecule(cell, disp)

        # depict original fragmentation in each header

        cellwidth, cellheight = report.GetHeaderWidth(), report.GetHeaderHeight()
        opts.SetDimensions(cellwidth, cellheight, oedepict.OEScale_AutoScale)
        opts.SetAtomColorStyle(oedepict.OEAtomColorStyle_WhiteMonochrome)

        bondlabel = LabelWibergBondOrder()
        opts.SetBondPropertyFunctor(bondlabel)
        disp = oedepict.OE2DMolDisplay(self.molecule, opts)
        #oegrapheme.OEAddGlyph(disp, bondglyph, oechem.IsTrueBond())

        headerpen = oedepict.OEPen(oechem.OEWhite, oechem.OELightGrey, oedepict.OEFill_Off, 2.0)
        for header in report.GetHeaders():
            oedepict.OERenderMolecule(header, disp)
            oedepict.OEDrawBorder(header, headerpen)

        return oedepict.OEWriteReport(fname, report)
