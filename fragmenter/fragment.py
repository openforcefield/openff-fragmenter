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
from .chemi import to_smi, get_charges, generate_conformers, LabelWibergBondOrder
from .torsions import find_torsion_around_bond


OPENEYE_VERSION = oe.__name__ + '-v' + oe.__version__


def expand_states(inp_molecule, protonation=True, tautomers=False, stereoisomers=True, max_states=200, level=0, reasonable=True,
                  carbon_hybridization=True, suppress_hydrogen=True, verbose=False, filename=None,
                  return_smiles_list=False, return_molecules=False, expand_internally=True, strict=False, mapped=False):
    """
    Expand molecule states (choice of protonation, tautomers and/or stereoisomers).
    Protonation states expands molecules to protonation of protonation sites (Some states might only be reasonable in
    very high or low pH. ToDo: Only keep reasonable protonation states)
    Tatutomers: Should expand to tautomer states but most of hte results are some resonance structures. Defualt if False
    for this reason
    Stereoisomers expands enantiomers and geometric isomers (cis/trans).
    Returns set of SMILES

    Parameters
    ----------
    molecule: OEMol
        Molecule to expand
    protonation: Bool, optional, default=True
        If True will enumerate protonation states.
    tautomers: Bool, optional, default=False
        If True, will enumerate tautomers.  (Note: Default is False because results usually give resonance structures
        which ins't needed for torsion scans
    stereoisomers: Bool, optional, default=True
        If True will enumerate stereoisomers (cis/trans and R/S).
    max_states: int, optional, default=True
        maximum states enumeration should find
    level: int, optional, Defualt=0
        The level for enumerating tautomers. It can go up until 7. The higher the level, the more tautomers will be
        generated but they will also be less reasonable.
    reasonable: bool, optional, default=True
        Will rank tautomers enumerated energetically (https://docs.eyesopen.com/toolkits/python/quacpactk/tautomerstheory.html#reasonable-ranking)
    carbon_hybridization: bool, optional, default=True
        If True will allow carbons to change hybridization
    suppress_hydrogen: bool, optional, default=True
        If true, will suppress explicit hydrogen. It's considered best practice to set this to True when enumerating tautomers.
    verbose: Bool, optional, default=True
    filename: str, optional, default=None
        Filename to save SMILES to. If None, SMILES will not be saved to file.
    return_smiles_list: bool, optional, default=False
        If True, will return a list of SMILES with numbered name of molecule. Use this if you want ot write out an
        smi file of all molecules processed with a unique numbered name for each state.
    return_molecules: bool, optional, default=False
        If true, will return list of OEMolecules instead of SMILES

    Returns
    -------
    states: set of SMILES for enumerated states

    """
    #ToDo Stereoisomers should probably go first to avoid adding molecules that do not have stereo defined
    #ToDo Separate each expansion to its own function
    title = inp_molecule.GetTitle()
    states = set()
    molecules = [inp_molecule]
    if verbose:
        logger().info("Enumerating states for {}".format(title))
    if protonation:
        logger().info("Enumerating protonation states for {}".format(title))
        molecules.extend(_expand_states(molecules, enumerate='protonation', max_states=max_states, verbose=verbose,
                                        level=level, suppress_hydrogen=suppress_hydrogen))
    if tautomers:
        logger().info("Enumerating tautomers for {}".format(title))
        molecules.extend(_expand_states(molecules, enumerate='tautomers', max_states=max_states, reasonable=reasonable,
                                        carbon_hybridization=carbon_hybridization, verbose=verbose, level=level,
                                        suppress_hydrogen=suppress_hydrogen))
    if stereoisomers:
        logger().info("Enumerating stereoisomers for {}".format(title))
        molecules.extend(_expand_states(molecules, enumerate='stereoisomers', max_states=max_states, verbose=verbose,
                                        suppress_hydrogen=suppress_hydrogen))

    # if stereoisomers:
    #     molecules.remove(inp_molecule)
    for molecule in molecules:
        #states.add(fragmenter.utils.create_mapped_smiles(molecule, tagged=False, explicit_hydrogen=False))
        # Not using create mapped SMILES because OEMol is needed but state is OEMolBase.
        #states.add(oechem.OEMolToSmiles(molecule))
        try:
            states.add(mol_to_smiles(molecule, isomeric=True, mapped=mapped, explicit_hydrogen=True))
        except ValueError:
            if stereoisomers:
                continue

            if expand_internally:
                logger().warn("Tautomer or protonation state has a chiral center. Expanding stereoisomers")
                max_stereo_states = int(max_states/len(molecules))
                stereo_states = _expand_states(molecule, enumerate='stereoisomers', max_states=max_stereo_states)
                for state in stereo_states:
                    if mapped:
                        try:
                            states.add(mol_to_smiles(state, isomeric=True, mapped=mapped, explicit_hydrogen=True))
                        except ValueError:
                            logger().warning('Cannot generate isomeric SMILES for {}. Skipping'.format(oechem.OEMolToSmiles(state)))
                            continue
                    else:
                        states.add(mol_to_smiles(state, isomeric=True, mapped=mapped, explicit_hydrogen=False))
            elif not strict:
                states.add(mol_to_smiles(molecule, isomeric=False, mapped=mapped, explicit_hydrogen=True ))

            else:
                raise ValueError("molecule {} is missing stereochemistry".format(mol_to_smiles(molecule, isomeric=False, mapped=False, explicit_hydorge=False)))


    logger().info("{} states were generated for {}".format(len(states), oechem.OEMolToSmiles(molecule)))

    if filename:
        count = 0
        smiles_list = []
        for molecule in states:
            molecule = molecule + ' ' + title + '_' + str(count)
            count += 1
            smiles_list.append(molecule)
        to_smi(smiles_list, filename)

    if return_smiles_list:
        return smiles_list

    if return_molecules:
        return molecules

    return states


def _expand_states(molecules, enumerate='protonation', max_states=200, suppress_hydrogen=True, reasonable=True,
                   carbon_hybridization=True, level=0, verbose=True):
    """
    Expand the state specified by enumerate variable

    Parameters
    ----------
    molecules: OEMol or list of OEMol
        molecule to expand states
    enumerate: str, optional, default='protonation'
        Kind of state to enumerate. Choice of protonation, tautomers, stereoiserms
    suppress_hydrogen: bool, optional, default=True
        If True, will suppress explicit hydrogen
    reasonable: bool, optional, default=True
        If True, will rank tautomers by the most reasonable energetically
    carbon_hybridization: bool, optional, default=True
        If True, will allow carbon to change hybridization
    max_states: int, optional, default=200
    verbose: Bool, optional, deault=TRue

    Returns
    -------
    states: list of OEMol
        enumerated states

    """
    if type(molecules) != type(list()):
        molecules = [molecules]

    states = list()
    for molecule in molecules:
        ostream = oechem.oemolostream()
        ostream.openstring()
        ostream.SetFormat(oechem.OEFormat_SDF)
        states_enumerated = 0
        if suppress_hydrogen:
            oechem.OESuppressHydrogens(molecule)
        if enumerate == 'protonation':
            formal_charge_options = oequacpac.OEFormalChargeOptions()
            formal_charge_options.SetMaxCount(max_states)
            formal_charge_options.SetVerbose(verbose)
            #if verbose:
            #    logger().debug("Enumerating protonation states...")
            for protonation_state in oequacpac.OEEnumerateFormalCharges(molecule, formal_charge_options):
                states_enumerated += 1
                oechem.OEWriteMolecule(ostream, protonation_state)
                states.append(protonation_state)
        if enumerate == 'tautomers':
            #max_zone_size = molecule.GetMaxAtomIdx()
            tautomer_options = oequacpac.OETautomerOptions()
            tautomer_options.SetMaxTautomersGenerated(max_states)
            tautomer_options.SetLevel(level)
            tautomer_options.SetRankTautomers(reasonable)
            tautomer_options.SetCarbonHybridization(carbon_hybridization)
            #tautomer_options.SetMaxZoneSize(max_zone_size)
            tautomer_options.SetApplyWarts(True)
            if verbose:
                logger().debug("Enumerating tautomers...")
            for tautomer in oequacpac.OEEnumerateTautomers(molecule, tautomer_options):
                states_enumerated += 1
                states.append(tautomer)
        if enumerate == 'stereoisomers':
            if verbose:
                logger().debug("Enumerating stereoisomers...")
            for enantiomer in oeomega.OEFlipper(molecule, max_states, True):
                states_enumerated += 1
                enantiomer = oechem.OEMol(enantiomer)
                oechem.OEWriteMolecule(ostream, enantiomer)
                states.append(enantiomer)

    return states


class Fragmenter(object):

    def __init__(self, molecule):
        self.molecule = molecule
        self.functional_groups = {}

         # Add explicit hydrogen
        if not cmiles.utils.has_explicit_hydrogen(self.molecule):
            oechem.OEAddExplicitHydrogens(self.molecule)
        # Check for atom maps - it is used for finding atoms to tag
        if not has_atom_map(self.molecule):
            cmiles.utils.add_atom_map(self.molecule)

        # Keep track of stereo to make sure it does not flip
        self.stereo = {}
        self._atom_stereo_map = {oechem.OECIPAtomStereo_S: 'S', oechem.OECIPAtomStereo_R: 'R'}
        self._bond_stereo_map = {oechem.OECIPBondStereo_E: 'E', oechem.OECIPBondStereo_Z: 'Z'}
        self._find_stereo()

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
        for a in self.molecule.GetAtoms():
            if a.IsChiral():
                self.stereo[a.GetMapIdx()] = self._atom_stereo_map[oechem.OEPerceiveCIPStereo(self.molecule, a)]
        for bond in self.molecule.GetBonds():
            if bond.IsChiral():
                m1 = bond.GetBgn().GetMapIdx()
                m2 = bond.GetEnd().GetMapIdx()
                self.stereo[(m1, m2)] = self._bond_stereo_map[oechem.OEPerceiveCIPStereo(self.molecule, bond)]

    def _check_stereo(self, fragment, use_parent=False):
        for a in fragment.GetAtoms():
            if a.IsChiral():
                if use_parent:
                    s = self._atom_stereo_map[oechem.OEPerceiveCIPStereo(self.molecule, a)]
                else:
                    s = self._atom_stereo_map[oechem.OEPerceiveCIPStereo(fragment, a)]
                if not s == self.stereo[a.GetMapIdx()]:
                    raise AssertionError('Stereochemistry for atom {} flipped from {} to {}'.format(a.GetMapIdx(), self.stereo[a.GetMapIdx()],
                                                                                            s))
        for b in fragment.GetBonds():
            if b.IsChiral():
                b_tuple = (b.GetBgh().GetMapIdx(), b.GetEnd().GetMapIdx())
                if use_parent:
                    s = oechem.OEPerceiveCIPStereo(self.molecule, b)
                else:
                    s = oechem.OEPerceiveCIPStereo(fragment, b)
                if not s == self._bond_stereo_map[self.stereo[b_tuple]]:
                    raise AssertionError('Stereochemistry fro bond {} flipped from {} to {}'.format(b_tuple, self.stereo[b_tuple],
                                                                                                    self._bond_stereo_map[s]))

    def _tag_functional_groups(self, functional_groups):
        """
        This function tags atoms and bonds of functional groups defined in fgroup_smarts. fgroup_smarts is a dictionary
        that maps functional groups to their smarts pattern. It can be user generated or from yaml file.

        Parameters
        ----------
        mol: Openeye OEMolGraph
        frgroups_smarts: dictionary of functional groups mapped to their smarts pattern.
            Default is None. It uses 'fgroup_smarts.yaml'

        Returns
        -------
        fgroup_tagged: dict
            a dictionary that maps indexed functional groups to corresponding atom and bond indices in mol

        """
        # Add explicit hydrogen
        if not cmiles.utils.has_explicit_hydrogen(self.molecule):
            oechem.OEAddExplicitHydrogens(self.molecule)
        # Check for atom maps - it is used for finding atoms to tag
        if not has_atom_map(self.molecule):
            cmiles.utils.add_atom_map(self.molecule)
        if is_missing_atom_map(self.molecule):
            warnings.warn("Some atoms are missing atom maps. This might cause a problem at several points during "
                          "the fragmentation process. Make sure you know what you are doing. ")
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

    def atom_bond_set_to_mol(self, frag, adjust_hcount=True, expand_stereoisomers=True):
        """
        Convert fragments (AtomBondSet) to OEMol
        Parameters
        ----------
        frag: OEAtomBondSet
        mol: OEMol
        OESMILESFlag: str
            Either 'ISOMERIC' or 'DEFAULT'. This flag determines which OE function to use to generate SMILES string

        Returns
        -------
        smiles: dict of smiles to frag

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
        if has_stereo_defined(fragment):
            self._check_stereo(fragment)

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


class CombinatorialFragmenter(Fragmenter):
    """

    """
    def __init__(self, molecule, functional_groups=None):
        super().__init__(molecule)

        if functional_groups is None:
            fn = resource_filename('fragmenter', os.path.join('data', 'fgroup_smarts_comb.yml'))
            with open(fn, 'r') as f:
                functional_groups = yaml.safe_load(f)
        self._tag_functional_groups(functional_groups)
        self._nx_graph = self._mol_to_graph()
        self._fragments = list()  # all possible fragments without breaking rings
        self._fragment_combinations = list() # AtomBondSets of combined fragments. Used internally to generate PDF
        self.fragments = {} # Dict that maps SMILES to all equal combined fragments

    def fragment(self, min_rotors=1, max_rotors=None, min_heavy_atoms=4):
        """

        Returns
        -------

        """
        # Remove atom map and store it in data. This is needed to keep track of the fragments
        if has_atom_map(self.molecule):
            remove_atom_map(self.molecule, keep_map_data=True)
        self.fragment_all_bonds_not_in_ring_systems()
        if max_rotors is None:
            max_rotors = self.n_rotors + 1
        self.combine_fragments(min_rotors=min_rotors, max_rotors=max_rotors, min_heavy_atoms=min_heavy_atoms)

    def _mol_to_graph(self):

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
        Fragment all bonds that are not in rings
        Parameters
        ----------
        G: NetworkX graph
        bondOrderThreshold: int
            thershold for fragmenting graph. Default 1.2

        Returns
        -------
        subgraphs: list of subgraphs
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
        graph: NetworkX graph
        subgraph: NetworkX subgraph
        oemol: Openeye OEMolGraph

        Returns
        ------
        atomBondSet: Openeye oechem atomBondSet
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
        subgraphs = self._fragment_graph()
        self._subgraphs_to_atom_bond_sets(subgraphs)

    def combine_fragments(self, max_rotors=3, min_rotors=1, min_heavy_atoms=5, **kwargs):
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
                            mol = self.atom_bond_set_to_mol(frag, **kwargs)
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
                                        smiles.append(mol_to_smiles(m, isomeric=True, mapped=False, explicit_hydrogen=True))
                                    except ValueError:
                                        # fragment lost stereo information.
                                        logger().warning("Fragment lost stereo information or now has a new sterecenter. Expanding to 2 stereoisomers")
                                        states = _expand_states(m, enumerate='stereoisoemrs', max_states=2)
                                        smiles.extend(list(states))
                                        #smiles.append(mol_to_smiles(m, isomeric=False, mapped=False, explicit_hydrogen=True))
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

    def depict_fragments(self, fname, line_width=0.75):

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

    """
    def __init__(self, molecule, functional_groups=None, verbose=False):
        super().__init__(molecule)

        self.verbose = verbose

        if functional_groups is None:
            fn = resource_filename('fragmenter', os.path.join('data', 'fgroup_smarts.yml'))
            with open(fn, 'r') as f:
                functional_groups = yaml.safe_load(f)
        self._tag_functional_groups(functional_groups)

        self.rotors_wbo = {}
        self.ring_systems = {}
        self.fragments = {} # Fragments from fragmentation scheme for each rotatable bond
        self._fragments = {} # Used internally for depiction

        # For provenance
        self._options = {'functional_groups': functional_groups}

    def fragment(self,  keep_non_rotor_ring_substituents=True, **kwargs):
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

        Parameters
        ----------
        ring_idx :
        rot_bond :

        Returns
        -------

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

    def build_fragment(self, bond_tuple, threshold=0.01, heuristic='path_length', **kwargs):
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
        if 'threshold' not in self._options:
            self._options['threshold'] = threshold
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
        fragment_mol = self.atom_bond_set_to_mol(atom_bond_set, expand_stereoisomers=False)
        # #return fragment_mol
        diff = self.compare_wbo(fragment_mol, bond_tuple, **kwargs)
        if diff <= threshold:
            self._fragments[bond_tuple] = atom_bond_set
        if diff > threshold:
            self._add_next_substituent(atom_map_idx, bond_tuples, target_bond=bond_tuple,
                                                      threshold=threshold, heuristic=heuristic, **kwargs)


    def _add_next_substituent(self, atoms, bonds, target_bond, threshold=0.01, heuristic='path_length', **kwargs):
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
            fragment_mol = self.atom_bond_set_to_mol(atom_bond_set, expand_stereoisomers=False)
            if self.compare_wbo(fragment_mol, target_bond, **kwargs) < threshold:
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
        for bond in self.fragments:
            molecule = self.fragments[bond]
            qcschema_fragments.append(self._to_qcschema_mol(molecule, **kwargs))

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
        provenance = {'creator': fragmenter.__package__,
                      'job_id': str(uuid.uuid4()),
                      'hostname': socket.gethostname(),
                      'username': getpass.getuser(),
                      'routine': {'fragment_molecule': {
                          'version': fragmenter_version,
                          'options': self._options,
                          'parent_molecule': cmiles.utils.mol_to_smiles(self.molecule, mapped=False,
                                                                        explicit_hydrogen=False)
                      }}}
        return provenance


    def depict_fragments(self, fname):
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

            bondlabel_2 = LabelBondOrder()
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


class ColorAtomByFragmentIndex(oegrapheme.OEAtomGlyphBase):
    """
    This class was taken from OpeneEye cookbook
    https://docs.eyesopen.com/toolkits/cookbook/python/depiction/enumfrags.html
    """
    def __init__(self, colorlist, tag):
        oegrapheme.OEAtomGlyphBase.__init__(self)
        self.colorlist = colorlist
        self.tag = tag

    def RenderGlyph(self, disp, atom):

        a_disp = disp.GetAtomDisplay(atom)
        if a_disp is None or not a_disp.IsVisible():
            return False

        if not atom.HasData(self.tag):
            return False

        linewidth = disp.GetScale() / 1.5
        color = self.colorlist[atom.GetData(self.tag)]
        radius = disp.GetScale() / 4.8
        pen = oedepict.OEPen(color, color, oedepict.OEFill_Off, linewidth)

        layer = disp.GetLayer(oedepict.OELayerPosition_Below)
        oegrapheme.OEDrawCircle(layer, oegrapheme.OECircleStyle_Default, a_disp.GetCoords(), radius, pen)

        return True

    def ColorAtomByFragmentIndex(self):
        return ColorAtomByFragmentIndex(self.colorlist, self.tag).__disown__()


class LabelBondOrder(oedepict.OEDisplayBondPropBase):
    def __init__(self):
        oedepict.OEDisplayBondPropBase.__init__(self)

    def __call__(self, bond):
        if 'WibergBondOrder_frag' in bond.GetData():
            bondOrder = bond.GetData('WibergBondOrder_frag')
            label = "{:.2f}".format(bondOrder)
            return label
        return ' '

    def CreateCopy(self):
        copy = LabelBondOrder()
        return copy.__disown__()