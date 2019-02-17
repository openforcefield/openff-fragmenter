from itertools import combinations
import openeye as oe
from openeye import oechem, oedepict, oegrapheme, oequacpac, oeomega, oeiupac
from cmiles.utils import mol_to_smiles, has_stereo_defined, has_atom_map, remove_atom_map, restore_atom_map

import yaml
import os
from pkg_resources import resource_filename
import copy
import itertools
import json
import warnings
import networkx as nx

from .utils import logger, make_python_identifier
from .chemi import to_smi, normalize_molecule, get_charges


OPENEYE_VERSION = oe.__name__ + '-v' + oe.__version__


def expand_states(inp_molecule, protonation=True, tautomers=False, stereoisomers=True, max_states=200, level=0, reasonable=True,
                  carbon_hybridization=True, suppress_hydrogen=True, verbose=False, filename=None,
                  return_smiles_list=False, return_molecules=False, expand_internally=True, strict=False):
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
                                        suppress_hydrogen=False))

    # if stereoisomers:
    #     molecules.remove(inp_molecule)
    for molecule in molecules:
        #states.add(fragmenter.utils.create_mapped_smiles(molecule, tagged=False, explicit_hydrogen=False))
        # Not using create mapped SMILES because OEMol is needed but state is OEMolBase.
        #states.add(oechem.OEMolToSmiles(molecule))
        try:
            states.add(mol_to_smiles(molecule, isomeric=True, mapped=False, explicit_hydrogen=True))
        except ValueError:
            if stereoisomers:
                continue
            elif not strict:
                states.add(mol_to_smiles(molecule, isomeric=False, mapped=False, explicit_hydrogen=True ))
            elif expand_internally:
                logger().warn("Tautomer or protonation state has a chiral center. Expanding stereoisomers")
                stereo_states = _expand_states(molecule, enumerate='steroisomers')
                for state in stereo_states:
                    states.add(mol_to_smiles(molecule, isomeric=True, mapped=False, explicit_hydrogen=False))
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
        if has_atom_map(molecule):
            remove_atom_map(molecule, keep_map_data=True)
        self.molecule = molecule
        self._tag_fgroups()
        self._nx_graph = self._mol_to_graph()
        self._fragments = list()  # all possible fragments without breaking rings

        self._fragment_combinations = list() # AtomBondSets of combined fragments. Used internally to generate PDF
        self.fragment_combinations = {} # Dict that maps SMILES to all equal combined fragments

    @property
    def n_rotors(self):
        return sum([bond.IsRotor() for bond in self.molecule.GetBonds()])

    def _tag_fgroups(self, fgroups_smarts=None):
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
        if not fgroups_smarts:
            # Load yaml file
            fn = resource_filename('fragmenter', os.path.join('data', 'fgroup_smarts_comb.yml'))
            f = open(fn, 'r')
            fgroups_smarts = yaml.safe_load(f)
            f.close()

        fgroup_tagged = {}
        for f_group in fgroups_smarts:
            qmol = oechem.OEQMol()
            if not oechem.OEParseSmarts(qmol, fgroups_smarts[f_group]):
                print('OEParseSmarts failed')
            ss = oechem.OESubSearch(qmol)
            oechem.OEPrepareSearch(self.molecule, ss)

            for i, match in enumerate(ss.Match(self.molecule, True)):
                fgroup_atoms = set()
                for ma in match.GetAtoms():
                    fgroup_atoms.add(ma.target.GetIdx())
                    tag = oechem.OEGetTag('fgroup')
                    ma.target.SetData(tag, '{}_{}'.format(f_group, str(i)))
                fgroup_bonds = set()
                for ma in match.GetBonds():
                    #if not ma.target.IsInRing():
                    fgroup_bonds.add(ma.target.GetIdx())
                    tag =oechem.OEGetTag('fgroup')
                    ma.target.SetData(tag, '{}_{}'.format(f_group, str(i)))

                fgroup_tagged['{}_{}'.format(f_group, str(i))] = (fgroup_atoms, fgroup_bonds)
        return fgroup_tagged

    def _mol_to_graph(self):

        G = nx.Graph()
        for atom in self.molecule.GetAtoms():
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
        #ToDo - add option for fgroups
        subgraphs = self._fragment_graph()
        self._subgraphs_to_atom_bond_sets(subgraphs)

    def combine_fragments(self, max_rotors=3, min_rotors=1, min_heavy_atoms=5, **kwargs):
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
                            mol = self.frag_to_mol(frag, **kwargs)
                            smiles = []
                            if not isinstance(mol, list):
                                mol = [mol]
                            for m in mol:
                                smiles.append(mol_to_smiles(m, isomeric=True, mapped=False, explicit_hydrogen=True))
                            for sm in smiles:
                                if sm not in self.fragment_combinations:
                                    self.fragment_combinations[sm] = []
                                self.fragment_combinations[sm].extend(mol)

    def frag_to_mol(self, frag, adjust_hcount=True, explicit_hydrogens=True, expand_stereoisomers=True,
                       restore_maps=False):
        """
        Convert fragments (AtomBondSet) to OEMol
        Parameters
        ----------
        frags: list
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

        if restore_maps:
            restore_atom_map(fragment)
        # check for stereo defined
        if not has_stereo_defined(fragment):
            # Try to convert to smiles and back. A molecule might look like it's missing stereo because of submol
            new_smiles = oechem.OEMolToSmiles(fragment)
            fragment = oechem.OEMol()
            oechem.OESmilesToMol(fragment, new_smiles)
            # add explicit H
            oechem.OEAddExplicitHydrogens(fragment)
            oechem.OEPerceiveChiral(fragment)
            # If it's still missing stereo, expand states
            if not has_stereo_defined(fragment):
                if expand_stereoisomers:
                    enantiomers = _expand_states(fragment, enumerate='stereoisomers')
                    return enantiomers

        return fragment

        # try:
        #     #ToDo generate smiles with restored map
        #     smiles = mol_to_smiles(fragment, mapped=False, explicit_hydrogen=True, isomeric=True)
        # except ValueError:
        #     if expand_stereoisomers:
        #         # Generate stereoisomers
        #         smiles = list()
        #         enantiomers = _expand_states(fragment, enumerate='stereoisomers')
        #         for mol in enantiomers:
        #             smiles.append(mol_to_smiles(mol, mapped=False, explicit_hydrogen=True, isomeric=True))
        #     else:
        #         # generate non isomeric smiles
        #         smiles = mol_to_smiles(fragment, mapped=False, explicit_hydrogen=explicit_hydrogens, isomeric=False)
        #
        # return smiles

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

    @staticmethod
    def _count_rotors_in_fragment(fragment):
        return sum([bond.IsRotor() for bond in fragment.GetBonds()])

    @staticmethod
    def _count_heavy_atoms_in_fragment(fragment):
        return sum([not atom.IsHydrogen() for atom in fragment.GetAtoms()])

    def depict_fragment_combinations(self, fname, line_width=0.75):

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


def generate_fragments(molecule, generate_visualization=False, strict_stereo=False, combinatorial=True, MAX_ROTORS=2,
                       remove_map=True, json_filename=None):
    """
    This function generates fragments from molecules. The output is a dictionary that maps SMILES of molecules to SMILES
     for fragments. The default SMILES are generated with openeye.oechem.OEMolToSmiles. These SMILES strings are canonical
     isomeric SMILES.
     The dictionary also includes a provenance field which defines how the fragments were generated.

    Parameters
    ----------
    molecule: OEMol to fragment
    generate_visualization: bool
        If true, visualization of the fragments will be written to pdf files. The pdf will be writtten in the directory
        where this function is run from.
    combinatorial: bool
        If true, find all connected fragments from fragments and add all new fragments that have less than MAX_ROTORS
    MAX_ROTORS: int
        rotor threshold for combinatorial
    strict_stereo: bool
        Note: This applies to the molecule being fragmented. Not the fragments.
        If True, omega will generate conformation with stereochemistry defined in the SMILES string for charging.
    remove_map: bool
        If True, the index tags will be removed. This will remove duplicate fragments. Defualt True
    json_filename: str
        filenmae for JSON. If provided, will save the returned dictionary to a JSON file. Default is None

    Returns
    -------
    fragments: dict
        mapping of SMILES from the parent molecule to the SMILES of the fragments
    """
    fragments = dict()

    try:
        molecules = list(molecule)
    except TypeError:
        molecules = [molecule]
    for molecule in molecules:
        # normalize molecule
        molecule = normalize_molecule(molecule, molecule.GetTitle())
        if remove_map:
            # Remove tags from smiles. This is done to make it easier to find duplicate fragments
            for a in molecule.GetAtoms():
                a.SetMapIdx(0)
        frags = _generate_fragments(molecule, strict_stereo=strict_stereo)
        if not frags:
            logger().warning('Skipping {}, SMILES: {}'.format(molecule.GetTitle(), oechem.OECreateSmiString(molecule)))
            continue
        charged = frags[0]
        frags = frags[-1]
        frag_list = list(frags.values())
        if combinatorial:
            smiles = smiles_with_combined(frag_list, charged, MAX_ROTORS)
        else:
            smiles = frag_to_smiles(frag_list, charged)

        parent_smiles = mol_to_smiles(molecule, isomeric=True, explicit_hydrogen=False, mapped=False)
        if smiles:
            fragments[parent_smiles] = list(smiles.keys())
        else:
            # Add molecule where no fragments were found for terminal torsions and / or rings and non rotatable bonds
            fragments[parent_smiles] = [mol_to_smiles(molecule, isomeric=True, explicit_hydrogen=True, mapped=False)]

        if generate_visualization:
            IUPAC = oeiupac.OECreateIUPACName(molecule)
            name = molecule.GetTitle()
            if IUPAC == name:
                name = make_python_identifier(oechem.OEMolToSmiles(molecule))[0]
            oname = '{}.pdf'.format(name)
            ToPdf(charged, oname, frags)
        del charged, frags
    if json_filename:
        f = open(json_filename, 'w')
        j = json.dump(fragments, f, indent=2, sort_keys=True)
        f.close()

    return fragments


def _generate_fragments(mol, strict_stereo=True):
    """
    This function generates fragments from a molecule.

    Parameters
    ----------
    mol: OEMol
    strict_stereo: bool
        If False, omega will generate conformer without the specific stereochemistry

    Returns
    -------
    charged: charged OEMOl
    frags: dict of AtomBondSet mapped to rotatable bond index the fragment was built up from.
    """

    try:
        charged = get_charges(mol, keep_confs=1, strict_stereo=strict_stereo)
    except RuntimeError:
        logger().warning("Could not charge molecule {} so no WBO were calculated. Cannot fragment molecule {}".format(mol.GetTitle(),
                                                                                                                      mol.GetTitle()))
        return False

    # Check if WBO were calculated
    bonds = [bond for bond in charged.GetBonds()]
    for bond in bonds[:1]:
        try:
            bond.GetData('WibergBondOrder')
        except ValueError:
            logger().warning("WBO were not calculate. Cannot fragment molecule {}".format(charged.GetTitle()))
            return False

    tagged_rings, tagged_fgroups = tag_molecule(charged)

    # Iterate over bonds
    frags = {}
    for bond in charged.GetBonds():
        if bond.IsRotor():
            atoms, bonds = _build_frag(bond=bond, mol=charged, tagged_fgroups=tagged_fgroups, tagged_rings=tagged_rings)
            atom_bond_set = _to_AtomBondSet(charged, atoms, bonds)
            frags[bond.GetIdx()] = atom_bond_set

    return charged, frags


def _tag_fgroups(mol, fgroups_smarts=None):
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
    if not fgroups_smarts:
        # Load yaml file
        fn = resource_filename('fragmenter', os.path.join('data', 'fgroup_smarts.yml'))
        f = open(fn, 'r')
        fgroups_smarts = yaml.safe_load(f)
        f.close()

    fgroup_tagged = {}
    for f_group in fgroups_smarts:
        qmol = oechem.OEQMol()
        if not oechem.OEParseSmarts(qmol, fgroups_smarts[f_group]):
            print('OEParseSmarts failed')
        ss = oechem.OESubSearch(qmol)
        oechem.OEPrepareSearch(mol, ss)

        for i, match in enumerate(ss.Match(mol, True)):
            fgroup_atoms = set()
            for ma in match.GetAtoms():
                fgroup_atoms.add(ma.target.GetIdx())
                tag = oechem.OEGetTag('fgroup')
                ma.target.SetData(tag, '{}_{}'.format(f_group, str(i)))
            fgroup_bonds = set()
            for ma in match.GetBonds():
                #if not ma.target.IsInRing():
                fgroup_bonds.add(ma.target.GetIdx())
                tag =oechem.OEGetTag('fgroup')
                ma.target.SetData(tag, '{}_{}'.format(f_group, str(i)))

            fgroup_tagged['{}_{}'.format(f_group, str(i))] = (fgroup_atoms, fgroup_bonds)
    return fgroup_tagged


def _tag_rings(mol):
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
    tagged_rings = {}
    nringsystems, parts = oechem.OEDetermineRingSystems(mol)
    for ringidx in range(1, nringsystems +1):
        ringidx_atoms = set()
        for atom in mol.GetAtoms():
            if parts[atom.GetIdx()] == ringidx:
                ringidx_atoms.add(atom.GetIdx())
                tag = oechem.OEGetTag('ringsystem')
                atom.SetData(tag, ringidx)
        # Find bonds in ring and tag
        ringidx_bonds = set()
        for a_idx in ringidx_atoms:
            atom = mol.GetAtom(oechem.OEHasAtomIdx(a_idx))
            for bond in atom.GetBonds():
                nbrAtom = bond.GetNbr(atom)
                nbrIdx = nbrAtom.GetIdx()
                if nbrIdx in ringidx_atoms and nbrIdx != a_idx:
                    ringidx_bonds.add(bond.GetIdx())
                    tag = oechem.OEGetTag('ringsystem')
                    bond.SetData(tag, ringidx)
        tagged_rings[ringidx] = (ringidx_atoms, ringidx_bonds)
    return tagged_rings


def _ring_fgroup_union(mol, tagged_rings, tagged_fgroups, wbo_threshold=1.2):
    """
    This function combines rings and fgroups that are conjugated (the bond between them has a Wiberg bond order > 1.2)

    Parameters
    ----------
    mol: OpenEye OEMolGraph
    tagged_rings: dict
        map of ringsystem indices to ring atom and bond indices
    tagged_fgroup: dict
        map of fgroup to fgroup atom and bond indices

    Returns
    -------
    tagged_fgroup: dict
        updated tagged_fgroup mapping with rings that shouldn't be fragmented from fgroups
    """
    ring_idxs = list(tagged_rings.keys())
    fgroups = list(tagged_fgroups.keys())
    tagged_fgroups = copy.deepcopy(tagged_fgroups)

    # Check if fgroups are overlapping. If they are - combine them.
    for func_group_1, func_group_2 in itertools.combinations(list(tagged_fgroups.keys()), 2):
        atoms_intersection = tagged_fgroups[func_group_1][0].intersection(tagged_fgroups[func_group_2][0])
        if len(atoms_intersection) > 1:
            # Combine fgroups
            atoms_union = tagged_fgroups[func_group_1][0].union(tagged_fgroups[func_group_2][0])
            bonds_union = tagged_fgroups[func_group_1][-1].union(tagged_fgroups[func_group_2][-1])
            tagged_fgroups[func_group_1] = (atoms_union, bonds_union)
            tagged_fgroups[func_group_2] = (atoms_union, bonds_union)
    for idx in ring_idxs:
        for fgroup in fgroups:
            atom_intersection = tagged_rings[idx][0].intersection(tagged_fgroups[fgroup][0])
            if len(atom_intersection) > 1:
                # Must include ring if including fgroup. Add ring atoms and bonds to fgroup
                atoms_union = tagged_rings[idx][0].union(tagged_fgroups[fgroup][0])
                bonds_union = tagged_rings[idx][-1].union(tagged_fgroups[fgroup][-1])
                tagged_fgroups[fgroup] = (atoms_union, bonds_union)
            elif len(atom_intersection) > 0:
                # Check Wiberg bond order of bond
                # First find bond connectiong fgroup and ring
                atom = mol.GetAtom(oechem.OEHasAtomIdx(atom_intersection.pop()))
                for a in atom.GetAtoms():
                    if a.GetIdx() in tagged_fgroups[fgroup][0]:
                        bond = mol.GetBond(a, atom)
                        if bond.GetData('WibergBondOrder') > wbo_threshold:
                            # Don't cut off ring.
                            atoms_union = tagged_rings[idx][0].union(tagged_fgroups[fgroup][0])
                            bonds_union = tagged_rings[idx][-1].union(tagged_fgroups[fgroup][-1])
                            tagged_fgroups[fgroup] = (atoms_union, bonds_union)
                        # Should I also combine non-rotatable rings? This will pick up the alkyn in ponatinib?
                        if not bond.IsRotor():
                            atoms_union = tagged_rings[idx][0].union(tagged_fgroups[fgroup][0])
                            bonds_union = tagged_rings[idx][-1].union(tagged_fgroups[fgroup][-1])
                            tagged_fgroups[fgroup] = (atoms_union, bonds_union)
                        # Do something when it's neither for edge cases (Afatinib)?
    return tagged_fgroups


def tag_molecule(mol, func_group_smarts=None):
    """
    Tags atoms and molecules in functional groups and ring systems. The molecule gets tagged and the function returns
    a 2 dictionaries that map
    1) ring system indices to corresponding atoms and bonds indices in molecule
    2) functional groups to corresponding atoms and bonds indices in molecule

    Parameters
    ----------
    mol: OEMol
    func_group_smarts: dict
        dictionary mapping functional groups to SMARTS. Default is None and uses shipped yaml file.

    Returns
    -------
    tagged_rings: dict
        mapping of ring system index to corresponding atoms and bonds in molecule. Each index maps to 2 sets. First set
        includes atom indices and second set includes bonds indices
    tagged_func_group: dict
        mapping of functional group to corresponding atoms and bonds in molecule. Each functional group maps to 2 sets.
        The first set is atom indices, the second set is bond indices.

    """
    tagged_func_group = _tag_fgroups(mol, func_group_smarts)
    tagged_rings = _tag_rings(mol)

    tagged_func_group = _ring_fgroup_union(mol=mol, tagged_fgroups=tagged_func_group, tagged_rings=tagged_rings)

    return tagged_rings, tagged_func_group


def _is_fgroup(fgroup_tagged, element):
    """
    This function checks if an atom or a bond is part of a tagged fgroup.

    Parameters
    ----------
    atom: Openeye Atom Base
    fgroup_tagged: dict of indexed functional group and corresponding atom and bond indices

    Returns
    -------
    atoms, bonds: sets of atom and bond indices if the atom is tagged, False otherwise

    """
    try:
        fgroup = element.GetData('fgroup')
        atoms, bonds = fgroup_tagged[fgroup]
        return atoms, bonds
    except ValueError:
        return False


def _to_AtomBondSet(mol, atoms, bonds):
    """
    Builds OpeneyeAtomBondet from atoms and bonds set of indices
    Parameters
    ----------
    mol: Openeye OEMolGraph
    atoms: Set of atom indices
    bonds: Set of bond indices

    Returns
    -------
    AtomBondSet: Openeye AtomBondSet of fragment
    """

    AtomBondSet = oechem.OEAtomBondSet()
    for a_idx in atoms:
        AtomBondSet.AddAtom(mol.GetAtom(oechem.OEHasAtomIdx(a_idx)))
    for b_idx in bonds:
        AtomBondSet.AddBond(mol.GetBond(oechem.OEHasBondIdx(b_idx)))
    return AtomBondSet


def _is_ortho(bond, rot_bond, next_bond):
    """
    This function checks if a bond is ortho to the rotatable bond
    Parameters
    ----------
    bond: OEBondBase
        current bond to check if it's ortho
    rot_bond: OEBondBase
        the rotatable bond the bond needs to be ortho to
    next_bond: OEBondBase
        The bond between rot_bond and bond if bond is ortho to rot_bond

    Returns
    -------
    bool: True if ortho, False if not
    """

    bond_attached = set()
    rot_attached = set()

    for bon in [bond.GetBgn(), bond.GetEnd()]:
        for b in bon.GetBonds():
            bond_attached.add(b.GetIdx())
    for bon in [rot_bond.GetBgn(), rot_bond.GetEnd()]:
        for b in bon.GetBonds():
            rot_attached.add(b.GetIdx())

    if not next_bond.IsInRing():
        next_attached = set()
        for bon in [next_bond.GetBgn(), next_bond.GetEnd()]:
            for b in bon.GetBonds():
                next_attached.add(b.GetIdx())

    intersection = (bond_attached & rot_attached)
    if not bool(intersection) and not next_bond.IsInRing():
        # Check if it's ortho to next bond
        intersection = (bond_attached & next_attached)
    return bool(intersection)


def _build_frag(bond, mol, tagged_fgroups, tagged_rings):
    """
    This functions builds a fragment around a rotatable bond. It grows out one bond in all directions
    If the next atoms is in a ring or functional group, it keeps that.
    If the next bond has a Wiberg bond order > 1.2, grow another bond and check next bond's Wiberg bond order.

    Parameters
    ----------
    bond: OpenEye bond
    mol: OpenEye OEMolGraph
    tagged_fgroups: dict
        maps functional groups to atoms and bond indices on mol
    tagged_rings: dict
        maps ringsystem index to atom and bond indices in mol

    Returns
    -------
    atoms, bonds: sets of atom and bond indices for fragment
    """

    atoms = set()
    bonds = set()
    b_idx = bond.GetIdx()
    bonds.add(b_idx)
    beg = bond.GetBgn()
    end = bond.GetEnd()
    beg_idx = beg.GetIdx()
    end_idx = end.GetIdx()

    atoms.add(beg_idx)
    atoms_nb, bonds_nb = iterate_nbratoms(mol=mol, rotor_bond=bond, atom=beg, pair=end, fgroup_tagged=tagged_fgroups,
                                          tagged_rings=tagged_rings)
    atoms = atoms.union(atoms_nb)
    bonds = bonds.union(bonds_nb)

    atoms.add(end_idx)
    atoms_nb, bonds_nb = iterate_nbratoms(mol=mol, rotor_bond=bond, atom=end, pair=beg, fgroup_tagged=tagged_fgroups,
                                          tagged_rings=tagged_rings)
    atoms = atoms.union(atoms_nb)
    bonds = bonds.union(bonds_nb)

    return atoms, bonds


def iterate_nbratoms(mol, rotor_bond, atom, pair, fgroup_tagged, tagged_rings, i=0, threshold=1.05):
    """
    This function iterates over neighboring atoms and checks if it's part of a functional group, ring, or if the next
    bond has a Wiberg bond order > 1.2.

    Parameters
    ----------
    mol: Openeye OEMolGraph
    atom: Openeye AtomBase
        atom that will iterate over
    paired: OpeneEye AtomBase
        atom that's bonded to this atom in rotor_bond
    fgroup_tagged: dict
        map of functional group and atom and bond indices in mol
    tagged_rings: dict
        map of ringsystem index and atom and bond indices in mol
    rotor_bond: Openeye Bond base
        rotatable bond that the fragment is being built on

    Returns
    -------
    atoms, bonds: sets of atom and bond indices of the fragment

    """
    def _iterate_nbratoms(mol, rotor_bond, atom, pair, fgroup_tagged, tagged_rings, atoms_2, bonds_2, i=0, threshold=threshold):

        for a in atom.GetAtoms():
            if a.GetIdx() == pair.GetIdx():
                continue
            a_idx = a.GetIdx()
            next_bond = mol.GetBond(a, atom)
            nb_idx = next_bond.GetIdx()
            # atoms_2.add(a_idx)
            # bonds_2.add(nb_idx)
            if nb_idx in bonds_2:
                FGROUP_RING = False
                try:
                    ring_idx = a.GetData('ringsystem')
                    fgroup = a.GetData('fgroup')
                    FGROUP_RING = True
                except ValueError:
                    try:
                        ring_idx = atom.GetData('ringsystem')
                        fgroup = atom.GetData('fgroup')
                        FGROUP_RING = True
                    except ValueError:
                        continue
                if FGROUP_RING:
                    # Add ring and continue
                    ratoms, rbonds = tagged_rings[ring_idx]
                    atoms_2 = atoms_2.union(ratoms)
                    bonds_2 = bonds_2.union(rbonds)
                    rs_atoms, rs_bonds = _ring_substiuents(mol=mol, bond=next_bond, rotor_bond=rotor_bond,
                                                          tagged_rings=tagged_rings, ring_idx=ring_idx,
                                                          fgroup_tagged=fgroup_tagged)
                    atoms_2 = atoms_2.union(rs_atoms)
                    bonds_2 = bonds_2.union(rs_bonds)
                    continue

            if i > 0:
                wiberg = next_bond.GetData('WibergBondOrder')
                if wiberg < threshold:
                    #atoms_2.remove(a_idx)
                    continue

            atoms_2.add(a_idx)
            bonds_2.add(nb_idx)
            if a.IsInRing():
                ring_idx = a.GetData('ringsystem')
                ratoms, rbonds = tagged_rings[ring_idx]
                atoms_2 = atoms_2.union(ratoms)
                bonds_2 = bonds_2.union(rbonds)
                # Find non-rotatable sustituents
                rs_atoms, rs_bonds = _ring_substiuents(mol=mol, bond=next_bond, rotor_bond=rotor_bond,
                                                      tagged_rings=tagged_rings, ring_idx=ring_idx,
                                                      fgroup_tagged=fgroup_tagged)
                atoms_2 = atoms_2.union(rs_atoms)
                bonds_2 = bonds_2.union(rs_bonds)
            fgroup = _is_fgroup(fgroup_tagged, element=a)
            if fgroup: # and i < 1:
                atoms_2 = atoms_2.union(fgroup[0])
                bonds_2 = bonds_2.union(fgroup[-1])
                # if something is in a ring - have a flag? Then use that to continue iterating and change the flag

            for nb_a in a.GetAtoms():
                nn_bond = mol.GetBond(a, nb_a)
                if (nn_bond.GetData('WibergBondOrder') > threshold) and (not nn_bond.IsInRing()) and (not nn_bond.GetIdx() in bonds_2):
                    # Check the degree of the atoms in the bond
                    deg_1 = a.GetDegree()
                    deg_2 = nb_a.GetDegree()
                    if deg_1 == 1 or deg_2 == 1:
                        continue

                    atoms_2.add(nb_a.GetIdx())
                    bonds_2.add(nn_bond.GetIdx())
                    i += 1
                    _iterate_nbratoms(mol, nn_bond, nb_a, pair, fgroup_tagged, tagged_rings, atoms_2, bonds_2, i=i)
        return atoms_2, bonds_2
    return _iterate_nbratoms(mol, rotor_bond, atom, pair, fgroup_tagged, tagged_rings, atoms_2=set(), bonds_2=set(), i=0)


def _ring_substiuents(mol, bond, rotor_bond, tagged_rings, ring_idx, fgroup_tagged):
    """
    This function finds ring substituents that shouldn't be cut off

    Parameters
    ----------
    mol: OpeneEye OEMolGraph
    bond: OpenEye Bond Base
        current bond that the iterator is looking at
    rotor_bond: OpeneEye Bond Base
        rotatable bond that fragment is being grown on
    tagged_rings: dict
        mapping of ring index and atom and bonds indices
    ring_idx: int
        ring index
    fgroup_tagged: dict
        mapping of functional group and atom and bond indices

    Returns
    -------
    rs_atoms, rs_bonds: sets of ring substituents atoms and bonds indices

    """
    rs_atoms = set()
    rs_bonds = set()
    r_atoms, r_bonds = tagged_rings[ring_idx]
    for a_idx in r_atoms:
        atom = mol.GetAtom(oechem.OEHasAtomIdx(a_idx))
        for a in atom.GetAtoms():
            if a.GetIdx() in rs_atoms:
                continue
            fgroup = False
            rs_bond = mol.GetBond(atom, a)
            if not a.IsInRing():
                if not rs_bond.IsRotor():
                    rs_atoms.add(a.GetIdx())
                    rs_bonds.add(rs_bond.GetIdx())
                    # Check for functional group
                    fgroup = _is_fgroup(fgroup_tagged, element=a)
                elif _is_ortho(rs_bond, rotor_bond, bond):
                    # Keep bond and attached atom.
                    rs_atoms.add(a.GetIdx())
                    rs_bonds.add(rs_bond.GetIdx())
                    # Check for functional group
                    fgroup = _is_fgroup(fgroup_tagged, element=a)
                if fgroup:
                    rs_atoms = rs_atoms.union(fgroup[0])
                    rs_bonds = rs_bonds.union(fgroup[-1])
            else:
                # Check if they are in the same ring
                r_idx2 = a.GetData('ringsystem')
                if r_idx2 != ring_idx:
                    if _is_ortho(rs_bond, rotor_bond, bond):
                        # Add ring system
                        rs_bonds.add(rs_bond.GetIdx())
                        r2_atoms, r2_bonds = tagged_rings[r_idx2]
                        rs_atoms = rs_atoms.union(r2_atoms)
                        rs_bonds = rs_bonds.union(r2_bonds)

    return rs_atoms, rs_bonds


def frag_to_smiles(frags, mol, adjust_hcount=True, expand_stereoisomers=True):
    """
    Convert fragments (AtomBondSet) to canonical isomeric SMILES string
    Parameters
    ----------
    frags: list
    mol: OEMol
    OESMILESFlag: str
        Either 'ISOMERIC' or 'DEFAULT'. This flag determines which OE function to use to generate SMILES string

    Returns
    -------
    smiles: dict of smiles to frag

    """

    smiles = {}
    for frag in frags:
        fragatompred = oechem.OEIsAtomMember(frag.GetAtoms())
        fragbondpred = oechem.OEIsBondMember(frag.GetBonds())

        #fragment = oechem.OEGraphMol()
        fragment = oechem.OEMol()
        adjustHCount = adjust_hcount
        oechem.OESubsetMol(fragment, mol, fragatompred, fragbondpred, adjustHCount)

        # Add explicit H
        oechem.OEAddExplicitHydrogens(fragment)
        oechem.OEPerceiveChiral(fragment)
        # sanity check that all atoms are bonded
        for atom in fragment.GetAtoms():
            if not list(atom.GetBonds()):
                warnings.warn("Yikes!!! An atom that is not bonded to any other atom in the fragment. "
                              "You probably ran into a bug. Please report the input molecule to the issue tracker")
        #s = oechem.OEMolToSmiles(fragment)
        #s2 = fragmenter.utils.create_mapped_smiles(fragment, tagged=False, explicit_hydrogen=False)
        if not has_stereo_defined(fragment):
            # Try to convert to smiles and back. A molecule might look like it's missing stereo because of submol
            new_smiles = oechem.OEMolToSmiles(fragment)
            fragment = oechem.OEMol()
            oechem.OESmilesToMol(fragment, new_smiles)
            # add explicit H
            oechem.OEAddExplicitHydrogens(fragment)
            oechem.OEPerceiveChiral(fragment)
        try:
            s = mol_to_smiles(fragment, mapped=False, explicit_hydrogen=True, isomeric=True)
        except ValueError:
            if expand_stereoisomers:
                # Generate stereoisomers
                s = list()
                enantiomers = _expand_states(fragment, enumerate='stereoisomers')
                for mol in enantiomers:
                    sm = mol_to_smiles(mol, mapped=False, explicit_hydrogen=True, isomeric=True)
                    s.append(sm)
            else:
                # generate non isomeric smiles
                s = mol_to_smiles(fragment, mapped=False, explicit_hydrogen=explicit_hydrogens, isomeric=False)

        #s = mol_to_smiles(fragment, mapped=False, explicit_hydrogen=True, isomeric=True)
        if not isinstance(s, list):
            if not s:
                continue
            s = [s]
        for sm in s:
            if sm not in smiles:
                smiles[sm] = []
            smiles[sm].append(frag)

    return smiles


def _sort_by_rotbond(ifs, outdir):
    """

    Parameters
    ----------
    ifs: str
        absolute path to molecule database
    outdir: str
        absolute path to where output files should be written.
    """

    nrotors_map = {}
    moldb = oechem.OEMolDatabase(ifs)
    mol = oechem.OEGraphMol()
    for idx in range(moldb.GetMaxMolIdx()):
        if moldb.GetMolecule(mol, idx):
            nrotors = sum([bond.IsRotor() for bond in mol.GetBonds()])
            if nrotors not in nrotors_map:
                nrotors_map[nrotors] = []
            nrotors_map[nrotors].append(idx)
    # Write out a separate database for each num rotors
    for nrotor in nrotors_map:
        size = len(nrotors_map[nrotor])
        ofname = os.path.join(outdir, 'nrotor_{}.smi'.format(nrotor))
        ofs = new_output_stream(ofname)
        write_oedatabase(moldb, ofs, nrotors_map[nrotor], size)


def smiles_with_combined(frag_list, mol, MAX_ROTORS=2):
    """
    Generates Smiles:frags mapping for fragments and fragment combinations with less than MAX_ROTORS rotatable bonds

    Parameters
    ----------
    frags: dict of rot band mapped to AtomBondSet
    mol: OpenEye Mol
    OESMILESFlag: str
        Either 'ISOMERIC' or 'DEFAULT'. This flag determines which OE function to use to generate SMILES string

    Returns
    -------
    smiles: dict of smiles sting to fragment

    """
    comb_list = GetFragmentAtomBondSetCombinations(frag_list, MAX_ROTORS=MAX_ROTORS)

    combined_list = comb_list + frag_list

    smiles = frag_to_smiles(combined_list, mol)

    return smiles


def ToPdf(mol, oname, frags):#, fragcombs):
    """
    Parameters
    ----------
    mol: charged OEMolGraph
    oname: str
        Output file name
    Returns
    -------

    """
    itf = oechem.OEInterface()
    oedepict.OEPrepareDepiction(mol)

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

    DepictMoleculeWithFragmentCombinations(report, mol, frags, opts)

    return oedepict.OEWriteReport(oname, report)

    #return 0


def SmilesToFragments(smiles, fgroup_smarts, bondOrderThreshold=1.2, chargesMol=True):
    """
    Fragment molecule at bonds below Bond Order Threshold

    Parameters
    ----------
    smiles: str
        smiles string of molecule to fragment

    Returns
    -------
    frags: list of OE AtomBondSets

    """
    # Charge molecule
    mol = oechem.OEGraphMol()
    oemol = openeye.smiles_to_oemol(smiles)
    charged = openeye.get_charges(oemol, keep_confs=1)

    # Tag functional groups
    _tag_fgroups(charged, fgroups_smarts=fgroup_smarts)

    # Generate fragments
    G = OeMolToGraph(charged)
    subraphs = FragGraph(G, bondOrderThreshold=bondOrderThreshold)

    frags = []
    for subraph in subraphs:
        frags.append(subgraphToAtomBondSet(G, subraph, charged))

    if chargesMol:
        return frags, charged
    else:
        return frags


def DepictMoleculeWithFragmentCombinations(report, mol, frags, opts): #fragcombs, opts):
    """ This function was taken from https://docs.eyesopen.com/toolkits/cookbook/python/depiction/enumfrags.html with some modification
    """
    stag = "fragment idx"
    itag = oechem.OEGetTag(stag)
    for fidx, frag in enumerate(frags):
        print(fidx, frag)
        for bond in frags[frag].GetBonds():
            bond.SetData(itag, fidx)

    # setup depiction styles

    nrfrags = len(frags)
    colors = [c for c in oechem.OEGetLightColors()]
    if len(colors) < nrfrags:
        colors = [c for c in oechem.OEGetColors(oechem.OEYellowTint, oechem.OEDarkOrange, nrfrags)]

    bondglyph = ColorBondByFragmentIndex(colors, itag)

    lineWidthScale = 0.75
    fadehighlight = oedepict.OEHighlightByColor(oechem.OEGrey, lineWidthScale)

    # depict each fragment combinations

    for frag in frags:

        cell = report.NewCell()
        disp = oedepict.OE2DMolDisplay(mol, opts)

        fragatoms = oechem.OEIsAtomMember(frags[frag].GetAtoms())
        fragbonds = oechem.OEIsBondMember(frags[frag].GetBonds())

        notfragatoms = oechem.OENotAtom(fragatoms)
        notfragbonds = oechem.OENotBond(fragbonds)

        oedepict.OEAddHighlighting(disp, fadehighlight, notfragatoms, notfragbonds)

        bond = mol.GetBond(oechem.OEHasBondIdx(frag))

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

    bondlabel = LabelBondOrder()
    opts.SetBondPropertyFunctor(bondlabel)
    disp = oedepict.OE2DMolDisplay(mol, opts)
    #oegrapheme.OEAddGlyph(disp, bondglyph, oechem.IsTrueBond())

    headerpen = oedepict.OEPen(oechem.OEWhite, oechem.OELightGrey, oedepict.OEFill_Off, 2.0)
    for header in report.GetHeaders():
        oedepict.OERenderMolecule(header, disp)
        oedepict.OEDrawBorder(header, headerpen)


class ColorBondByFragmentIndex(oegrapheme.OEBondGlyphBase):
    """
    This class was taken from OpeneEye cookbook
    https://docs.eyesopen.com/toolkits/cookbook/python/depiction/enumfrags.html
    """
    def __init__(self, colorlist, tag):
        oegrapheme.OEBondGlyphBase.__init__(self)
        self.colorlist = colorlist
        self.tag = tag

    def RenderGlyph(self, disp, bond):

        bdisp = disp.GetBondDisplay(bond)
        if bdisp is None or not bdisp.IsVisible():
            return False

        if not bond.HasData(self.tag):
            return False

        linewidth = disp.GetScale() / 2.0
        color = self.colorlist[bond.GetData(self.tag)]
        pen = oedepict.OEPen(color, color, oedepict.OEFill_Off, linewidth)

        adispB = disp.GetAtomDisplay(bond.GetBgn())
        adispE = disp.GetAtomDisplay(bond.GetEnd())

        layer = disp.GetLayer(oedepict.OELayerPosition_Below)
        layer.DrawLine(adispB.GetCoords(), adispE.GetCoords(), pen)

        return True

    def ColorBondByFragmentIndex(self):
        return ColorBondByFragmentIndex(self.colorlist, self.tag).__disown__()


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
        bondOrder = bond.GetData('WibergBondOrder')
        label = "{:.2f}".format(bondOrder)
        return label

    def CreateCopy(self):
        copy = LabelBondOrder()
        return copy.__disown__()

def IsAdjacentAtomBondSets(fragA, fragB):
    """
    This function was taken from Openeye cookbook
    https://docs.eyesopen.com/toolkits/cookbook/python/cheminfo/enumfrags.html
    """
    for atomA in fragA.GetAtoms():
        for atomB in fragB.GetAtoms():
            if atomA.GetBond(atomB) is not None:
                return True
    return False

def CountRotorsInFragment(fragment):
    return sum([bond.IsRotor() for bond in fragment.GetBonds()])

def IsAdjacentAtomBondSetCombination(fraglist):
    """
    This function was taken from Openeye cookbook
    https://docs.eyesopen.com/toolkits/cookbook/python/cheminfo/enumfrags.html
    """

    parts = [0] * len(fraglist)
    nrparts = 0

    for idx, frag in enumerate(fraglist):
        if parts[idx] != 0:
            continue

        nrparts += 1
        parts[idx] = nrparts
        TraverseFragments(frag, fraglist, parts, nrparts)

    return (nrparts == 1)


def TraverseFragments(actfrag, fraglist, parts, nrparts):
    """
    This function was taken from openeye cookbook
    https://docs.eyesopen.com/toolkits/cookbook/python/cheminfo/enumfrags.html
    """
    for idx, frag in enumerate(fraglist):
        if parts[idx] != 0:
            continue

        if not IsAdjacentAtomBondSets(actfrag, frag):
            continue

        parts[idx] = nrparts
        TraverseFragments(frag, fraglist, parts, nrparts)


def CombineAndConnectAtomBondSets(fraglist):
    """
    This function was taken from OpeneEye Cookbook
    https://docs.eyesopen.com/toolkits/cookbook/python/cheminfo/enumfrags.html
    """

    # combine atom and bond sets

    combined = oechem.OEAtomBondSet()
    for frag in fraglist:
        for atom in frag.GetAtoms():
            combined.AddAtom(atom)
        for bond in frag.GetBonds():
            combined.AddBond(bond)

    # add connecting bonds

    for atomA in combined.GetAtoms():
        for atomB in combined.GetAtoms():
            if atomA.GetIdx() < atomB.GetIdx():
                continue

            bond = atomA.GetBond(atomB)
            if bond is None:
                continue
            if combined.HasBond(bond):
                continue

            combined.AddBond(bond)

    return combined


def GetFragmentAtomBondSetCombinations(fraglist, MAX_ROTORS=2, MIN_ROTORS=1):
    """
    This function was taken from OpeneEye cookbook
    https://docs.eyesopen.com/toolkits/cookbook/python/cheminfo/enumfrags.html
    Enumerate connected combinations from list of fragments
    Parameters
    ----------
    mol: OEMolGraph
    fraglist: list of OE AtomBondSet
    MAX_ROTORS: int
        min rotors in each fragment combination
    MIN_ROTORS: int
        max rotors in each fragment combination
    Returns
    -------
    fragcombs: list of connected combinations (OE AtomBondSet)
    """

    fragcombs = []

    nrfrags = len(fraglist)
    for n in range(1, nrfrags):

        for fragcomb in combinations(fraglist, n):

            if IsAdjacentAtomBondSetCombination(fragcomb):

                frag = CombineAndConnectAtomBondSets(fragcomb)

                if (CountRotorsInFragment(frag) <= MAX_ROTORS) and (CountRotorsInFragment(frag) >= MIN_ROTORS):
                    fragcombs.append(frag)

    return fragcombs