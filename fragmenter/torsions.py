__author__ = 'Chaya D. Stern'

try:
    import openeye.oechem as oechem
except ImportError:
    pass
import warnings
import numpy as np
import os
import json
import itertools

from . import utils, chemi
from cmiles import to_canonical_smiles_oe

warnings.simplefilter('always')


def find_torsions(molecule, restricted=True):
    """
    This function takes an OEMol (atoms must be tagged with index map) and finds the map indices for torsion that need
    to be driven.

    Parameters
    ----------
    molecule : OEMol
        The atoms in the molecule need to be tagged with map indices

    Returns
    -------
    needed_torsion_scans: dict
        a dictionary that maps internal and terminal torsions to map indices of torsion atoms

    """
    # Check if molecule has map
    is_mapped = chemi.is_mapped(molecule)
    if not is_mapped:
        utils.logger().warning('Molecule does not have atom map. A new map will be generated. You might need a new tagged SMARTS if the ordering was changed')
        tagged_smiles = to_canonical_smiles_oe(molecule, isomeric=True, mapped=True, explicit_hydrogen=True)
        # Generate new molecule with tags
        molecule = chemi.smiles_to_oemol(tagged_smiles)
        utils.logger().warning('If you already have a tagged SMARTS, compare it with the new one to ensure the ordering did not change')
        utils.logger().warning('The new tagged SMARTS is: {}'.format(tagged_smiles))
        # ToDo: save the new tagged SMILES somewhere. Maybe return it?

    if restricted:
        return _find_restricted_torsions(molecule)

    mid_tors = [[tor.a, tor.b, tor.c, tor.d ] for tor in oechem.OEGetTorsions(molecule)]

    # This smarts should match terminal torsions such as -CH3, -NH2, -NH3+, -OH, and -SH

    smarts = '[*]~[*]-[X2H1,X3H2,X4H3]-[#1]'
    qmol=oechem.OEQMol()
    if not oechem.OEParseSmarts(qmol, smarts):
        utils.logger().warning('OEParseSmarts failed')
    ss = oechem.OESubSearch(qmol)
    mol = oechem.OEMol(molecule)
    h_tors = []
    oechem.OEPrepareSearch(mol, ss)
    unique = True
    for match in ss.Match(mol, unique):
        tor = []
        for ma in match.GetAtoms():
            tor.append(ma.target)
        h_tors.append(tor)

    needed_torsion_scans = {'internal': {}, 'terminal': {}}
    if h_tors:
        h_tors_min = one_torsion_per_rotatable_bond(h_tors)
        for i, tor in enumerate(h_tors_min):
            tor_name = ((tor[0].GetMapIdx() -1 ), (tor[1].GetMapIdx() - 1), (tor[2].GetMapIdx() - 1), (tor[3].GetMapIdx() - 1))
            needed_torsion_scans['terminal']['torsion_{}'.format(str(i))] = tor_name
    if mid_tors:
        mid_tors_min = one_torsion_per_rotatable_bond(mid_tors)
        for i, tor in enumerate(mid_tors_min):
            tor_name = ((tor[0].GetMapIdx() - 1), (tor[1].GetMapIdx() - 1), (tor[2].GetMapIdx() - 1), (tor[3].GetMapIdx() - 1))
            needed_torsion_scans['internal']['torsion_{}'.format(str(i))] = tor_name

    # Check that there are no duplicate torsions in mid and h_torsions
    list_tor = list(needed_torsion_scans['internal'].values()) + list(needed_torsion_scans['terminal'].values())
    set_tor = set(list_tor)

    if not len(set_tor) == len(list_tor):
        raise Warning("There is a torsion defined in both mid and terminal torsions. This should not happen. Check "
                      "your molecule and the atom mapping")

    return needed_torsion_scans


def _find_restricted_torsions(molecule):

    restricted_smarts = '[*]~[C,c]=,@[C,c]~[*]' # This should capture double and triple bonds
    qmol=oechem.OEQMol()
    if not oechem.OEParseSmarts(qmol, restricted_smarts):
        utils.logger().warning('OEParseSmarts failed')
    ss = oechem.OESubSearch(qmol)
    mol = oechem.OEMol(molecule)
    restricted_tors = []
    oechem.OEPrepareSearch(mol, ss)
    unique = True
    for match in ss.Match(mol, unique):
        tor = []
        for ma in match.GetAtoms():
            tor.append(ma.target)
        restricted_tors.append(tor)

    restricted_torsion_scans = {'restricted_rotors': {}}
    if restricted_tors:
        non_rotor_min = one_torsion_per_rotatable_bond(restricted_tors)
        for i, tor in enumerate(non_rotor_min):
            tor_name = ((tor[0].GetMapIdx() - 1), (tor[1].GetMapIdx() - 1), (tor[2].GetMapIdx() - 1), (tor[3].GetMapIdx() - 1))
            restricted_torsion_scans['non_rotors']['torsion_{}'.format(str(i))] = tor_name
    return restricted_torsion_scans


def one_torsion_per_rotatable_bond(torsion_list):
    """
    Keep only one torsion per rotatable bond
    Parameters
    ----------
    torsion_list: list
        list of torsion in molecule

    Returns
    -------
    list of only one torsion per rotatable bonds

    """

    central_bonds = np.zeros((len(torsion_list), 3), dtype=int)
    for i, tor in enumerate(torsion_list):
        central_bonds[i][0] = i
        central_bonds[i][1] = tor[1].GetIdx()
        central_bonds[i][2] = tor[2].GetIdx()

    grouped = central_bonds[central_bonds[:, 2].argsort()]
    sorted_tors = [torsion_list[i] for i in grouped[:, 0]]

    # Keep only one torsion per rotatable bond
    tors = []
    best_tor = [sorted_tors[0][0], sorted_tors[0][0], sorted_tors[0][0], sorted_tors[0][0]]
    best_tor_order = best_tor[0].GetAtomicNum() + best_tor[3].GetAtomicNum()
    first_pass = True
    for tor in sorted_tors:
        utils.logger().debug("Map Idxs: {} {} {} {}".format(tor[0].GetMapIdx(), tor[1].GetMapIdx(), tor[2].GetMapIdx(), tor[3].GetMapIdx()))
        utils.logger().debug("Atom Numbers: {} {} {} {}".format(tor[0].GetAtomicNum(), tor[1].GetAtomicNum(), tor[2].GetAtomicNum(), tor[3].GetAtomicNum()))
        if tor[1].GetMapIdx() != best_tor[1].GetMapIdx() or tor[2].GetMapIdx() != best_tor[2].GetMapIdx():
            new_tor = True
            if not first_pass:
                utils.logger().debug("Adding to list: {} {} {} {}".format(best_tor[0].GetMapIdx(), best_tor[1].GetMapIdx(), best_tor[2].GetMapIdx(), best_tor[3].GetMapIdx()))
                tors.append(best_tor)
            first_pass = False
            best_tor = tor
            best_tor_order = tor[0].GetAtomicNum() + tor[3].GetAtomicNum()
            utils.logger().debug("new_tor with central bond across atoms: {} {}".format(tor[1].GetMapIdx(), tor[2].GetMapIdx()))
        else:
            utils.logger().debug("Not a new_tor but now with end atoms: {} {}".format(tor[0].GetMapIdx(), tor[3].GetMapIdx()))
            tor_order = tor[0].GetAtomicNum() + tor[3].GetAtomicNum()
            if tor_order > best_tor_order:
                best_tor = tor
                best_tor_order = tor_order
    utils.logger().debug("Adding to list: {} {} {} {}".format(best_tor[0].GetMapIdx(), best_tor[1].GetMapIdx(), best_tor[2].GetMapIdx(), best_tor[3].GetMapIdx()))
    tors.append(best_tor)

    utils.logger().info("List of torsion to drive:")
    for tor in tors:
        utils.logger().info("Idx: {} {} {} {}".format(tor[0].GetMapIdx(), tor[1].GetMapIdx(), tor[2].GetMapIdx(), tor[3].GetMapIdx()))
        utils.logger().info("Atom numbers: {} {} {} {}".format(tor[0].GetAtomicNum(), tor[1].GetAtomicNum(), tor[2].GetAtomicNum(), tor[3].GetAtomicNum()))

    return tors


def define_torsiondrive_jobs(needed_torsion_drives, internal_torsion_resolution=30, terminal_torsion_resolution=0,
                     scan_internal_terminal_combination=0, scan_dimension=2):
    """
    define crank jobs with torsions to drive and resolution to drive them at.

    Parameters
    ----------
    fragment_data: dict
        dictionary that maps fragment to needed torsions
    internal_torsion_resolution: int, optional. Default 15
        interval in degrees for torsion scan. If 0, internal torsions will not be driven
    terminal_torsion_resolution: int, optional. Default 0
        interval in degrees for torsion scans. If 0, terminal torsions will not be driven
    scan_internal_terminal_combination: int, optional. Default 0
        flag if internal and terminal torsions should be combined for higher dimension. If 0, only internal torsions will
        be driven. If 1, terminal and internal torsions will be scanned together.
    scan_dimension: int, optional. Default 2
        dimension of torsion scan. Combinations of torsions at the specified dimension will be generated as separate crank jobs
    qc_program: str, optional. Default Psi4
    method: str, optional. Default B3LYP
    basis: str, optional. Default aug-cc-pVDZ
    kwargs: optional keywords for psi4

    Returns
    -------
    fragment_data: dict
        dictionary that maps fragment to crank torsion jobs specifications.

    """

    if not internal_torsion_resolution and not terminal_torsion_resolution:
        utils.logger().warning("Resolution for internal and terminal torsions are 0. No torsions will be driven", Warning)

    if scan_internal_terminal_combination and (not internal_torsion_resolution or not terminal_torsion_resolution):
        raise Warning("If you are not scanning internal or terminal torsions, you must set scan_internal_terminal_"
                      "combinations to 0")

    internal_torsions = needed_torsion_drives['internal']
    terminal_torsions = needed_torsion_drives['terminal']
    internal_dimension = len(internal_torsions)
    terminal_dimension = len(terminal_torsions)
    torsion_dimension = internal_dimension + terminal_dimension

    crank_job = 0
    crank_jobs = dict()

    if not scan_internal_terminal_combination:
        if internal_torsion_resolution:
            for comb in itertools.combinations(internal_torsions, scan_dimension):
                dihedrals = [internal_torsions[torsion] for torsion in comb]
                grid = [internal_torsion_resolution]*len(dihedrals)
                crank_jobs['crank_job_{}'.format(crank_job)] = {'dihedrals': dihedrals, 'grid_spacing': grid}
                crank_job +=1
            if internal_dimension < scan_dimension and internal_dimension > 0:
                dihedrals = [internal_torsions[torsion] for torsion in internal_torsions]
                grid = [internal_torsion_resolution]*len(dihedrals)
                crank_jobs['crank_job_{}'.format(crank_job)] = {'dihedrals': dihedrals, 'grid_spacing': grid}
                crank_job +=1

        if terminal_torsion_resolution:
            for comb in itertools.combinations(terminal_torsions, scan_dimension):
                dihedrals = [terminal_torsions[torsion] for torsion in comb]
                grid = [terminal_torsion_resolution]*scan_dimension
                crank_jobs['crank_job_{}'.format(crank_job)] = {'dihedrals': dihedrals, 'grid_spacing': grid}
                crank_job +=1
            if terminal_dimension < scan_dimension and terminal_dimension > 0:
                dihedrals = [terminal_torsions[torsion] for torsion in terminal_torsions]
                grid = [terminal_torsion_resolution]*len(dihedrals)
                crank_jobs['crank_job_{}'.format(crank_job)] = {'dihedrals': dihedrals, 'grid_spacing': grid}
                crank_job +=1
    else:
        # combine both internal and terminal torsions
        all_torsion_idx = np.arange(0, torsion_dimension)
        for comb in itertools.combinations(all_torsion_idx, scan_dimension):
            internal_torsions = [internal_torsions['torsion_{}'.format(i)] for i in comb if i < internal_dimension]
            terminal_torsions = [terminal_torsions['torsion_{}'.format(i-internal_dimension)] for i in comb if i >= internal_dimension]
            grid = [internal_torsion_resolution]*len(internal_torsions)
            grid.extend([terminal_torsion_resolution]*len(terminal_torsions))
            dihedrals = internal_torsions + terminal_torsions
            crank_jobs['crank_job_{}'.format(crank_job)] = {'diherals': dihedrals, 'grid_spacing': grid}
            crank_job += 1
        if torsion_dimension < scan_dimension:
            internal_torsions = [internal_torsions['torsion_{}'.format(i)] for i in all_torsion_idx if i < internal_dimension]
            terminal_torsions = [terminal_torsions['torsion_{}'.format(i-internal_dimension)] for i in all_torsion_idx if i >= internal_dimension]
            grid = [internal_torsion_resolution]*len(internal_torsions)
            grid.extend([terminal_torsion_resolution]*len(terminal_torsions))
            dihedrals = internal_torsions + terminal_torsions
            crank_jobs['crank_job_{}'.format(crank_job)] = {'diherals': dihedrals, 'grid_spacing': grid}
            crank_job += 1

    return crank_jobs


def get_initial_crank_state(fragment):
    """
    Generate initial crank state JSON for each crank job in fragment
    Parameters
    ----------
    fragment: dict
        A fragment from JSON crank jobs

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
        needed_mid_torsions = needed_torsions['internal']
        for mid_torsion in crank_jobs[job]['internal_torsions']:
            # convert 1-based indexing to 0-based indexing
            dihedrals.append([j-1 for j in needed_mid_torsions[mid_torsion]])
            grid_spacing.append(crank_jobs[job]['internal_torsions'][mid_torsion])
        needed_terminal_torsions = needed_torsions['terminal']
        for terminal_torsion in crank_jobs[job]['terminal_torsions']:
            # convert 1-based indexing to 0-based indexing
            dihedrals.append([j-1 for j in needed_terminal_torsions[terminal_torsion]])
            grid_spacing.append(crank_jobs[job]['terminal_torsions'][terminal_torsion])

        crank_state = {}
        crank_state['dihedrals'] = dihedrals
        crank_state['grid_spacing'] = grid_spacing
        crank_state['elements'] = fragment['molecule']['symbols']

        #ToDo add ability to start with many geomotries
        crank_state['init_coords'] = [init_geometry]
        crank_state['grid_status'] = {}

        crank_initial_states[job] = crank_state
    return crank_initial_states
