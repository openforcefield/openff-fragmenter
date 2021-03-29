import numpy as np

from fragmenter.utils import logger


def _find_torsions_from_smarts(molecule, smarts):
    """
    Do a substrcutre search on provided SMARTS to find torsions that match the SAMRTS

    Parameters
    ----------
    molecule: OEMol
        molecule to search on
    smarts: str
        SMARTS pattern to search for

    Returns
    -------
    tors: list
        list of torsions that match the SMARTS string

    """
    from openeye import oechem

    #ToDO use MDL aromaticity model
    qmol=oechem.OEQMol()
    if not oechem.OEParseSmarts(qmol, smarts):
        logger().warning('OEParseSmarts failed')
    ss = oechem.OESubSearch(qmol)
    tors = []
    oechem.OEPrepareSearch(molecule, ss)
    unique = True
    for match in ss.Match(molecule, unique):
        tor = []
        for ma in match.GetAtoms():
            tor.append(ma.target)
        tors.append(tor)

    return tors


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
        logger().debug("Map Idxs: {} {} {} {}".format(tor[0].GetMapIdx(), tor[1].GetMapIdx(), tor[2].GetMapIdx(), tor[3].GetMapIdx()))
        logger().debug("Atom Numbers: {} {} {} {}".format(tor[0].GetAtomicNum(), tor[1].GetAtomicNum(), tor[2].GetAtomicNum(), tor[3].GetAtomicNum()))
        if tor[1].GetMapIdx() != best_tor[1].GetMapIdx() or tor[2].GetMapIdx() != best_tor[2].GetMapIdx():
            #new_tor = True
            if not first_pass:
                logger().debug("Adding to list: {} {} {} {}".format(best_tor[0].GetMapIdx(), best_tor[1].GetMapIdx(), best_tor[2].GetMapIdx(), best_tor[3].GetMapIdx()))
                tors.append(best_tor)
            first_pass = False
            best_tor = tor
            best_tor_order = tor[0].GetAtomicNum() + tor[3].GetAtomicNum()
            logger().debug("new_tor with central bond across atoms: {} {}".format(tor[1].GetMapIdx(), tor[2].GetMapIdx()))
        else:
            logger().debug("Not a new_tor but now with end atoms: {} {}".format(tor[0].GetMapIdx(), tor[3].GetMapIdx()))
            tor_order = tor[0].GetAtomicNum() + tor[3].GetAtomicNum()
            if tor_order > best_tor_order:
                best_tor = tor
                best_tor_order = tor_order
    logger().debug("Adding to list: {} {} {} {}".format(best_tor[0].GetMapIdx(), best_tor[1].GetMapIdx(), best_tor[2].GetMapIdx(), best_tor[3].GetMapIdx()))
    tors.append(best_tor)

    logger().debug("List of torsion to drive:")
    for tor in tors:
        logger().debug("Idx: {} {} {} {}".format(tor[0].GetMapIdx(), tor[1].GetMapIdx(), tor[2].GetMapIdx(), tor[3].GetMapIdx()))
        logger().debug("Atom numbers: {} {} {} {}".format(tor[0].GetAtomicNum(), tor[1].GetAtomicNum(), tor[2].GetAtomicNum(), tor[3].GetAtomicNum()))

    return tors


def find_torsion_around_bond(molecule, bond):
    """
    Find the torsion around a given bond
    Parameters
    ----------
    molecule : molecule with atom maps
    bond : tuple of map idx of bond atoms

    Returns
    -------
    list of 4 atom map idx (-1)

    Note:
    This returns the map indices of the torsion -1, not the atom indices.

    """
    from openeye import oechem
    #if not has_atom_map(molecule):
    #    raise ValueError("Molecule must have atom maps")

    terminal_smarts = '[*]~[*]-[X2H1,X3H2,X4H3]-[#1]'
    terminal_torsions = _find_torsions_from_smarts(molecule, terminal_smarts)
    mid_torsions = [[tor.a, tor.b, tor.c, tor.d] for tor in oechem.OEGetTorsions(molecule)]
    all_torsions = terminal_torsions + mid_torsions

    tors = one_torsion_per_rotatable_bond(all_torsions)

    tor_idx = [tuple(i.GetMapIdx() for i in tor) for tor in tors]
    central_bonds = [(tor[1], tor[2]) for tor in tor_idx]
    try:
        idx = central_bonds.index(bond)
    except ValueError:
        idx = central_bonds.index(tuple(reversed(bond)))

    torsion = [i-1 for i in tor_idx[idx]]
    return torsion
