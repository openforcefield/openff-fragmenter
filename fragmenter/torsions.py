import logging
from typing import Tuple

from openff.toolkit.topology import Molecule

from fragmenter.utils import get_map_index

logger = logging.getLogger(__name__)


def find_torsion_around_bond(molecule: Molecule, bond: Tuple[int, int]) -> Tuple[int]:
    """Find the torsion around a given central bond. When multiple torsions are found,
    the torsion with the heaviest end groups (i.e. with the largest mass_0 + mass_3)
    is returned.

    Parameters
    ----------
    molecule
        The molecule containing the central bond.
    bond
        The *map* indices of the atoms involved in the central bond.

    Returns
    -------
        The map indices of the torsion -1, *not* the atom indices.
    """

    # Sort the bond indices to make matching torsions easier.
    central_bond = sorted(bond)

    largest_rank = -1
    found_torsion = None

    for atoms in molecule.propers:

        atom_indices = tuple(atom.molecule_atom_index for atom in atoms)

        try:
            map_indices = tuple(get_map_index(molecule, i) for i in atom_indices)
        except KeyError:
            # Skip torsions involving un-mapped (i.e. cap) atoms.
            continue

        if sorted(map_indices[1:3]) != central_bond:
            continue

        rank = sum(molecule.atoms[atom_indices[i]].atomic_number for i in [0, 3])

        if rank <= largest_rank:
            continue

        largest_rank = rank
        found_torsion = [i - 1 for i in map_indices]

    assert found_torsion is not None, "map indices should always be found."

    return tuple(found_torsion)
