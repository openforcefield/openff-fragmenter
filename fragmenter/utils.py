import os
from typing import Dict

import yaml
from openff.toolkit.topology import Molecule
from pkg_resources import resource_filename


def get_fgroup_smarts() -> Dict[str, str]:
    """Returns a dictionary containing the SMARTS representations of different
    functional groups loaded from the internal ``fgroup_smarts.yml`` file.

    Returns
    -------
        A dictionary where each key is the name of a functional group and each value
        the corresponding SMARTS pattern.
    """

    file_name = resource_filename(
        "fragmenter", os.path.join("data", "fgroup_smarts.yml")
    )

    with open(file_name, "r") as file:
        functional_groups = yaml.safe_load(file)

    return functional_groups


def get_map_index(
    molecule: Molecule, atom_index: int, error_on_missing: bool = True
) -> int:
    """Returns the map index of a particular atom in a molecule.

    Parameters
    ----------
    molecule
        The molecule containing the atom.
    atom_index
        The index of the atom in the molecule.
    error_on_missing
        Whether an error should be raised if the atom does not have a corresponding
        map index

    Returns
    -------
        The map index if found, otherwise 0.
    """
    atom_map = molecule.properties.get("atom_map", {})
    atom_map_index = atom_map.get(atom_index, None)

    if atom_map_index is None and error_on_missing:
        raise KeyError(f"{atom_index} is not in the atom map ({atom_map}).")

    return 0 if atom_map_index is None else atom_map_index


def get_atom_index(molecule: Molecule, map_index: int) -> int:
    """Returns the atom index of the atom in a molecule which has the specified map
    index.

    Parameters
    ----------
    molecule
        The molecule containing the atom.
    map_index
        The map index of the atom in the molecule.

    Returns
    -------
        The corresponding atom index
    """
    inverse_atom_map = {
        j: i for i, j in molecule.properties.get("atom_map", {}).items()
    }

    atom_index = inverse_atom_map.get(map_index, None)

    assert (
        atom_index is not None
    ), f"{map_index} does not correspond to an atom in the molecule."

    return atom_index
