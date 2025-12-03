"""functions to manipulate, read and write OpenEye and Psi4 molecules"""

import copy
import logging

import networkx
import numpy
from openff.toolkit.topology import Molecule
from openff.toolkit.utils import LicenseError, ToolkitUnavailableException
from openff.utilities import MissingOptionalDependencyError, requires_oe_module

from openff.fragmenter.utils import get_atom_index, get_map_index

logger = logging.getLogger(__name__)


def assign_elf10_am1_bond_orders(molecule: Molecule, max_confs: int = 800, rms_threshold: float = 1.0) -> Molecule:
    """Generate ELF10 AM1 WBOs for a molecule.

    Parameters
    ----------
    molecule
        The molecule to compute WBOs for.
    max_confs : int, optional, default=800
        Max number of conformers to generate
    rms_threshold
        The minimum RMS value [Angstrom] at which two conformers are considered redundant
        and one is deleted.

    Returns
    -------
        A new molecule which contains ELF10 AM1 Wiberg Bond orders.
    """

    molecule = copy.deepcopy(molecule)

    # Store the atom map separately in case it gets removed by a TK.
    atom_map = molecule.properties.pop("atom_map", None)

    # Generate a set of ELF10 conformers.
    molecule = _generate_conformers(molecule, max_confs, rms_threshold)
    molecule.apply_elf_conformer_selection()

    per_conformer_bond_orders = []

    for conformer in molecule.conformers:
        molecule.assign_fractional_bond_orders("am1-wiberg", use_conformers=[conformer])

        per_conformer_bond_orders.append([bond.fractional_bond_order for bond in molecule.bonds])

    bond_orders = [*numpy.mean(per_conformer_bond_orders, axis=0)]

    for bond, bond_order in zip(molecule.bonds, bond_orders):
        bond.fractional_bond_order = bond_order

    # Restore the atom map.
    if atom_map is not None:
        molecule.properties["atom_map"] = atom_map

    return molecule


def _generate_conformers(molecule: Molecule, max_confs: int = 800, rms_threshold: float = 1.0) -> Molecule:
    """Generate conformations for the supplied molecule.

    Parameters
    ----------
    molecule
        Molecule for which to generate conformers
    max_confs
        Max number of conformers to generate.
    rms_threshold
        The minimum RMS value [Angstrom] at which two conformers are considered redundant
        and one is deleted.

    Returns
    -------
        A new multi-conformer molecule with up to max_confs conformers.
    """

    from openff.units import unit

    molecule = copy.deepcopy(molecule)

    # Store the atom map separately in case it gets removed / mangled by a TK.
    atom_map = molecule.properties.pop("atom_map", None)

    # Canonically order the atoms in the molecule before generating the conformer.
    # This helps ensure the same conformers are generated for the same molecules
    # independently of their atom order.
    canonical_molecule = molecule.canonical_order_atoms()

    canonical_molecule.generate_conformers(n_conformers=max_confs, rms_cutoff=rms_threshold * unit.angstrom)

    _, canonical_map = Molecule.are_isomorphic(canonical_molecule, molecule, return_atom_map=True)

    molecule = canonical_molecule.remap(canonical_map)

    # Restore the atom map.
    if atom_map is not None:
        molecule.properties["atom_map"] = atom_map

    return molecule


def find_ring_systems(molecule: Molecule) -> dict[int, int]:
    """This function attempts to find all ring systems (see [1] for more details) in
    a given molecule.

    The method first attempts to determine which atoms and bonds are part of rings
    by matching the `[*:1]@[*:2]` SMIRKS pattern against the molecule.

    The matched bonds are then used to construct a graph (using ``networkx``), from which
    the ring systems are identified as those sets of atoms which are 'connected'
    together (using ``connected_components``) by at least one path.

    Parameters
    ----------
    molecule:
        The molecule to search for ring systems.

    Notes
    -----
    * Two molecular rings with only one common atom (i.e. spiro compounds) are
      considered to be part of the same ring system.

    References
    ----------
    [1] `Ring Perception <https://docs.eyesopen.com/toolkits/python/oechemtk/ring.html>`_

    Returns
    -------
        The index of which ring system each atom belongs to. Only ring atoms are
        included in the returned dictionary.
    """

    # Find the ring atoms
    ring_atom_index_pairs = {tuple(sorted(pair)) for pair in molecule.chemical_environment_matches("[*:1]@[*:2]")}

    # Construct a networkx graph from the found ring bonds.
    graph = networkx.Graph()

    for atom_index_pair in ring_atom_index_pairs:
        graph.add_edge(*atom_index_pair)

    ring_systems = {}

    for i, ring_system in enumerate(networkx.connected_components(graph)):
        for atom_index in ring_system:
            assert atom_index not in ring_systems
            ring_systems[atom_index] = i + 1

    return ring_systems


@requires_oe_module("oechem")
def _find_oe_stereocenters(
    molecule: Molecule,
) -> tuple[list[int], list[tuple[int, int]]]:
    """A method which returns the of the stereogenic atom and bonds of a molecule
    using the OpenEye toolkit.

    Parameters
    ----------
    molecule
        The molecule whose stereocenters should be returned.

    Notes
    -----
    * This method currently deals with OE directly and should be removed once
      the API points suggested in openff-toolkit issue #903 have been added.

    Returns
    -------
        The indices of the stereogenic atoms in the molecule and the indices of the
        atoms involved in stereogenic bonds in the molecule.
    """

    from openeye import oechem

    oe_molecule = molecule.to_openeye()
    oechem.OEPerceiveChiral(oe_molecule)

    stereogenic_atoms = {atom.GetIdx() for atom in oe_molecule.GetAtoms() if atom.IsChiral()}

    stereogenic_bonds = {(bond.GetBgnIdx(), bond.GetEndIdx()) for bond in oe_molecule.GetBonds() if bond.IsChiral()}

    return sorted(stereogenic_atoms), sorted(stereogenic_bonds)


def _find_rd_stereocenters(
    molecule: Molecule,
) -> tuple[list[int], list[tuple[int, int]]]:
    """A method which returns the of the stereogenic atom and bonds of a molecule
    using the RDKit.

    Parameters
    ----------
    molecule
        The molecule whose stereocenters should be returned.

    Notes
    -----
    * This method currently deals with RDKit directly and should be removed once
      the API points suggested in openff-toolkit issue #903 have been added.

    Returns
    -------
        The indices of the stereogenic atoms in the molecule and the indices of the
        atoms involved in stereogenic bonds in the molecule.
    """

    from rdkit import Chem

    rd_molecule = molecule.to_rdkit()

    Chem.AssignStereochemistry(rd_molecule, cleanIt=True, force=True, flagPossibleStereoCenters=True)

    stereogenic_atoms = {atom.GetIdx() for atom in rd_molecule.GetAtoms() if atom.HasProp("_ChiralityPossible")}

    # Clear any previous assignments on the bonds, since FindPotentialStereo
    # may not overwrite it
    for bond in rd_molecule.GetBonds():
        bond.SetStereo(Chem.BondStereo.STEREONONE)

    Chem.FindPotentialStereoBonds(rd_molecule, cleanIt=True)

    stereogenic_bonds = {
        (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        for bond in rd_molecule.GetBonds()
        if bond.GetStereo() == Chem.BondStereo.STEREOANY
    }

    return sorted(stereogenic_atoms), sorted(stereogenic_bonds)


def find_stereocenters(molecule: Molecule) -> tuple[list[int], list[tuple[int, int]]]:
    """A method which returns the of the stereogenic atom and bonds of a molecule.

    Parameters
    ----------
    molecule
        The molecule whose stereocenters should be returned.

    Notes
    -----
    * This method currently deals with OE and RDKit molecules and should be updated once
      the API points suggested in openff-toolkit issue #903 have been added.

    Returns
    -------
        The indices of the stereogenic atoms in the molecule and the indices of the
        atoms involved in stereogenic bonds in the molecule.
    """

    try:
        stereogenic_atoms, stereogenic_bonds = _find_oe_stereocenters(molecule)
    except (
        ModuleNotFoundError,
        ToolkitUnavailableException,
        LicenseError,
        MissingOptionalDependencyError,
    ):
        stereogenic_atoms, stereogenic_bonds = _find_rd_stereocenters(molecule)

    return stereogenic_atoms, stereogenic_bonds


def _extract_rd_fragment(molecule: Molecule, atom_indices: set[int], bond_indices: set[tuple[int, int]]) -> Molecule:
    from rdkit import Chem

    rd_molecule = Chem.RWMol(molecule.to_rdkit())
    rd_atoms_by_map: dict[int, Chem.Atom] = {}

    # Restore the map indices as to_rdkit does not automatically add them.
    for atom in rd_molecule.GetAtoms():
        atom.SetAtomMapNum(get_map_index(molecule, atom.GetIdx()))

        rd_atoms_by_map[atom.GetAtomMapNum()] = atom

    atoms_to_use = [get_atom_index(molecule, i) for i in atom_indices]
    bonds_to_use = [
        rd_molecule.GetBondBetweenAtoms(get_atom_index(molecule, pair[0]), get_atom_index(molecule, pair[1])).GetIdx()
        for pair in bond_indices
    ]

    # Make sure to include any Hs bonded to the included atom set otherwise radicals
    # will form.
    for map_index in atom_indices:
        for neighbour in rd_atoms_by_map[map_index].GetNeighbors():
            if (
                neighbour.GetAtomicNum() != 1
                or neighbour.GetAtomMapNum() < 1
                or neighbour.GetAtomMapNum() in atom_indices
            ):
                continue

            atoms_to_use.append(neighbour.GetIdx())
            bonds_to_use.append(
                rd_molecule.GetBondBetweenAtoms(rd_atoms_by_map[map_index].GetIdx(), neighbour.GetIdx()).GetIdx()
            )

    # Add additional hydrogens to atoms where the total valence will change likewise to
    # ensure the valence does not change.
    rd_atoms_by_index = {atom.GetIdx(): atom for atom in rd_molecule.GetAtoms()}

    for atom_index in [*atoms_to_use]:
        atom = rd_atoms_by_index[atom_index]

        old_valence = atom.GetTotalValence()
        new_valence = atom.GetTotalValence()

        for neighbour_bond in rd_atoms_by_index[atom_index].GetBonds():
            if neighbour_bond.GetBeginAtomIdx() in atoms_to_use and neighbour_bond.GetEndAtomIdx() in atoms_to_use:
                continue

            new_valence -= neighbour_bond.GetValenceContrib(atom)

        if numpy.isclose(old_valence, new_valence):
            # Skip the cases where the valence won't change
            continue

        if (
            atom.GetAtomicNum() == 6
            and atom.GetIsAromatic()
            and sum(1 for bond_tuple in bond_indices if atom.GetAtomMapNum() in bond_tuple) == 1
        ):
            # This is likely a cap carbon which was retained from an existing ring. It's
            # aromaticity needs to be cleared before calling ``MolFragmentToSmiles``
            # otherwise will (understandably) be confused and throw an exception.
            atom.SetIsAromatic(False)

        # Add a hydrogen to the atom whose valence will change.
        for _ in range(int(numpy.rint(old_valence - new_valence))):
            new_atom = Chem.Atom(1)
            new_atom_index = rd_molecule.AddAtom(new_atom)

            rd_molecule.AddBond(atom_index, new_atom_index)

            new_bond = rd_molecule.GetBondBetweenAtoms(atom_index, new_atom_index)
            new_bond.SetBondType(Chem.BondType.SINGLE)
            new_bond.SetIsAromatic(False)

            atoms_to_use.append(new_atom_index)
            bonds_to_use.append(new_bond.GetIdx())

    fragment_smiles = Chem.MolFragmentToSmiles(rd_molecule, atoms_to_use, bonds_to_use)
    fragment = Molecule.from_smiles(fragment_smiles, allow_undefined_stereo=True)

    return fragment


@requires_oe_module("oechem")
def _extract_oe_fragment(molecule: Molecule, atom_indices: set[int], bond_indices: set[tuple[int, int]]) -> Molecule:
    from openeye import oechem

    oe_molecule = molecule.to_openeye()

    # Restore the map indices as to_openeye does not automatically add them.
    for atom_index, map_index in molecule.properties["atom_map"].items():
        oe_atom = oe_molecule.GetAtom(oechem.OEHasAtomIdx(atom_index))
        oe_atom.SetMapIdx(map_index)

    # Include any Hs bonded to the included atom set so we can retain their map
    # indices.
    for map_index in {*atom_indices}:
        oe_atom = oe_molecule.GetAtom(oechem.OEHasMapIdx(map_index))

        for neighbour in oe_atom.GetAtoms():
            if neighbour.GetAtomicNum() != 1 or neighbour.GetMapIdx() < 1 or neighbour.GetMapIdx() in atom_indices:
                continue

            atom_indices.add(neighbour.GetMapIdx())
            bond_indices.add((map_index, neighbour.GetMapIdx()))

    atom_bond_set = oechem.OEAtomBondSet()

    for map_index in atom_indices:
        atom = oe_molecule.GetAtom(oechem.OEHasMapIdx(map_index))
        atom_bond_set.AddAtom(atom)

    for map_index_1, map_index_2 in bond_indices:
        atom_1 = oe_molecule.GetAtom(oechem.OEHasMapIdx(map_index_1))
        atom_2 = oe_molecule.GetAtom(oechem.OEHasMapIdx(map_index_2))

        bond = oe_molecule.GetBond(atom_1, atom_2)

        if not bond:
            raise ValueError(f"{(map_index_1, map_index_2)} is a disconnected bond")

        atom_bond_set.AddBond(bond)

    atom_predicate = oechem.OEIsAtomMember(atom_bond_set.GetAtoms())
    bond_predicate = oechem.OEIsBondMember(atom_bond_set.GetBonds())

    fragment = oechem.OEMol()
    oechem.OESubsetMol(fragment, oe_molecule, atom_predicate, bond_predicate, True)

    oechem.OEAddExplicitHydrogens(fragment)
    oechem.OEPerceiveChiral(fragment)

    # Always restore map?
    # if restore_maps:
    # In some cases (symmetric molecules) this changes the atom map so skip it
    # restore_atom_map(fragment)
    # atom map should be restored for combinatorial fragmentation
    # Perceive stereo and check that defined stereo did not change
    oechem.OEPerceiveChiral(fragment)
    oechem.OE3DToAtomStereo(fragment)
    oechem.OE3DToBondStereo(fragment)

    return Molecule.from_openeye(fragment, allow_undefined_stereo=True)


def extract_fragment(molecule: Molecule, atom_indices: set[int], bond_indices: set[tuple[int, int]]) -> Molecule:
    """Returns a fragment which contains the specified atoms and bonds of a parent
    molecule

    Parameters
    ----------
    molecule
        The parent molecule being fragmented.
    atom_indices
        The *map* indices of the subset of atoms to retain.
    bond_indices
        The *map* indices of the subset of bonds to retain.
    """

    # Make sure that the bond and atom indices are self consistent and that there
    # are no disconnected bonds.
    if not all(i in atom_indices for map_tuple in bond_indices for i in map_tuple):
        raise ValueError("The ``bond_indices`` set includes atoms not in the ``atom_indices`` set.")

    try:
        fragment = _extract_oe_fragment(molecule, atom_indices, bond_indices)
    except (
        ModuleNotFoundError,
        ToolkitUnavailableException,
        LicenseError,
        MissingOptionalDependencyError,
    ):
        fragment = _extract_rd_fragment(molecule, atom_indices, bond_indices)

    # Sanity check that all atoms are still bonded
    fragment_smiles = fragment.to_smiles()

    assert "." not in fragment_smiles, (
        "An atom that is not bonded to any other atom in the fragment. "
        "You probably ran into a bug. Please report the input molecule to the "
        "issue tracker"
    )

    return fragment


def smiles_to_molecule(smiles, add_atom_map: bool = False) -> Molecule:
    """Create a molecule object from an input SMILES pattern.

    Parameters
    ----------
    smiles : str
        SMILES representation of desired molecule.
    add_atom_map
        Whether to create a canonical atom map for the molecule.

    Returns
    -------
        A normalized molecule with desired smiles string.
    """

    molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=True)

    # Add canonical ordered atom maps
    if add_atom_map:
        molecule = molecule.canonical_order_atoms()
        molecule.properties["atom_map"] = {i: i + 1 for i in range(molecule.n_atoms)}

    return molecule
