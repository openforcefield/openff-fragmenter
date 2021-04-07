"""functions to manipulate, read and write OpenEye and Psi4 molecules"""
import copy
from typing import Dict, List, Tuple

import cmiles
import networkx
import numpy
from openff.toolkit.topology import Molecule
from openff.toolkit.utils import LicenseError, ToolkitUnavailableException

from fragmenter.utils import to_off_molecule


def assign_elf10_am1_bond_orders(molecule, max_confs=800):
    """Generate ELF10 AM1 WBOs for an OpenEye OEMol molecule.

    Parameters
    ----------
    molecule : OEMol
        Molecule for which to generate conformers.
        Omega will be used to generate max_confs conformations.
    max_confs : int, optional, default=800
        Max number of conformers to generate

    Returns
    -------
        A molecule which contains ELF10 AM1 Wiberg Bond orders.
    """

    from openeye import oechem

    off_molecule = to_off_molecule(oechem.OEMol(molecule))

    # Store the atom map separately in case it gets removed by a TK.
    atom_map = off_molecule.properties.pop("atom_map", None)

    # Generate a set of ELF10 conformers.
    off_molecule = generate_conformers(off_molecule, max_confs)
    off_molecule.apply_elf_conformer_selection()

    per_conformer_bond_orders = []

    for conformer in off_molecule.conformers:

        off_molecule.assign_fractional_bond_orders(
            "am1-wiberg", use_conformers=[conformer]
        )

        per_conformer_bond_orders.append(
            [bond.fractional_bond_order for bond in off_molecule.bonds]
        )

    bond_orders = [*numpy.mean(per_conformer_bond_orders, axis=0)]

    for bond, bond_order in zip(off_molecule.bonds, bond_orders):
        bond.fractional_bond_order = bond_order

    # Restore the atom map.
    if atom_map is not None:
        off_molecule.properties["atom_map"] = atom_map

    # TODO: Return an OFF molecule object once the refactor is further along. Below
    #       this line is temporary code.
    oe_molecule = oechem.OEMol(molecule)

    if atom_map is not None:

        for oe_atom in oe_molecule.GetAtoms():
            oe_atom.SetMapIdx(atom_map[oe_atom.GetIdx()])

    for oe_bond in oe_molecule.GetBonds():

        bond = off_molecule.get_bond_between(oe_bond.GetBgnIdx(), oe_bond.GetEndIdx())
        oe_bond.SetData("WibergBondOrder", bond.fractional_bond_order)

    return oe_molecule


def generate_conformers(
    molecule: Molecule, max_confs: int = 800, rms_threshold: float = 1.0
) -> Molecule:
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
        A multi-conformer molecule with up to max_confs conformers.
    """

    from simtk import unit

    off_molecule = copy.deepcopy(molecule)

    # Store the atom map separately in case it gets removed / mangled by a TK.
    atom_map = off_molecule.properties.pop("atom_map", None)

    # Canonically order the atoms in the molecule before generating the conformer.
    # This helps ensure the same conformers are generated for the same molecules
    # independently of their atom order.
    canonical_molecule = off_molecule.canonical_order_atoms()

    canonical_molecule.generate_conformers(
        n_conformers=max_confs, rms_cutoff=rms_threshold * unit.angstrom
    )

    _, canonical_map = Molecule.are_isomorphic(
        canonical_molecule, off_molecule, return_atom_map=True
    )

    off_molecule = canonical_molecule.remap(canonical_map)

    # Restore the atom map.
    if atom_map is not None:
        off_molecule.properties["atom_map"] = atom_map

    return off_molecule


def find_ring_systems(molecule: Molecule) -> Dict[int, int]:
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
    ring_atom_index_pairs = {
        tuple(sorted(pair))
        for pair in molecule.chemical_environment_matches("[*:1]@[*:2]")
    }

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


def _find_oe_stereocenters(
    molecule: Molecule,
) -> Tuple[List[int], List[Tuple[int, int]]]:
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

    stereogenic_atoms = {
        atom.GetIdx() for atom in oe_molecule.GetAtoms() if atom.IsChiral()
    }

    stereogenic_bonds = {
        (bond.GetBgnIdx(), bond.GetEndIdx())
        for bond in oe_molecule.GetBonds()
        if bond.IsChiral()
    }

    return sorted(stereogenic_atoms), sorted(stereogenic_bonds)


def _find_rd_stereocenters(
    molecule: Molecule,
) -> Tuple[List[int], List[Tuple[int, int]]]:
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

    Chem.AssignStereochemistry(
        rd_molecule, cleanIt=True, force=True, flagPossibleStereoCenters=True
    )

    stereogenic_atoms = {
        atom.GetIdx()
        for atom in rd_molecule.GetAtoms()
        if atom.HasProp("_ChiralityPossible")
    }

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


def find_stereocenters(molecule: Molecule) -> Tuple[List[int], List[Tuple[int, int]]]:
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
    except (ModuleNotFoundError, ToolkitUnavailableException, LicenseError):
        stereogenic_atoms, stereogenic_bonds = _find_rd_stereocenters(molecule)

    return stereogenic_atoms, stereogenic_bonds


def smiles_to_oemol(smiles, add_atom_map: bool = False):
    """Create an OE molecule object from an input SMILES pattern.
    Parameters
    ----------
    smiles : str
        SMILES representation of desired molecule.
    add_atom_map
        Whether to create a canonical atom map for the molecule.

    Returns
    -------
    molecule : OEMol
        A normalized molecule with desired smiles string.
    """

    from openeye import oechem

    off_molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    oe_molecule = off_molecule.to_openeye()

    oechem.OEPerceiveChiral(oe_molecule)

    # Add canonical ordered atom maps
    if add_atom_map:
        cmiles.utils.add_atom_map(oe_molecule)

    return oe_molecule
