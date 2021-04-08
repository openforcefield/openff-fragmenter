import abc
import logging
import time
import warnings
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple, Union

import cmiles
from cmiles.utils import (
    add_atom_map,
    has_explicit_hydrogen,
    has_stereo_defined,
    is_missing_atom_map,
    mol_to_map_ordered_qcschema,
    mol_to_smiles,
    restore_atom_map,
    to_canonical_label,
)
from openff.toolkit.topology import Atom, Molecule

from fragmenter import torsions

from .chemi import (
    assign_elf10_am1_bond_orders,
    find_ring_systems,
    find_stereocenters,
    generate_conformers,
)
from .states import _enumerate_stereoisomers
from .utils import get_fgroup_smarts, get_map_index, to_off_molecule

logger = logging.getLogger(__name__)

BondTuple = Tuple[int, int]
AtomAndBondSet = Tuple[Set[int], Set[BondTuple]]


class Fragmenter(abc.ABC):
    """Base fragmenter class."""

    def __init__(self, molecule, functional_groups):
        from openeye import oechem

        self.molecule = molecule

        # Add explicit hydrogen
        if not has_explicit_hydrogen(self.molecule):
            oechem.OEAddExplicitHydrogens(self.molecule)

        # Add canonical atom map
        add_atom_map(self.molecule)

        if is_missing_atom_map(self.molecule):
            warnings.warn(
                "Some atoms are missing atom maps. This might cause a problem at several points during "
                "the fragmentation process. Make sure you know what you are doing. "
            )

        # Keep track of stereo to make sure it does not flip
        self.stereo = {}
        self._atom_stereo_map = {
            oechem.OECIPAtomStereo_S: "S",
            oechem.OECIPAtomStereo_R: "R",
            oechem.OECIPAtomStereo_NotStereo: None,
            oechem.OECIPAtomStereo_UnspecStereo: "unspecified",
        }
        self._bond_stereo_map = {
            oechem.OECIPBondStereo_E: "E",
            oechem.OECIPBondStereo_Z: "Z",
        }
        self._find_stereo()

        # keep track of fragments that form new stereocenters
        self.new_stereo = []

        # For provenance
        self._options = {}

        # Find the functional groups.
        self.functional_groups = {}

        self._options["functional_groups"] = functional_groups
        self._tag_functional_groups(functional_groups)

        # Track any ring systems.
        self.ring_systems = {}

        # Fragments from fragmentation scheme for each rotatable bond
        self.fragments = {}

    def _find_stereo(self):
        """Find chiral atoms and bonds, store the chirality. This is needed to check if
        fragments flipped chirality. Currently this can happen and it is a bug
        """

        off_molecule = to_off_molecule(self.molecule)

        atom_stereo = {
            get_map_index(off_molecule, atom.molecule_atom_index): atom.stereochemistry
            for atom in off_molecule.atoms
            if atom.stereochemistry is not None
        }

        bond_stereo = {
            (
                get_map_index(off_molecule, bond.atom1_index),
                get_map_index(off_molecule, bond.atom2_index),
            ): bond.stereochemistry
            for bond in off_molecule.bonds
            if bond.stereochemistry is not None
        }

        self.stereo = {**atom_stereo, **bond_stereo}

    def _check_stereo(self, fragment) -> bool:
        """Check if stereo in fragment is different than stereo in parent.

        Parameters
        ----------
        fragment : oemol
        """

        off_fragment = to_off_molecule(fragment)

        atom_stereocenters, bond_stereocenters = find_stereocenters(off_fragment)

        # Check for new / flipped chiral centers.
        for atom_index in atom_stereocenters:

            map_index = get_map_index(off_fragment, atom_index)

            if map_index not in self.stereo:

                logger.warning(f"A new stereocenter formed at atom {map_index}")
                return False

            fragment_stereo = off_fragment.atoms[atom_index].stereochemistry
            parent_stereo = self.stereo[map_index]

            if fragment_stereo != parent_stereo:

                logger.warning(
                    f"Stereochemistry for atom {map_index} flipped from "
                    f"{parent_stereo} to {fragment_stereo}"
                )

                return False

        for index_tuple in bond_stereocenters:

            map_tuple = tuple(get_map_index(off_fragment, i) for i in index_tuple)

            map_tuple = (
                map_tuple if map_tuple in self.stereo else tuple(reversed(map_tuple))
            )

            if map_tuple not in self.stereo:

                logger.warning(f"A new chiral bond formed at bond {map_tuple}")
                return False

            fragment_stereo = off_fragment.get_bond_between(
                *index_tuple
            ).stereochemistry

            parent_stereo = self.stereo[map_tuple]

            if fragment_stereo != parent_stereo:

                logger.warning(
                    f"Stereochemistry for bond {map_tuple} flipped from "
                    f"{parent_stereo} to {fragment_stereo}"
                )

                return False

        return True

    def _fix_stereo(self, fragment):
        """Flip all stereocenters and find the stereoisomer that matches the parent

        Parameters
        ----------
        fragment : oemol

        Returns
        -------
            oemol with same stereochemistry as parent molecule
        """

        from openeye import oechem

        for stereoisomer in _enumerate_stereoisomers(fragment, 200, True, True):

            if not self._check_stereo(stereoisomer):
                continue

            return stereoisomer

        raise RuntimeError(
            f"The stereochemistry of {oechem.OEMolToSmiles(fragment)} could not be "
            f"fixed."
        )

    def _tag_functional_groups(
        self, functional_groups: Optional[Union[bool, Dict[str, str]]]
    ):
        """Tags atoms and bonds of functional groups specified by ``functional_groups``.

        Parameters
        ----------
        functional_groups
            A dictionary of SMARTS of functional groups that should not be fragmented.
            If it is None, ``fragmenter/fragmenter/data/fgroup_smarts.yaml`` will be
            used. If False, no functional groups will be tagged and they will all be
            fragmented.
        """

        if functional_groups is False:
            return

        if functional_groups is None:
            functional_groups = get_fgroup_smarts()

        off_molecule = to_off_molecule(self.molecule)

        for functional_group, smarts in functional_groups.items():

            unique_matches = {
                tuple(sorted(match))
                for match in off_molecule.chemical_environment_matches(smarts)
            }

            for i, match in enumerate(unique_matches):

                atoms = set(get_map_index(off_molecule, index) for index in match)
                bonds = set(
                    (
                        get_map_index(off_molecule, bond.atom1_index),
                        get_map_index(off_molecule, bond.atom2_index),
                    )
                    for bond in off_molecule.bonds
                    if bond.atom1_index in match and bond.atom2_index in match
                )

                self.functional_groups[f"{functional_group}_{i}"] = (atoms, bonds)

    def _find_rotatable_bonds(self) -> List[BondTuple]:
        """Finds the rotatable bonds in a molecule *including* rotatable double
        bonds.

        Notes
        -----
        * This does not find terminal rotatable bonds such as -OH, -NH2 -CH3.

        Todos
        -----
        * Add the option to build fragments around terminal torsions (-OH, -NH2, -CH3)

        Returns
        -------
        rotatable_bonds: list of tuples
            list of rotatable bonds map indices [(m1, m2),...]

        """

        off_molecule = to_off_molecule(self.molecule)

        matches = off_molecule.chemical_environment_matches(
            "[!$(*#*)&!D1:1]-,=;!@[!$(*#*)&!D1:2]"
        )
        unique_matches = {tuple(sorted(match)) for match in matches}

        # Drop bonds without a heavy degree of at least 2 on each end to avoid finding
        # terminal bonds
        def heavy_degree(atom_index: int) -> int:
            atom = off_molecule.atoms[atom_index]
            return sum(1 for atom in atom.bonded_atoms if atom.atomic_number != 1)

        unique_matches = {
            match for match in unique_matches if all(heavy_degree(i) > 1 for i in match)
        }

        return [
            (
                get_map_index(off_molecule, match[0]),
                get_map_index(off_molecule, match[1]),
            )
            for match in unique_matches
        ]

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
        from openeye import oechem

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

        return atom_bond_set

    def _atom_bond_set_to_mol(self, atoms, bonds, adjust_hcount=True):
        """
        Convert fragments (AtomBondSet) to OEMol
        Parameters
        ----------
        atoms : set of ints
            set of map indices
        bonds : set of tuples of ints
            set of bond tuples (m1, m2)
        adjust_hcount: bool, optional, default True
            If False, hydrogen counts will not be adjusted. Not recommended.

        Returns
        -------
        fragment: OEMol
        """
        from openeye import oechem

        frag = self._to_atom_bond_set(atoms, bonds)

        fragatompred = oechem.OEIsAtomMember(frag.GetAtoms())
        fragbondpred = oechem.OEIsBondMember(frag.GetBonds())

        fragment = oechem.OEMol()
        adjustHCount = adjust_hcount
        oechem.OESubsetMol(
            fragment, self.molecule, fragatompred, fragbondpred, adjustHCount
        )

        oechem.OEAddExplicitHydrogens(fragment)
        oechem.OEPerceiveChiral(fragment)
        # sanity check that all atoms are bonded
        for atom in fragment.GetAtoms():
            if not list(atom.GetBonds()):
                warnings.warn(
                    "Yikes!!! An atom that is not bonded to any other atom in the fragment. "
                    "You probably ran into a bug. Please report the input molecule to the issue tracker"
                )
        # Always restore map?
        # if restore_maps:
        # In some cases (symmetric molecules) this changes the atom map so skip it
        # restore_atom_map(fragment)
        # atom map should be restored for combinatorial fragmentation
        # Perceive stereo and check that defined stereo did not change
        oechem.OEPerceiveChiral(fragment)
        oechem.OE3DToAtomStereo(fragment)
        oechem.OE3DToBondStereo(fragment)
        if has_stereo_defined(fragment):
            if not self._check_stereo(fragment):
                fragment = self._fix_stereo(fragment)

        return fragment

    def _get_torsion_quartet(self, bond):
        """
        Get all atoms bonded to the torsion quartet around rotatable bond

        Parameters
        ----------
        bond: tuple of ints
            map indices of atoms in bond

        Returns
        -------
        atoms, bonds: set of ints (map indices of atoms in quartet), set of tuples of ints (bonds in quartet)

        """
        from openeye import oechem

        atom_map_idx = set()
        bond_tuples = set()
        atoms = [self.molecule.GetAtom(oechem.OEHasMapIdx(i)) for i in bond]

        m1, m2 = atoms[0].GetMapIdx(), atoms[1].GetMapIdx()
        atom_map_idx.update((m1, m2))
        bond_tuples.add((m1, m2))
        for atom in atoms:
            for a in atom.GetAtoms():
                m = a.GetMapIdx()
                atom_map_idx.add(m)
                bond_tuples.add((atom.GetMapIdx(), m))
                for nbr in a.GetAtoms():
                    m_nbr = nbr.GetMapIdx()
                    atom_map_idx.add(m_nbr)
                    bond_tuples.add((m, m_nbr))

        return atom_map_idx, bond_tuples

    def _find_ring_systems(self, keep_non_rotor_ring_substituents: bool = False):
        """This function finds all ring systems and stores them in the `ring_systems`
        field.

        Parameters
        ----------
        keep_non_rotor_ring_substituents
            If True, keep all non rotatable ring substituents. According to the
            benchmark, it is not necessary.
        """

        off_molecule = to_off_molecule(self.molecule)

        atom_to_ring_indices = find_ring_systems(off_molecule)

        # Find the map indices of the atoms involved in each ring system.
        ring_system_atoms = {
            ring_index: {
                get_map_index(off_molecule, i)
                for i in atom_to_ring_indices
                if atom_to_ring_indices[i] == ring_index
            }
            for ring_index in {*atom_to_ring_indices.values()}
        }

        # Find the map indices of the bonds involved in each ring system.
        ring_system_bonds = defaultdict(set)

        for bond in off_molecule.bonds:

            ring_index_1 = atom_to_ring_indices.get(bond.atom1_index, -1)
            ring_index_2 = atom_to_ring_indices.get(bond.atom2_index, -2)

            if ring_index_1 != ring_index_2:
                continue

            ring_system_bonds[ring_index_1].add(
                (
                    get_map_index(off_molecule, bond.atom1_index),
                    get_map_index(off_molecule, bond.atom2_index),
                )
            )

        # Scan the neighbours of the ring system atoms for any functional groups
        # / non-rotor substituents which should be included in the ring systems.
        map_index_to_functional_group = {
            map_index: f_group
            for f_group in self.functional_groups
            for map_index in self.functional_groups[f_group][0]
        }

        for ring_index in ring_system_atoms:

            # If any atoms are part of a functional group, include the other atoms in the
            # group in the ring system lists
            functional_groups = {
                map_index_to_functional_group[map_index]
                for map_index in ring_system_atoms[ring_index]
                if map_index in map_index_to_functional_group
            }

            ring_system_atoms[ring_index].update(
                map_index
                for functional_group in functional_groups
                for map_index in self.functional_groups[functional_group][0]
            )
            ring_system_bonds[ring_index].update(
                map_tuple
                for functional_group in functional_groups
                for map_tuple in self.functional_groups[functional_group][1]
            )

            if not keep_non_rotor_ring_substituents:
                continue

            non_rotor_atoms, non_rotor_bonds = self._find_non_rotor_ring_substituents(
                off_molecule, ring_index, ring_system_atoms
            )

            ring_system_atoms[ring_index].update(non_rotor_atoms)
            ring_system_bonds[ring_index].update(non_rotor_bonds)

        for ring_index in ring_system_atoms:

            self.ring_systems[ring_index] = (
                ring_system_atoms[ring_index],
                ring_system_bonds[ring_index],
            )

    def _find_non_rotor_ring_substituents(
        self, off_molecule, ring_index, ring_system_atoms
    ):
        """Find the non-rotor substituents attached to a particular ring system."""

        rotatable_bonds = off_molecule.find_rotatable_bonds()

        def heavy_degree(atom: Atom) -> int:
            return sum(1 for atom in atom.bonded_atoms if atom.atomic_number != 1)

        rotor_bonds = [
            bond
            for bond in rotatable_bonds
            if heavy_degree(bond.atom1) >= 2 and heavy_degree(bond.atom2) >= 2
        ]

        non_rotor_atoms = set()
        non_rotor_bonds = set()

        for bond in off_molecule.bonds:

            # Check if the bond is a rotor.
            if bond in rotor_bonds:
                continue

            if bond.atom1.atomic_number == 1 or bond.atom2.atomic_number == 1:
                continue

            map_index_1 = get_map_index(off_molecule, bond.atom1_index)
            map_index_2 = get_map_index(off_molecule, bond.atom2_index)

            in_system_1 = map_index_1 in ring_system_atoms[ring_index]
            in_system_2 = map_index_2 in ring_system_atoms[ring_index]

            if (in_system_1 and in_system_2) or (not in_system_1 and not in_system_2):
                continue

            non_rotor_atoms.update((map_index_1, map_index_2))
            non_rotor_bonds.add((map_index_1, map_index_2))

        return non_rotor_atoms, non_rotor_bonds

    def _get_ring_and_fgroups(
        self, atoms: Set[int], bonds: Set[BondTuple]
    ) -> AtomAndBondSet:
        """Keep ortho substituents

        Parameters
        ----------
        atoms
            map indices of atom in fragment
        bonds
            map indices of bonds in fragment

        Returns
        -------
            The updated set of atom and bond map indices to to retain.
        """

        off_molecule = to_off_molecule(molecule=self.molecule)

        # Find the sets of atoms which are located ortho to one of the bonds being
        # fragmented.
        ortho_atoms, ortho_bonds = self._find_ortho_substituents(off_molecule, bonds)

        # Include the rings systems and functional groups connected to the current
        # atom sets.
        atoms.update(ortho_atoms)
        bonds.update(ortho_bonds)

        map_index_to_functional_group = {
            map_index: functional_group
            for functional_group in self.functional_groups
            for map_index in self.functional_groups[functional_group][0]
        }
        map_index_to_ring_system = {
            map_index: ring_system
            for ring_system in self.ring_systems
            for map_index in self.ring_systems[ring_system][0]
        }

        new_atoms = set()
        new_bonds = set()

        for map_index in atoms:

            if map_index in map_index_to_functional_group:

                functional_group = map_index_to_functional_group[map_index]

                new_atoms.update(self.functional_groups[functional_group][0])
                new_bonds.update(self.functional_groups[functional_group][1])

            if map_index in map_index_to_ring_system:

                ring_index = map_index_to_ring_system[map_index]

                new_atoms.update(self.ring_systems[ring_index][0])
                new_bonds.update(self.ring_systems[ring_index][1])

        atoms.update(new_atoms)
        bonds.update(new_bonds)

        # Ensure the matched bonds doesn't include duplicates.
        bonds = {tuple(sorted(bond)) for bond in bonds}

        return atoms, bonds

    def _find_ortho_substituents(
        self, molecule: Molecule, bonds: Set[BondTuple]
    ) -> AtomAndBondSet:
        """Find ring substituents that are ortho to one of the rotatable bonds specified
        in a list of bonds.

        Parameters
        ----------
        molecule
            The molecule being fragmented.
        bonds
            The map indices of the rotatable bonds.

        Returns
        -------
            The set of map indices of atoms in ortho group and of bond tuples in ortho
            group.
        """

        matched_atoms = set()
        matched_bonds = set()

        for match in molecule.chemical_environment_matches(
            "[!#1:1]~&!@[*:2]@[*:3]~&!@[!#1*:4]"
        ):

            map_tuple = tuple(get_map_index(molecule, i) for i in match)

            if map_tuple[:2] not in bonds and map_tuple[:2][::-1] not in bonds:
                continue

            matched_atoms.update(map_tuple[::3])
            matched_bonds.update((map_tuple[i], map_tuple[i + 1]) for i in [0, 2])

        # Ensure the matched bonds doesn't include duplicates.
        matched_bonds = {tuple(sorted(bond)) for bond in matched_bonds}

        return matched_atoms, matched_bonds

    def _cap_open_valence(self, atoms, bonds, target_bond):
        """
        Cap with methyl for fragments that ends with N, O or S. Otherwise cap with H
        Parameters
        ----------
        atoms: set of ints
            map indices of atom in fragment
        bpnds: set of tuples of ints
            map indices of bonds in fragment
        target_bond: tuple of ints
            bond map indices the fragment is for
        """
        from openeye import oechem

        map_index_to_functional_group = {
            map_index: functional_group
            for functional_group in self.functional_groups
            for map_index in self.functional_groups[functional_group][0]
        }

        atoms_to_add = set()
        bonds_to_add = set()
        for m in atoms:
            a = self.molecule.GetAtom(oechem.OEHasMapIdx(m))
            for nbr in a.GetAtoms():
                nbr_m_idx = nbr.GetMapIdx()
                if not nbr.IsHydrogen() and nbr_m_idx not in atoms:
                    # This is a cut. If atom is N, O or S, it needs to be capped
                    if (
                        a.GetAtomicNum() in (7, 8, 16)
                        or m in map_index_to_functional_group
                    ):
                        # Add all carbons bonded to this atom
                        for nbr_2 in a.GetAtoms():
                            if nbr_2.IsCarbon():
                                atoms_to_add.add(nbr_2.GetMapIdx())
                                bonds_to_add.add((a.GetMapIdx(), nbr_2.GetMapIdx()))
        atoms.update(atoms_to_add)
        bonds.update(bonds_to_add)
        return atoms, bonds

    def _get_bond(self, bond_tuple):
        """
        Get bond in molecule by atom indices of atom A and atom B

        Parameters
        ----------
        bond_tuple : tuple
            (mapidx, mapidx)

        Returns
        -------
        bond: oechem.OEBondBase
            The bond in the molecule given by the bond tuple

        """
        from openeye import oechem

        if is_missing_atom_map(self.molecule):
            restore_atom_map(self.molecule)
        a1 = self.molecule.GetAtom(oechem.OEHasMapIdx(bond_tuple[0]))
        a2 = self.molecule.GetAtom(oechem.OEHasMapIdx(bond_tuple[-1]))
        bond = self.molecule.GetBond(a1, a2)
        if not bond:
            raise ValueError("({}) atoms are not connected".format(bond_tuple))
        return bond

    @abc.abstractmethod
    def fragment(self):
        """Fragment molecules."""
        raise NotImplementedError()

    def get_provenance(self):
        """
        Get version of fragmenter and options used

        """
        import getpass
        import socket
        import uuid

        from openeye import oechem

        import fragmenter

        fragmenter_version = fragmenter.__version__
        # Restore map to parent
        restore_atom_map(self.molecule)
        provenance = {
            "creator": fragmenter.__package__,
            "job_id": str(uuid.uuid4()),
            "hostname": socket.gethostname(),
            "username": getpass.getuser(),
            "routine": {
                "fragment_molecule": {
                    "version": fragmenter_version,
                    "options": self._options,
                    "parent_molecule": mol_to_smiles(
                        self.molecule, mapped=False, explicit_hydrogen=False
                    ),
                    "parent_name": self.molecule.GetTitle(),
                    "mapped_parent_smiles": oechem.OEMolToSmiles(self.molecule),
                }
            },
        }
        return provenance

    def _to_qcschema_mol(self, molecule, **kwargs):
        """

        Parameters
        ----------
        molecule : OEMol
        kwargs :

        Returns
        -------
        qcschema initial molecules, cmiles identifiers and provenance

        """
        from openeye import oechem

        self._options.update(kwargs)
        mol_copy = oechem.OEMol(molecule)
        oechem.OEAddExplicitHydrogens(mol_copy)
        explicit_h_smiles = mol_to_smiles(mol_copy, mapped=False)
        cmiles_identifiers = cmiles.get_molecule_ids(explicit_h_smiles)
        can_mapped_smiles = cmiles_identifiers[
            "canonical_isomeric_explicit_hydrogen_mapped_smiles"
        ]

        conformers = generate_conformers(
            to_off_molecule(mol_copy), **kwargs
        ).to_openeye()
        qcschema_mols = [
            mol_to_map_ordered_qcschema(conf, can_mapped_smiles)
            for conf in conformers.GetConfs()
        ]

        return {
            "initial_molecule": qcschema_mols,
            "identifiers": cmiles_identifiers,
            "provenance": self.get_provenance(),
        }

    def to_torsiondrive_json(self, **kwargs):
        """
        Generates torsiondrive input JSON for QCArchive

        Returns
        -------
        torsiondrive_json_dict: dict
            dictionary with the QCArchive job label as keys that maps to the torsiondrive input for each fragment

        """
        from openeye import oechem

        # capture options
        self._options.update(kwargs)
        torsiondrive_json_dict = {}
        for bond in self.fragments:
            # find torsion around bond in fragment
            molecule = self.fragments[bond]
            mapped_smiles = oechem.OEMolToSmiles(molecule)
            torsion_map_idx = torsions.find_torsion_around_bond(molecule, bond)
            torsiondrive_job_label = to_canonical_label(mapped_smiles, torsion_map_idx)
            torsiondrive_json_dict[torsiondrive_job_label] = {}

            # prepare torsiondrive input
            mol_copy = oechem.OEMol(molecule)

            # Map torsion to canonical ordered molecule
            # First canonical order the molecule
            # Note: potential bug - add explicit hydrogen before order atoms
            oechem.OEAddExplicitHydrogens(mol_copy)
            cmiles._cmiles_oe.canonical_order_atoms(mol_copy)
            dih = []
            for m_idx in torsion_map_idx:
                atom = mol_copy.GetAtom(oechem.OEHasMapIdx(m_idx + 1))
                dih.append(atom.GetIdx())

            torsiondrive_json_dict[torsiondrive_job_label] = self._to_qcschema_mol(
                mol_copy, **kwargs
            )
            torsiondrive_json_dict[torsiondrive_job_label].update(
                {"dihedral": [dih], "grid": [15], "provenance": self.get_provenance()}
            )
        return torsiondrive_json_dict


class WBOFragmenter(Fragmenter):
    """
    Fragment engine for fragmenting molecules using Wiberg Bond Order

    Parameters
    ----------
    molecule : OEMol
        Molecule to fragment.
    functional_groups : dict, optional, default None
        `{f_group: SMARTS}`. Dictionary that maps the name of a functional group to its SMARTS pattern.
        These functional groups, if they exist in the molecule, will be tagged so they are not fragmented.
        If None, will use internal list of functional group. If False, will not tag any functional groups.

    """

    def __init__(self, molecule, functional_groups=None, verbose=False):

        if functional_groups is None:
            functional_groups = get_fgroup_smarts()

        super().__init__(molecule, functional_groups)

        self.rotors_wbo = {}

        self.verbose = verbose
        self.threshold = None

    def fragment(
        self, threshold=0.03, keep_non_rotor_ring_substituents=False, **kwargs
    ):
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
        self._options[
            "keep_non_rotor_ring_substituents"
        ] = keep_non_rotor_ring_substituents
        if "threshold" not in self._options:
            self._options["threshold"] = threshold

        # Add threshold as attribute because it is used in more than one function
        setattr(self, "threshold", threshold)
        self._options.update(kwargs)
        # Calculate WBO for molecule
        self.calculate_wbo(**kwargs)
        self._get_rotor_wbo()
        # Find ring systems
        self._find_ring_systems(
            keep_non_rotor_ring_substituents=keep_non_rotor_ring_substituents
        )

        # Build fragments
        for bond in self.rotors_wbo:
            self._build_fragment(bond, **kwargs)

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
            self.molecule = assign_elf10_am1_bond_orders(self.molecule, **kwargs)
            time2 = time.time()
            if self.verbose:
                logger.info("WBO took {} seconds to calculate".format(time2 - time1))
        else:
            time1 = time.time()
            fragment = assign_elf10_am1_bond_orders(fragment, **kwargs)
            time2 = time.time()
            if self.verbose:
                logger.info("WBO took {} seconds to calculate".format(time2 - time1))
        return fragment

    def _get_rotor_wbo(self):
        """
        Cache WBO for each rotatable bond
        """

        if (
            "nrotor" not in self.molecule.GetData()
            and "cputime" not in self.molecule.GetData()
        ):
            logger.info("WBO was not calculated for this molecule. Calculating WBO...")
            self.calculate_wbo()
        rotatable_bonds = self._find_rotatable_bonds()
        for bond in rotatable_bonds:
            b = self._get_bond(bond)
            self.rotors_wbo[bond] = b.GetData("WibergBondOrder")

    def _compare_wbo(self, fragment, bond_tuple, **kwargs):
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
        from openeye import oechem

        # Create new oemol because sometimes the molecule created from atom bond set is wonky and then the WBOs are not reproducible
        smiles = oechem.OEMolToSmiles(fragment)
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, smiles)
        try:
            charged_fragment = self.calculate_wbo(fragment=mol, **kwargs)
        except RuntimeError:
            logger.warn(
                "Cannot calculate WBO for fragment {}. Continue growing fragment".format(
                    oechem.OEMolToSmiles(mol)
                )
            )
            # Most of the time it fails because it is either missing parameters or a functional group that should not
            # be fragmented was fragmented
            # ToDo:  hanlde different kinds of failures instead of just continuing to grow until the fialure goes away
            # Some fail because there are functional groups that should not be fragmented.
            return 1

        # Get new WBO
        restore_atom_map(charged_fragment)
        a1 = charged_fragment.GetAtom(oechem.OEHasMapIdx(bond_tuple[0]))
        a2 = charged_fragment.GetAtom(oechem.OEHasMapIdx(bond_tuple[-1]))
        bond = charged_fragment.GetBond(a1, a2)
        if bond is None:
            raise RuntimeError(
                "{} with _idx {} and {} with map_idx {} are not bonded".format(
                    a1, bond_tuple[0], a2, bond_tuple[1]
                )
            )
        if "WibergBondOrder" not in bond.GetData():
            logger.warn(
                "Cannot calculate WBO for fragment {}. Continue growing fragment".format(
                    oechem.OEMolToSmiles(charged_fragment)
                )
            )
            return 1

        frag_wbo = bond.GetData("WibergBondOrder")
        mol_wbo = self.rotors_wbo[bond_tuple]

        return abs(mol_wbo - frag_wbo), charged_fragment

    def _build_fragment(self, bond_tuple, heuristic="path_length", cap=True, **kwargs):
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
        if "heuristic" not in self._options:
            self._options["heuristic"] = heuristic
        atoms, bonds = self._get_torsion_quartet(bond_tuple)
        atom_map_idx, bond_tuples = self._get_ring_and_fgroups(atoms, bonds)

        # Cap open valence
        if cap:
            atom_map_idx, bond_tuples = self._cap_open_valence(
                atom_map_idx, bond_tuples, bond_tuple
            )

        fragment_mol = self._atom_bond_set_to_mol(atom_map_idx, bond_tuples)
        # #return fragment_mol
        diff, fragment_mol = self._compare_wbo(fragment_mol, bond_tuple, **kwargs)
        if diff <= self.threshold:
            self.fragments[bond_tuple] = fragment_mol
        if diff > self.threshold:
            self._add_next_substituent(
                atom_map_idx,
                bond_tuples,
                target_bond=bond_tuple,
                heuristic=heuristic,
                **kwargs,
            )

    def _add_next_substituent(
        self, atoms, bonds, target_bond, heuristic="path_length", cap=True, **kwargs
    ):
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
        from openeye import oechem

        map_index_to_ring_system = {
            map_index: ring_system
            for ring_system in self.ring_systems
            for map_index in self.ring_systems[ring_system][0]
        }
        map_index_to_functional_group = {
            map_index: functional_group
            for functional_group in self.functional_groups
            for map_index in self.functional_groups[functional_group][0]
        }

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
                if not nbr.IsHydrogen() and nbr_m_idx not in atoms:
                    b = self.molecule.GetBond(a, nbr)
                    atoms_to_add.append((nbr_m_idx, (m_idx, nbr_m_idx)))
                    if heuristic == "wbo":
                        sort_by.append(b.GetData("WibergBondOrder"))
                        reverse = True
                    elif heuristic == "path_length":
                        sort_by_1.append(oechem.OEGetPathLength(bond_atom_1, nbr))
                        sort_by_2.append(oechem.OEGetPathLength(bond_atom_2, nbr))
                        reverse = False
                    else:
                        raise ValueError(
                            "Only wbo and path_lenght are supported heuristics"
                        )

        # A work around for a known bug where if stereochemistry changes or gets removed, the WBOs can change more than
        # the threshold (this will sometimes happen if a very small threshold is chosen) and even the parent will have
        # a wBO difference greater than the threshold. In this case, return the molecule
        if heuristic == "wbo" and len(sort_by) == 0:
            fragment_mol = self._atom_bond_set_to_mol(atoms, bonds)
            return fragment_mol
        if heuristic == "path_length" and len(sort_by_1) == 0 and len(sort_by_2) == 0:
            fragment_mol = self._atom_bond_set_to_mol(atoms, bonds)
            return fragment_mol

        if heuristic == "path_length":
            min_1 = min(sort_by_1)
            min_2 = min(sort_by_2)
            if min_1 < min_2:
                sort_by = sort_by_1
            elif min_2 < min_1:
                sort_by = sort_by_2
            elif min_1 == min_2:
                indices = []
                for length_index in [sort_by_1, sort_by_2]:
                    indices.extend(
                        [i for i, x in enumerate(length_index) if x == min_1]
                    )
                atoms_to_add = [atoms_to_add[i] for i in indices]
                for a, b in atoms_to_add:
                    bond = self._get_bond(b)
                    sort_by.append(bond.GetData("WibergBondOrder"))
                    reverse = True

        sorted_atoms = [
            a for _, a in sorted(zip(sort_by, atoms_to_add), reverse=reverse)
        ]
        for atom, bond in sorted_atoms:
            a = self.molecule.GetAtom(oechem.OEHasMapIdx(atom))
            if a.GetMapIdx() in map_index_to_ring_system:
                ring_idx = map_index_to_ring_system[a.GetMapIdx()]
                atoms.update(self.ring_systems[ring_idx][0])
                bonds.update(self.ring_systems[ring_idx][1])
            if a.GetMapIdx() in map_index_to_functional_group:
                fgroup = map_index_to_functional_group[a.GetMapIdx()]
                atoms.update(self.functional_groups[fgroup][0])
                bonds.update(self.functional_groups[fgroup][1])
            atoms.add(atom)
            bonds.add(bond)
            # Check new WBO
            # cap open valence
            # if cap:
            #     atoms, bonds = self._cap_open_valence(atoms, bonds, target_bond)
            fragment_mol = self._atom_bond_set_to_mol(atoms, bonds)

            diff, fragment_mol = self._compare_wbo(fragment_mol, target_bond, **kwargs)

            if diff < self.threshold:
                self.fragments[target_bond] = fragment_mol
                return fragment_mol
            else:
                return self._add_next_substituent(
                    atoms=atoms,
                    bonds=bonds,
                    target_bond=target_bond,
                    heuristic=heuristic,
                    **kwargs,
                )


class PfizerFragmenter(Fragmenter):
    """
    Fragment engine for fragmenting molecules using Pfizer's protocol (doi: 10.1021/acs.jcim.9b00373)

    Parameters
    ----------
    molecule : OEMol
        Molecule to fragment.
    """

    def __init__(self, molecule):
        super().__init__(molecule, get_fgroup_smarts())
        self._find_ring_systems(keep_non_rotor_ring_substituents=False)

    def fragment(self):
        """
        Fragment molecules according to Pfizer protocol
        """
        rotatable_bonds = self._find_rotatable_bonds()
        for bond in rotatable_bonds:

            atoms, bonds = self._get_torsion_quartet(bond)
            atoms, bonds = self._get_ring_and_fgroups(atoms, bonds)
            atoms, bonds = self._cap_open_valence(atoms, bonds, bond)

            self.fragments[bond] = self._atom_bond_set_to_mol(
                atoms, bonds, adjust_hcount=True
            )
