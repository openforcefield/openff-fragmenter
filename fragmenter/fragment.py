import abc
import logging
import time
from collections import defaultdict
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import networkx
from openff.toolkit.topology import Atom, Molecule
from typing_extensions import Literal

from fragmenter.chemi import (
    assign_elf10_am1_bond_orders,
    extract_fragment,
    find_ring_systems,
    find_stereocenters,
)
from fragmenter.states import _enumerate_stereoisomers
from fragmenter.utils import get_atom_index, get_fgroup_smarts, get_map_index

logger = logging.getLogger(__name__)

BondTuple = Tuple[int, int]
AtomAndBondSet = Tuple[Set[int], Set[BondTuple]]

Heuristic = Literal["path_length", "wbo"]


class Fragmenter(abc.ABC):
    """Base fragmenter class."""

    def __init__(self, molecule: Molecule, functional_groups):

        # Canonically order the molecule to try and make the fragmentation more
        # deterministic.
        self.molecule: Molecule = molecule.canonical_order_atoms()
        self.molecule.properties["atom_map"] = {
            i: i + 1 for i in range(self.molecule.n_atoms)
        }

        # Keep track of stereo to make sure it does not flip
        self.stereo = {}
        self._find_stereo()

        # For provenance
        self._options = {}

        # Find the functional groups.
        self.functional_groups = {}

        self._options["functional_groups"] = functional_groups
        self._tag_functional_groups(functional_groups)

        # Track any ring systems.
        self.ring_systems = {}

        # Fragments from fragmentation scheme for each rotatable bond
        self.fragments: Dict[BondTuple, Molecule] = {}

    def _find_stereo(self):
        """Find chiral atoms and bonds, store the chirality. This is needed to check if
        fragments flipped chirality. Currently this can happen and it is a bug
        """

        atom_stereo = {
            get_map_index(self.molecule, atom.molecule_atom_index): atom.stereochemistry
            for atom in self.molecule.atoms
            if atom.stereochemistry is not None
        }

        bond_stereo = {
            (
                get_map_index(self.molecule, bond.atom1_index),
                get_map_index(self.molecule, bond.atom2_index),
            ): bond.stereochemistry
            for bond in self.molecule.bonds
            if bond.stereochemistry is not None
        }

        self.stereo = {**atom_stereo, **bond_stereo}

    def _check_stereo(self, fragment: Molecule) -> bool:
        """Check if stereo in fragment is different than stereo in parent.

        Parameters
        ----------
        fragment
            The fragment whose stereo should be compared to the parent.

        Returns
        -------
            Whether the fragment has the same stereo as the parent.
        """

        atom_stereocenters, bond_stereocenters = find_stereocenters(fragment)

        # Check for new / flipped chiral centers.
        for atom_index in atom_stereocenters:

            map_index = get_map_index(fragment, atom_index)

            if map_index not in self.stereo:

                logger.warning(f"A new stereocenter formed at atom {map_index}")
                return False

            fragment_stereo = fragment.atoms[atom_index].stereochemistry
            parent_stereo = self.stereo[map_index]

            if fragment_stereo != parent_stereo:

                logger.warning(
                    f"Stereochemistry for atom {map_index} flipped from "
                    f"{parent_stereo} to {fragment_stereo}"
                )

                return False

        for index_tuple in bond_stereocenters:

            map_tuple = tuple(get_map_index(fragment, i) for i in index_tuple)

            map_tuple = (
                map_tuple if map_tuple in self.stereo else tuple(reversed(map_tuple))
            )

            if map_tuple not in self.stereo:

                logger.warning(f"A new chiral bond formed at bond {map_tuple}")
                return False

            fragment_stereo = fragment.get_bond_between(*index_tuple).stereochemistry

            parent_stereo = self.stereo[map_tuple]

            if fragment_stereo != parent_stereo:

                logger.warning(
                    f"Stereochemistry for bond {map_tuple} flipped from "
                    f"{parent_stereo} to {fragment_stereo}"
                )

                return False

        return True

    def _fix_stereo(self, fragment: Molecule) -> Molecule:
        """Flip all stereocenters and find the stereoisomer that matches the parent

        Parameters
        ----------
        fragment
            The fragment whose stereochemistry should be flipped to match the parent.

        Returns
        -------
            A fragment with the same stereochemistry as parent molecule
        """

        for stereoisomer in _enumerate_stereoisomers(fragment, 200, True):

            if not self._check_stereo(stereoisomer):
                continue

            return stereoisomer

        raise RuntimeError(
            f"The stereochemistry of {fragment.to_smiles()} could not be fixed."
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

        for functional_group, smarts in functional_groups.items():

            unique_matches = {
                tuple(sorted(match))
                for match in self.molecule.chemical_environment_matches(smarts)
            }

            for i, match in enumerate(unique_matches):

                atoms = set(get_map_index(self.molecule, index) for index in match)
                bonds = set(
                    (
                        get_map_index(self.molecule, bond.atom1_index),
                        get_map_index(self.molecule, bond.atom2_index),
                    )
                    for bond in self.molecule.bonds
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

        matches = self.molecule.chemical_environment_matches(
            "[!$(*#*)&!D1:1]-,=;!@[!$(*#*)&!D1:2]"
        )
        unique_matches = {tuple(sorted(match)) for match in matches}

        # Drop bonds without a heavy degree of at least 2 on each end to avoid finding
        # terminal bonds
        def heavy_degree(atom_index: int) -> int:
            atom = self.molecule.atoms[atom_index]
            return sum(1 for atom in atom.bonded_atoms if atom.atomic_number != 1)

        unique_matches = {
            match for match in unique_matches if all(heavy_degree(i) > 1 for i in match)
        }

        return [
            (
                get_map_index(self.molecule, match[0]),
                get_map_index(self.molecule, match[1]),
            )
            for match in unique_matches
        ]

    def _atom_bond_set_to_mol(self, atoms: Set[int], bonds: Set[BondTuple]) -> Molecule:
        """Extracts a subset of a molecule based on a set of atom and bond indices.

        Parameters
        ----------
        atoms
            set of map indices
        bonds
            set of bond tuples (m1, m2)

        Returns
        -------
            The subset molecule.
        """

        fragment = extract_fragment(self.molecule, atoms, bonds)

        if not self._check_stereo(fragment):
            fragment = self._fix_stereo(fragment)

        return fragment

    def _get_torsion_quartet(self, bond: BondTuple) -> AtomAndBondSet:
        """Get all atoms bonded to the torsion quartet around rotatable bond

        Parameters
        ----------
        bond: tuple of ints
            map indices of atoms in bond

        Returns
        -------
            The map indices of atoms in quartet and the bonds in quartet)
        """

        atom_map_indices = {*bond}
        bond_map_indices = {bond}

        atoms = [
            self.molecule.atoms[i]
            for i, j in self.molecule.properties["atom_map"].items()
            if j in bond
        ]

        for atom in atoms:

            map_index = get_map_index(self.molecule, atom.molecule_atom_index)

            for neighbor in atom.bonded_atoms:

                neighbour_map_index = get_map_index(
                    self.molecule, neighbor.molecule_atom_index
                )

                atom_map_indices.add(neighbour_map_index)
                bond_map_indices.add((map_index, neighbour_map_index))

                for next_neighbour in neighbor.bonded_atoms:

                    next_neighbour_map_index = get_map_index(
                        self.molecule, next_neighbour.molecule_atom_index
                    )

                    atom_map_indices.add(next_neighbour_map_index)
                    bond_map_indices.add(
                        (neighbour_map_index, next_neighbour_map_index)
                    )

        return atom_map_indices, bond_map_indices

    def _find_ring_systems(self, keep_non_rotor_ring_substituents: bool = False):
        """This function finds all ring systems and stores them in the `ring_systems`
        field.

        Parameters
        ----------
        keep_non_rotor_ring_substituents
            If True, keep all non rotatable ring substituents. According to the
            benchmark, it is not necessary.
        """

        atom_to_ring_indices = find_ring_systems(self.molecule)

        # Find the map indices of the atoms involved in each ring system.
        ring_system_atoms = {
            ring_index: {
                get_map_index(self.molecule, i)
                for i in atom_to_ring_indices
                if atom_to_ring_indices[i] == ring_index
            }
            for ring_index in {*atom_to_ring_indices.values()}
        }

        # Find the map indices of the bonds involved in each ring system.
        ring_system_bonds = defaultdict(set)

        for bond in self.molecule.bonds:

            ring_index_1 = atom_to_ring_indices.get(bond.atom1_index, -1)
            ring_index_2 = atom_to_ring_indices.get(bond.atom2_index, -2)

            if ring_index_1 != ring_index_2:
                continue

            ring_system_bonds[ring_index_1].add(
                (
                    get_map_index(self.molecule, bond.atom1_index),
                    get_map_index(self.molecule, bond.atom2_index),
                )
            )

        # Scan the neighbours of the ring system atoms for any functional groups
        # / non-rotor substituents which should be included in the ring systems.
        for ring_index in ring_system_atoms:

            # If any atoms are part of a functional group, include the other atoms in the
            # group in the ring system lists
            functional_groups = {
                functional_group
                for map_index in ring_system_atoms[ring_index]
                for functional_group in self.functional_groups
                if map_index in self.functional_groups[functional_group][0]
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
                ring_system_atoms[ring_index]
            )

            ring_system_atoms[ring_index].update(non_rotor_atoms)
            ring_system_bonds[ring_index].update(non_rotor_bonds)

        for ring_index in ring_system_atoms:

            self.ring_systems[ring_index] = (
                ring_system_atoms[ring_index],
                ring_system_bonds[ring_index],
            )

    def _find_non_rotor_ring_substituents(
        self, ring_system_atoms: Set[int]
    ) -> AtomAndBondSet:
        """Find the non-rotor substituents attached to a particular ring system."""

        rotatable_bonds = self.molecule.find_rotatable_bonds()

        def heavy_degree(atom: Atom) -> int:
            return sum(1 for atom in atom.bonded_atoms if atom.atomic_number != 1)

        rotor_bonds = [
            bond
            for bond in rotatable_bonds
            if heavy_degree(bond.atom1) >= 2 and heavy_degree(bond.atom2) >= 2
        ]

        non_rotor_atoms = set()
        non_rotor_bonds = set()

        for bond in self.molecule.bonds:

            # Check if the bond is a rotor.
            if bond in rotor_bonds:
                continue

            if bond.atom1.atomic_number == 1 or bond.atom2.atomic_number == 1:
                continue

            map_index_1 = get_map_index(self.molecule, bond.atom1_index)
            map_index_2 = get_map_index(self.molecule, bond.atom2_index)

            in_system_1 = map_index_1 in ring_system_atoms
            in_system_2 = map_index_2 in ring_system_atoms

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
            The updated set of atom and bond map indices to retain.
        """

        # Find the sets of atoms which are located ortho to one of the bonds being
        # fragmented.
        ortho_atoms, ortho_bonds = self._find_ortho_substituents(bonds)

        atoms.update(ortho_atoms)
        bonds.update(ortho_bonds)

        # Include the rings systems and functional groups connected to the current
        # atom sets.
        new_atoms = set()
        new_bonds = set()

        functional_groups = {
            group
            for group in self.functional_groups
            if any(atom in self.functional_groups[group][0] for atom in atoms)
        }

        for functional_group in functional_groups:
            new_atoms.update(self.functional_groups[functional_group][0])
            new_bonds.update(self.functional_groups[functional_group][1])

        ring_systems = {
            ring_index
            for ring_index in self.ring_systems
            if any(atom in self.ring_systems[ring_index][0] for atom in atoms)
        }

        for ring_system in ring_systems:
            new_atoms.update(self.ring_systems[ring_system][0])
            new_bonds.update(self.ring_systems[ring_system][1])

        atoms.update(new_atoms)
        bonds.update(new_bonds)

        # Ensure the matched bonds doesn't include duplicates.
        bonds = {tuple(sorted(bond)) for bond in bonds}

        return atoms, bonds

    def _find_ortho_substituents(self, bonds: Set[BondTuple]) -> AtomAndBondSet:
        """Find ring substituents that are ortho to one of the rotatable bonds specified
        in a list of bonds.

        Parameters
        ----------
        bonds
            The map indices of the rotatable bonds.

        Returns
        -------
            The set of map indices of atoms in ortho group and of bond tuples in ortho
            group.
        """

        matched_atoms = set()
        matched_bonds = set()

        for match in self.molecule.chemical_environment_matches(
            "[!#1:1]~&!@[*:2]@[*:3]~&!@[!#1*:4]"
        ):

            map_tuple = tuple(get_map_index(self.molecule, i) for i in match)

            if map_tuple[:2] not in bonds and map_tuple[:2][::-1] not in bonds:
                continue

            matched_atoms.update(map_tuple[::3])
            matched_bonds.update((map_tuple[i], map_tuple[i + 1]) for i in [0, 2])

        # Ensure the matched bonds doesn't include duplicates.
        matched_bonds = {tuple(sorted(bond)) for bond in matched_bonds}

        return matched_atoms, matched_bonds

    def _cap_open_valence(
        self, atoms: Set[int], bonds: Set[BondTuple]
    ) -> AtomAndBondSet:
        """Cap with methyl for fragments that ends with N, O or S. Otherwise cap with H

        Parameters
        ----------
        atoms
            The map indices of the atoms in the fragment being constructed.
        bonds
            The map indices of the bonds in the fragment being constructed.
        """

        map_index_to_functional_group = {
            map_index: functional_group
            for functional_group in self.functional_groups
            for map_index in self.functional_groups[functional_group][0]
        }

        atoms_to_add = set()
        bonds_to_add = set()

        for map_index in atoms:

            atom_index = get_atom_index(self.molecule, map_index)
            atom = self.molecule.atoms[atom_index]

            if (
                atom.atomic_number not in (7, 8, 16)
                and map_index not in map_index_to_functional_group
            ):
                continue

            # If atom is N, O or S, it needs to be capped
            should_cap = False

            for neighbour in atom.bonded_atoms:

                neighbour_map_index = get_map_index(
                    self.molecule, neighbour.molecule_atom_index
                )

                if neighbour.atomic_number == 1 or neighbour_map_index in atoms:
                    continue

                should_cap = True
                break

            if not should_cap:
                continue

            for neighbour in atom.bonded_atoms:

                if neighbour.atomic_number != 6:
                    continue

                neighbour_map_index = get_map_index(
                    self.molecule, neighbour.molecule_atom_index
                )

                atoms_to_add.add(neighbour_map_index)
                bonds_to_add.add((map_index, neighbour_map_index))

        atoms.update(atoms_to_add)
        bonds.update(bonds_to_add)

        return atoms, bonds

    @abc.abstractmethod
    def fragment(self):
        """Fragment molecules."""
        raise NotImplementedError()

    def get_provenance(self) -> Dict[str, Any]:
        """
        Get version of fragmenter and options used

        """
        import getpass
        import socket
        import uuid

        import fragmenter

        fragmenter_version = fragmenter.__version__

        provenance = {
            "creator": fragmenter.__package__,
            "job_id": str(uuid.uuid4()),
            "hostname": socket.gethostname(),
            "username": getpass.getuser(),
            "routine": {
                "fragment_molecule": {
                    "version": fragmenter_version,
                    "options": self._options,
                    "parent_molecule": self.molecule.to_smiles(
                        mapped=False, explicit_hydrogens=False
                    ),
                    "parent_name": self.molecule.name,
                    "mapped_parent_smiles": self.molecule.to_smiles(mapped=True),
                }
            },
        }
        return provenance


class WBOFragmenter(Fragmenter):
    """
    Fragment engine for fragmenting molecules using Wiberg Bond Order

    Parameters
    ----------
    molecule : OEMol
        Molecule to fragment.
    functional_groups : dict, optional, default None
        `{f_group: SMARTS}`. Dictionary that maps the name of a functional group to its
        SMARTS pattern. These functional groups, if they exist in the molecule, will be
        tagged so they are not fragmented. If None, will use internal list of functional
        group. If False, will not tag any functional groups.

    """

    def __init__(self, molecule, functional_groups=None, verbose=False):

        if functional_groups is None:
            functional_groups = get_fgroup_smarts()

        super().__init__(molecule, functional_groups)

        self.rotors_wbo = {}

        self.verbose = verbose
        self.threshold = None

    def fragment(
        self,
        threshold: float = 0.03,
        keep_non_rotor_ring_substituents: bool = False,
        **kwargs,
    ):
        """
        Fragment molecules using the Wiberg Bond Order as a surrogate

        Parameters
        ----------
        threshold
            The threshold for the central bond WBO. If the fragment WBO is below this
            threshold, fragmenter will grow out the fragment one bond at a time via the
            path specified by the heuristic option
        keep_non_rotor_ring_substituents: bool
            If True, will always keep all non rotor substituents on ring. If False, will
            only add them if they are ortho to rotatable bond or if it's needed for WBO
            to be within the threshold
        **heuristic : str, optional, default 'path_length'
            The path fragmenter should take when fragment needs to be grown out. The
            other option is 'wbo'
        """
        # Capture options used
        self._options[
            "keep_non_rotor_ring_substituents"
        ] = keep_non_rotor_ring_substituents
        if "threshold" not in self._options:
            self._options["threshold"] = threshold

        # Add threshold as attribute because it is used in more than one function
        self.threshold = threshold
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

    def calculate_wbo(
        self, fragment: Optional[Molecule] = None, **kwargs
    ) -> Optional[Molecule]:
        """Calculate the WBOs on either a parent or fragment molecule.

        Parameters
        ----------
        fragment
            The fragment to recalculate the WBO for. When fragment is ``None``,
            fragmenter assumes it's the full molecule and saves the calculated values
            in self.molecule

        Returns
        -------
            The fragment with WBOs.
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
        """Cache the WBO of each rotatable bond"""

        if any(bond.fractional_bond_order is None for bond in self.molecule.bonds):

            logger.info("WBO was not calculated for this molecule. Calculating WBO...")
            self.calculate_wbo()

        rotatable_bonds = self._find_rotatable_bonds()

        for bond_indices in rotatable_bonds:

            bond = self.molecule.get_bond_between(
                get_atom_index(self.molecule, bond_indices[0]),
                get_atom_index(self.molecule, bond_indices[1]),
            )

            self.rotors_wbo[bond_indices] = bond.fractional_bond_order

    def _compare_wbo(
        self, fragment: Molecule, bond_tuple: BondTuple, **kwargs
    ) -> float:
        """Compare Wiberg Bond order of rotatable bond in a fragment to the parent.

        Parameters
        ----------
        fragment
            The fragment containing the rotatable bond.
        bond_tuple
            The map indices of the rotatable bond.

        Returns
        -------
            The absolute difference between the fragment and parent WBOs.
        """

        # Create new fragment object because sometimes the molecule created from atom
        # bond set is wonky and then the WBOs are not reproducible
        fragment = Molecule.from_smiles(
            fragment.to_smiles(mapped=True), allow_undefined_stereo=True
        )

        fragment_map = fragment.properties.pop("atom_map", None)

        try:
            fragment = self.calculate_wbo(fragment=fragment, **kwargs)

        except RuntimeError:

            # Most of the time it fails because it is either missing parameters or a
            # functional group that should not be fragmented was fragmented
            logger.warning(
                f"Cannot calculate WBO for fragment {fragment.to_smiles()}. Continue "
                f"growing fragment"
            )

            # TODO: handle different kinds of failures instead of just continuing to
            #      grow until the failure goes away. Some fail because there are
            #      functional groups that should not be fragmented.
            return 1.0

        if fragment_map is not None:
            fragment.properties["atom_map"] = fragment_map

        bond = fragment.get_bond_between(
            get_atom_index(fragment, bond_tuple[0]),
            get_atom_index(fragment, bond_tuple[1]),
        )

        fragment_wbo = bond.fractional_bond_order
        parent_wbo = self.rotors_wbo[bond_tuple]

        return abs(parent_wbo - fragment_wbo)

    def _build_fragment(
        self,
        bond_tuple: BondTuple,
        heuristic: Heuristic = "path_length",
        cap: bool = True,
        **kwargs,
    ):
        """Build a fragment around a specified bond.

        Parameters
        ----------
        bond_tuple
            The map indices specifying which bond to build the fragment around.
        heuristic
            The heuristic to use when building the fragment.
        cap
            Whether to cap open valences.
        """

        # Capture options
        if "heuristic" not in self._options:
            self._options["heuristic"] = heuristic

        atoms, bonds = self._get_torsion_quartet(bond_tuple)
        atoms, bonds = self._get_ring_and_fgroups(atoms, bonds)

        # Cap open valence
        if cap:
            atoms, bonds = self._cap_open_valence(atoms, bonds)

        fragment = self._atom_bond_set_to_mol(atoms, bonds)

        wbo_difference = self._compare_wbo(fragment, bond_tuple, **kwargs)

        while fragment is not None and wbo_difference > self.threshold:

            fragment = self._add_next_substituent(
                self.molecule,
                atoms,
                bonds,
                target_bond=bond_tuple,
                heuristic=heuristic,
            )

            if fragment is None:
                break

            wbo_difference = self._compare_wbo(fragment, bond_tuple, **kwargs)

        # A work around for a known bug where if stereochemistry changes or gets removed,
        # the WBOs can change more than the threshold (this will sometimes happen if a
        # very small threshold is chosen) and even the parent will have a WBO difference
        # greater than the threshold. In this case, return the molecule
        if fragment is None:
            fragment = self._atom_bond_set_to_mol(atoms, bonds)

        self.fragments[bond_tuple] = fragment

    @classmethod
    def _select_neighbour_by_path_length(
        cls, molecule: Molecule, atoms: Set[int], target_bond: BondTuple
    ) -> Optional[Tuple[int, BondTuple]]:

        atom_indices = {get_atom_index(molecule, atom) for atom in atoms}

        atoms_to_add = [
            (atom_index, neighbour.molecule_atom_index)
            for atom_index in atom_indices
            for neighbour in molecule.atoms[atom_index].bonded_atoms
            if neighbour.atomic_number != 1
            and neighbour.molecule_atom_index not in atom_indices
        ]
        map_atoms_to_add = [
            (
                get_map_index(molecule, j),
                (get_map_index(molecule, i), get_map_index(molecule, j)),
            )
            for i, j in atoms_to_add
        ]

        # Compute the distance from each neighbouring atom to each of the atoms in the
        # target bond.
        nx_molecule = molecule.to_networkx()

        target_indices = [get_atom_index(molecule, atom) for atom in target_bond]

        path_lengths_1, path_lengths_2 = zip(
            *(
                (
                    networkx.shortest_path_length(
                        nx_molecule, target_index, neighbour_index
                    )
                    for target_index in target_indices
                )
                for atom_index, neighbour_index in atoms_to_add
            )
        )

        if len(path_lengths_1) == 0 and len(path_lengths_2) == 0:
            return None

        reverse = False

        min_path_length_1 = min(path_lengths_1)
        min_path_length_2 = min(path_lengths_2)

        if min_path_length_1 < min_path_length_2:
            sort_by = path_lengths_1
        elif min_path_length_2 < min_path_length_1:
            sort_by = path_lengths_2

        else:

            # If there are multiple neighbouring atoms the same path length away
            # from the target bond fall back to sorting by the WBO.
            map_atoms_to_add = [
                map_tuple
                for map_tuple, *path_length_tuple in zip(
                    map_atoms_to_add, path_lengths_1, path_lengths_2
                )
                if min_path_length_1 in path_length_tuple
            ]

            sort_by = [
                molecule.get_bond_between(
                    get_atom_index(molecule, neighbour_bond[0]),
                    get_atom_index(molecule, neighbour_bond[1]),
                ).fractional_bond_order
                for _, neighbour_bond in map_atoms_to_add
            ]

            reverse = True

        sorted_atoms = [
            a for _, a in sorted(zip(sort_by, map_atoms_to_add), reverse=reverse)
        ]

        return None if len(sorted_atoms) == 0 else sorted_atoms[0]

    @classmethod
    def _select_neighbour_by_wbo(
        cls, molecule: Molecule, atoms: Set[int]
    ) -> Optional[Tuple[int, BondTuple]]:
        """A function which return those atoms which neighbour those in the ``atoms``
        list sorted by the WBO of the bond between the input atom and the neighbouring
        atom from largest to smallest.

        Parameters
        ----------
        molecule
            The original molecule being fragmented.
        atoms
            The map indices of the atoms currently in the fragment.

        Returns
        -------
            The indices of the atoms to be added to the fragment sorted into the
            order that they should be added in.
        """

        map_bond_orders = {
            (
                get_map_index(molecule, bond.atom1_index),
                get_map_index(molecule, bond.atom2_index),
            ): bond.fractional_bond_order
            for bond in molecule.bonds
            if bond.atom1.atomic_number != 1 and bond.atom2.atomic_number != 1
        }

        neighbour_bond_orders = {
            (bond_order, (map_tuple[1 - i], map_tuple))
            for i in range(2)
            for map_tuple, bond_order in map_bond_orders.items()
            if map_tuple[i] in atoms and map_tuple[1 - i] not in atoms
        }

        sorted_atoms = [
            atom_to_add
            for _, atom_to_add in sorted(
                neighbour_bond_orders, key=lambda x: x[0], reverse=True
            )
        ]

        return None if len(sorted_atoms) == 0 else sorted_atoms[0]

    def _add_next_substituent(
        self,
        molecule: Molecule,
        atoms: Set[int],
        bonds: Set[BondTuple],
        target_bond: BondTuple,
        heuristic: Heuristic = "path_length",
    ) -> Molecule:
        """Expand the fragment to include the next set of substituents / ring systems.

        Parameters
        ----------
        molecule
            The original molecule being fragmented.
        atoms
            The map indices of the atoms currently in the fragment.
        bonds
            The map indices of the bonds currently in the fragment.
        target_bond
            The bond the fragment is being built around.
        heuristic
            How to add substituents. The choices are `'path_length'` or `'wbo'`

        Returns
        -------
            The expanded fragment.
        """

        # Select the next atom neighbour (and the groups / rings that it is part of)
        # that should be added to the fragment.
        if heuristic == "wbo":
            neighbour_atom_and_bond = self._select_neighbour_by_wbo(molecule, atoms)
        elif heuristic == "path_length":
            neighbour_atom_and_bond = self._select_neighbour_by_path_length(
                molecule, atoms, target_bond
            )
        else:
            raise NotImplementedError(
                "Only `'wbo'` and `'path_length'` are supported heuristics."
            )

        if neighbour_atom_and_bond is None:
            return None

        neighbour_atom, neighbour_bond = neighbour_atom_and_bond

        # If the neighbour to include is part of a functional group / ring system
        # the entire group should be included in the fragment.
        for group, group_atoms in self.functional_groups.items():

            if neighbour_atom not in group_atoms[0]:
                continue

            atoms.update(self.functional_groups[group][0])
            bonds.update(self.functional_groups[group][1])

        for ring_index, ring_atoms in self.ring_systems.items():

            if neighbour_atom not in ring_atoms[0]:
                continue

            atoms.update(self.ring_systems[ring_index][0])
            bonds.update(self.ring_systems[ring_index][1])

        atoms.add(neighbour_atom)
        bonds.add(neighbour_bond)

        # Check new WBO
        return self._atom_bond_set_to_mol(atoms, bonds)


class PfizerFragmenter(Fragmenter):
    """Fragment engine for fragmenting molecules using Pfizer's protocol
    (doi: 10.1021/acs.jcim.9b00373)
    """

    def __init__(self, molecule: Molecule):
        """

        Parameters
        ----------
        molecule
            The molecule to fragment.
        """
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
            atoms, bonds = self._cap_open_valence(atoms, bonds)

            self.fragments[bond] = self._atom_bond_set_to_mol(atoms, bonds)
