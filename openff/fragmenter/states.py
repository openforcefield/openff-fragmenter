import logging

from openff.toolkit.topology import Molecule

logger = logging.getLogger(__name__)


def _enumerate_stereoisomers(
    molecule: Molecule,
    max_states: int = 200,
    force_flip: bool = True,
) -> list[Molecule]:
    """Enumerate stereoisomers

    Parameters
    ----------
    molecule
        The molecule to enumerate the stereoisomers of.
    max_states
        The max number of states to enumerate
    force_flip
        If True, will flip all stereocenters. If False, will only flip centers that are
        undefined

    Returns
    -------
        The stereoisomers including the input molecule if it had full stereochemistry
        defined.
    """

    if max_states > 200:
        raise NotImplementedError(
            "The max states must currently be less than 200 due to a hard coded maximum "
            "value in the OpenFF toolkit"
        )

    # Check if the input molecule has any undefined stereochemistry. If it does not
    # then it should be included in the list of returned stereoisomers.
    stereoisomers = molecule.enumerate_stereoisomers(
        undefined_only=not force_flip,
        max_isomers=max_states,
        rationalise=False,
    )

    if len(stereoisomers) == 0:
        # Handle the case where the input molecule does not have any stereoisomers.
        return [molecule]

    # Attach the atom map to any stereoisomers.
    if "atom_map" in molecule.properties:
        for stereoisomer in stereoisomers:
            stereoisomer.properties["atom_map"] = molecule.properties["atom_map"]

    # Check to see if the original molecule had undefined stereochemistry.
    undefined_atom_stereochemistry = any(
        atom_old.stereochemistry is None and atom_new.stereochemistry is not None
        for atom_old, atom_new in zip(molecule.atoms, stereoisomers[0].atoms)
    )
    undefined_bond_stereochemistry = any(
        bond_old.stereochemistry is None and bond_new.stereochemistry is not None
        for bond_old, bond_new in zip(molecule.bonds, stereoisomers[0].bonds)
    )

    # If the input molecule had all of the input stereochemistry defined make sure to
    # return it in addition to the newly found stereoisomers.
    if not undefined_atom_stereochemistry and not undefined_bond_stereochemistry:
        stereoisomers = [molecule] + stereoisomers[: max_states - 1]

    return stereoisomers
