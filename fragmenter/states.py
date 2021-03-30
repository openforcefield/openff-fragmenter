import logging

from openff.toolkit.topology import Molecule

from fragmenter.utils import copy_molecule

logger = logging.getLogger(__name__)


def _enumerate_stereoisomers(
    molecule,
    max_states=200,
    force_flip=True,
    enum_nitrogen=True,
    warts=False,
    verbose=True,
):
    """Enumerate stereoisomers

    Parameters
    ----------
    molecule : OEMol
    max_states : int, optional, default 200
        max number of states to enumerate
    force_flip : bool, optional, default True
        If True, will flip all steocenters. If False, will only flip centers that are undefined
    enum_nitrogen : bool, optional, default True
        Invert non-planar nitrogen
    warts : bool, optional, default True
        If True, add int to molecule name
    verbose : bool, optional, default True

    Returns
    -------
    stereoisomers: list of oemols

    """

    # Make a copy of the input molecule so we don't accidentally change it. This is
    # mainly needed due to OFF TK issue #890.
    molecule = copy_molecule(molecule)

    off_molecule = Molecule.from_openeye(molecule, allow_undefined_stereo=True)

    if max_states > 200:

        raise NotImplementedError(
            "The max states must currently be less than 200 due to a hard coded maximum "
            "value in the OpenFF toolkit"
        )

    if not enum_nitrogen:

        raise NotImplementedError(
            "Non-planer nitrogens will always be enumerated where possible."
        )

    if warts is True:
        raise NotImplementedError("The ``warts`` argument is no longer supported.")

    if verbose:
        logger.debug("Enumerating stereoisomers...")

    # Check if the input molecule has any undefined stereochemistry. If it does not
    # then it should be included in the list of returned stereoisomers.
    off_stereoisomers = off_molecule.enumerate_stereoisomers(
        undefined_only=not force_flip,
        max_isomers=max_states,
        rationalise=False,
    )

    if len(off_stereoisomers) == 0:
        # Handle the case where the input molecule does not have any stereoisomers.
        return [molecule]

    # Check to see if the original molecule had undefined stereochemistry.
    undefined_atom_stereochemistry = any(
        atom_old.stereochemistry is None and atom_new.stereochemistry is not None
        for atom_old, atom_new in zip(off_molecule.atoms, off_stereoisomers[0].atoms)
    )
    undefined_bond_stereochemistry = any(
        bond_old.stereochemistry is None and bond_new.stereochemistry is not None
        for bond_old, bond_new in zip(off_molecule.bonds, off_stereoisomers[0].bonds)
    )

    # If the input molecule had all of the input stereochemistry defined make sure to
    # return it in addition to the newly found stereoisomers.
    if not undefined_atom_stereochemistry and not undefined_bond_stereochemistry:
        off_stereoisomers = [off_molecule] + off_stereoisomers[: max_states - 1]

    return [off_stereoisomer.to_openeye() for off_stereoisomer in off_stereoisomers]
