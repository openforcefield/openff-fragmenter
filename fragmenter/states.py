import logging

logger = logging.getLogger(__name__)


def _enumerate_stereoisomers(
    molecule,
    max_states=200,
    force_flip=True,
    enum_nitrogen=True,
    warts=True,
    verbose=True,
):
    """
    Enumerate stereoisomers
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
    from openeye import oechem, oeomega

    stereoisomers = []
    if verbose:
        logger.debug("Enumerating stereoisomers...")
    i = 0
    for enantiomer in oeomega.OEFlipper(
        molecule, max_states, force_flip, enum_nitrogen, warts
    ):
        i += 1
        enantiomer = oechem.OEMol(enantiomer)
        stereoisomers.append(enantiomer)
    return stereoisomers
