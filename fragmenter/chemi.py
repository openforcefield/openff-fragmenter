"""functions to manipulate OpenEye and RDKit molecules"""

try:
    from openeye import oechem, oeomega
except ImportError:
    raise Warning("Need license for OpenEye!")


def generate_conformers(molecule, max_confs=800, strict_stereo=True, ewindow=15.0, rms_threshold=1.0, strict_types=True,
                        copy=True):
    """Generate conformations for the supplied molecule
    Parameters
    ----------
    molecule : OEMol
        Molecule for which to generate conformers
    max_confs : int, optional, default=800
        Max number of conformers to generate.  If None, use default OE Value.
    strict_stereo : bool, optional, default=True
        If False, permits smiles strings with unspecified stereochemistry.
    strict_types : bool, optional, default=True
        If True, requires that Omega have exact MMFF types for atoms in molecule; otherwise, allows the closest atom
        type of the same element to be used.
    Returns
    -------
    molcopy : OEMol
        A multi-conformer molecule with up to max_confs conformers.
    Notes
    -----
    Roughly follows
    http://docs.eyesopen.com/toolkits/cookbook/python/modeling/am1-bcc.html
    """
    if copy:
        molcopy = oechem.OEMol(molecule)
    else:
        molcopy = molecule
    omega = oeomega.OEOmega()

    # These parameters were chosen to match http://docs.eyesopen.com/toolkits/cookbook/python/modeling/am1-bcc.html
    omega.SetMaxConfs(max_confs)
    omega.SetIncludeInput(True)
    omega.SetCanonOrder(False)

    omega.SetSampleHydrogens(True)  # Word to the wise: skipping this step can lead to significantly different charges!
    omega.SetEnergyWindow(ewindow)
    omega.SetRMSThreshold(rms_threshold)  # Word to the wise: skipping this step can lead to significantly different charges!

    omega.SetStrictStereo(strict_stereo)
    omega.SetStrictAtomTypes(strict_types)

    omega.SetIncludeInput(False)  # don't include input
    if max_confs is not None:
        omega.SetMaxConfs(max_confs)

    status = omega(molcopy)  # generate conformation
    if not status:
        raise(RuntimeError("omega returned error code %d" % status))

    return molcopy
