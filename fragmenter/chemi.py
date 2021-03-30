"""functions to manipulate, read and write OpenEye and Psi4 molecules"""
import logging

import cmiles

logger = logging.getLogger(__name__)


def carboxylic_acid_hack(mol):
    """
    This is a workaround for a known bug in Openeye where carboxylic acid fails to get charged
    Parameters
    ----------
    mol : oemol

    Returns
    -------

    """
    from openeye import oechem

    # Check for Carboxylic Acid patterns in the molecule
    smarts = "(O=)[C][O,S][H]"
    ss = oechem.OESubSearch(smarts)

    oechem.OEPrepareSearch(mol, ss)
    unique_match = True

    a_match_list = []
    for match in ss.Match(mol, unique_match):

        for ma in match.GetAtoms():
            a_match_list.append(ma.target)

    # Set the Carboxylic Acid torsion to zero for each generated conformers
    if a_match_list:

        if len(a_match_list) % 4 != 0:
            raise ValueError("The atom matching list must be multiple of 4")

        for i in range(0, len(a_match_list), 4):

            chunk = a_match_list[i : i + 4]

            for conf in mol.GetConfs():
                conf.SetTorsion(chunk[0], chunk[1], chunk[2], chunk[3], 0.0)


def smiles_to_oemol(smiles, normalize=True, **kwargs):
    """Create a OEMolBuilder from a smiles string.
    Parameters
    ----------
    smiles : str
        SMILES representation of desired molecule.
    Returns
    -------
    molecule : OEMol
        A normalized molecule with desired smiles string.
    """
    from openeye import oechem

    molecule = oechem.OEMol()
    if not oechem.OESmilesToMol(molecule, smiles):
        raise ValueError("The supplied SMILES '%s' could not be parsed." % smiles)

    if normalize:
        molecule = normalize_molecule(molecule, **kwargs)

    return molecule


def normalize_molecule(molecule, name="", add_atom_map=False):
    """Normalize a copy of the molecule by checking aromaticity, adding explicit hydrogens and renaming by IUPAC name
    or given title

    Parameters
    ----------
    molecule: OEMol
        The molecule to be normalized:
    title: str
        Name of molecule. If the string is empty, will use IUPAC name

    Returns
    -------
    molcopy: OEMol
        A (copied) version of the normalized molecule

    """
    from openeye import oechem, oeiupac

    molcopy = oechem.OEMol(molecule)

    # Assign aromaticity.
    oechem.OEAssignAromaticFlags(molcopy, oechem.OEAroModelOpenEye)

    # Add hydrogens.
    oechem.OEAddExplicitHydrogens(molcopy)

    # Set title to IUPAC name.
    title = name
    if not title:
        title = oeiupac.OECreateIUPACName(molcopy)
    molcopy.SetTitle(title)

    # Check for any missing atom names, if found reassign all of them.
    if any([atom.GetName() == "" for atom in molcopy.GetAtoms()]):
        oechem.OETriposAtomNames(molcopy)

    # Add canonical ordered atom maps
    if add_atom_map:
        cmiles.utils.add_atom_map(molcopy)
    return molcopy


def get_charges(
    molecule,
    max_confs=800,
    strict_stereo=True,
    normalize=True,
    keep_confs=None,
    legacy=True,
    **kwargs
):
    """Generate charges for an OpenEye OEMol molecule.
    Parameters
    ----------
    molecule : OEMol
        Molecule for which to generate conformers.
        Omega will be used to generate max_confs conformations.
    max_confs : int, optional, default=800
        Max number of conformers to generate
    strictStereo : bool, optional, default=True
        If False, permits smiles strings with unspecified stereochemistry.
        See https://docs.eyesopen.com/omega/usage.html
    normalize : bool, optional, default=True
        If True, normalize the molecule by checking aromaticity, adding
        explicit hydrogens, and renaming by IUPAC name.
    keep_confs : int, optional, default=None
        If None, apply the charges to the provided conformation and return
        this conformation, unless no conformation is present.
        Otherwise, return some or all of the generated
        conformations. If -1, all generated conformations are returned.
        Otherwise, keep_confs = N will return an OEMol with up to N
        generated conformations.  Multiple conformations are still used to
        *determine* the charges.
    legacy : bool, default=True
        If False, uses the new OpenEye charging engine.
        See https://docs.eyesopen.com/toolkits/python/quacpactk/OEProtonFunctions/OEAssignCharges.html#
    Returns
    -------
    charged_copy : OEMol
        A molecule with OpenEye's recommended AM1BCC charge selection scheme.
    Notes
    -----
    Roughly follows
    http://docs.eyesopen.com/toolkits/cookbook/python/modeling/am1-bcc.html
    """
    from openeye import oechem, oequacpac

    # If there is no geometry, return at least one conformation.
    if molecule.GetConfs() == 0:
        keep_confs = 1

    if not oechem.OEChemIsLicensed():
        raise (ImportError("Need License for OEChem!"))
    if not oequacpac.OEQuacPacIsLicensed():
        raise (ImportError("Need License for oequacpac!"))

    if normalize:
        molecule = normalize_molecule(molecule, molecule.GetTitle())
    else:
        molecule = oechem.OEMol(molecule)

    charged_copy = generate_conformers(
        molecule, max_confs=max_confs, strict_stereo=strict_stereo, **kwargs
    )  # Generate up to max_confs conformers

    # fix issue that causes carboxylic acid to fail charging
    carboxylic_acid_hack(charged_copy)

    if not legacy:
        # 2017.2.1 OEToolkits new charging function
        status = oequacpac.OEAssignCharges(charged_copy, oequacpac.OEAM1BCCCharges())
        if not status:
            raise (RuntimeError("OEAssignCharges failed."))
    else:
        # AM1BCCSym recommended by Chris Bayly to KAB+JDC, Oct. 20 2014.
        status = oequacpac.OEAssignPartialCharges(
            charged_copy, oequacpac.OECharges_AM1BCCSym
        )
        if not status:
            raise (
                RuntimeError("OEAssignPartialCharges returned error code %d" % status)
            )

    # Determine conformations to return
    if keep_confs is None:
        # If returning original conformation
        original = molecule.GetCoords()
        # Delete conformers over 1
        for k, conf in enumerate(charged_copy.GetConfs()):
            if k > 0:
                charged_copy.DeleteConf(conf)
        # Copy coordinates to single conformer
        charged_copy.SetCoords(original)
    elif keep_confs > 0:
        logger.debug(
            "keep_confs was set to %s. Molecule positions will be reset." % keep_confs
        )

        # Otherwise if a number is provided, return this many confs if available
        for k, conf in enumerate(charged_copy.GetConfs()):
            if k > keep_confs - 1:
                charged_copy.DeleteConf(conf)
    elif keep_confs == -1:
        # If we want all conformations, continue
        pass
    else:
        # Not a valid option to keep_confs
        raise (ValueError("Not a valid option to keep_confs in get_charges."))

    return charged_copy


def generate_conformers(
    molecule,
    max_confs=800,
    dense=False,
    strict_stereo=True,
    ewindow=15.0,
    rms_threshold=1.0,
    strict_types=True,
    can_order=True,
    copy=True,
    timeout=100,
):
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
    from openeye import oechem, oeomega

    if copy:
        molcopy = oechem.OEMol(molecule)
    else:
        molcopy = molecule

    if dense:
        omega_opts = oeomega.OEOmegaOptions(oeomega.OEOmegaSampling_Dense)
        omega = oeomega.OEOmega(omega_opts)
    else:
        omega = oeomega.OEOmega()

    atom_map = False
    if cmiles.utils.has_atom_map(molcopy):
        atom_map = True
        cmiles.utils.remove_atom_map(molcopy)

    # These parameters were chosen to match http://docs.eyesopen.com/toolkits/cookbook/python/modeling/am1-bcc.html
    omega.SetMaxConfs(max_confs)
    omega.SetIncludeInput(True)
    omega.SetCanonOrder(can_order)

    omega.SetSampleHydrogens(
        True
    )  # Word to the wise: skipping this step can lead to significantly different charges!
    omega.SetEnergyWindow(ewindow)
    omega.SetRMSThreshold(
        rms_threshold
    )  # Word to the wise: skipping this step can lead to significantly different charges!

    omega.SetStrictStereo(strict_stereo)
    omega.SetStrictAtomTypes(strict_types)

    omega.SetIncludeInput(False)  # don't include input
    omega.SetMaxSearchTime(timeout)
    if max_confs is not None:
        omega.SetMaxConfs(max_confs)

    status = omega(molcopy)  # generate conformation
    if not status:
        raise (RuntimeError("omega returned error code %d" % status))

    if atom_map:
        cmiles.utils.restore_atom_map(molcopy)

    return molcopy
