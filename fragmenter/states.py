from cmiles.utils import mol_to_smiles

from .utils import logger


def enumerate_states(molecule, tautomers=True, stereoisomers=True, verbose=False, return_mols=False,
                     explicit_h=True, return_names=False, max_stereo_returns=1, filter_nitro=True, **kwargs):
    """
    Expand tautomeric state and stereoisomers for molecule.

    Parameters
    ----------
    molecule : OEMol
        Molecule to enumerate states
    tautomers : bool, optional, default True
        If False, will not generate tautomers
    stereoisomers : bool, optional, default True
        If False, will not generate all stereoisomers.
    verbose : bool, optional, default False
        If True, output will be verbose
    return_mols : bool, optional, default False
        If True, will return oemols instead of SMILES. Some molecules might be duplicate states
    explicit_h : bool, optional, default True
        If True, SMILES of states will have explicit hydrogen
    return_names : bool, optional, default True
        If True, will return names of molecules with SMILES
    max_stereo_returns : int, optional, default 1
        If stereoisomers is set to False, and the incoming molecule is missing stereo information, OEFlipper will
        generate stereoisomers for missing stereo center. max_stereo_returns controls how many of those will be returned
    ** max_states: int, optional, default 200
        This gets passed to `_enumerate_tautomers` and `_enumerate_stereoisomers`
        max number of states `_enumerate_tautomers` and `_enumerate_stereoisomers` generate
    ** pka_norm: bool, optional, default True
        This gets passed to `_enumerate_tautomers`. If True, ionization state of each tautomer will be assigned to a predominate
        state at pH ~7.4
    ** warts: bool, optional, default True
        This gets passed to `_enumerate_tautomers` and _enumerate_stereoisomers`
        If True, adds a wart to each new state. A 'wart' is a systematic
    ** force_flip: bool, optional, default True
        This gets passed to `_enumerate_stereoisomers`
        Force flipping all stereocenters. If False, will only generate stereoisomers for stereocenters that are undefined
    ** enum_nitorgen: bool, optional, default True
        This gets passed to `_enumerate_stereoisomers`
        If true, invert non-planer nitrogens

    Returns
    -------
    states: list
        list of oemols or SMILES of states generated for molecule

    """
    from openeye import oechem

    # If incoming molecule has nitro in form ([NX3](=O)=O), do not filter out later
    if _check_nitro(molecule):
        filter_nitro = False
    title = molecule.GetTitle()
    states = []
    if return_names:
        names = []

    if verbose:
        logger().info("Enumerating states for {}".format(title))

    if stereoisomers:
        if verbose:
            logger().info("Enumerating stereoisomers for {}".format(title))
        stereo_mols = (_enumerate_stereoisomers(molecule, **kwargs))
        if verbose:
            logger().info('Enumerated {} stereoisomers'.format(len(stereo_mols)))

    if tautomers:
        if not stereoisomers:
            stereo_mols = [molecule]
        tau_mols = []
        if verbose:
            logger().info("Enumerating tautomers states for {}".format(title))
        for mol in stereo_mols:
            tau_mols.extend(_enumerate_tautomers(mol, **kwargs))
        if verbose:
            logger().info('Enumerated {} tautomers'.format(len(tau_mols)))

        # check for nitro in ([NX3](=O)=O) form
        if filter_nitro:
            tau_mols[:] = [mol for mol in tau_mols if not _check_nitro(mol)]

    if stereoisomers and tautomers:
        all_mols = stereo_mols + tau_mols
    elif stereoisomers and not tautomers:
        all_mols = stereo_mols
    elif not stereoisomers and tautomers:
        all_mols = tau_mols
        all_mols.append(molecule)
    else:
        all_mols = [molecule]

    if return_mols:
        return all_mols

    for mol in all_mols:
        try:
            smiles = mol_to_smiles(mol, isomeric=True, mapped=False, explicit_hydrogen=explicit_h)
            if smiles not in states:
                states.append(smiles)
                if return_names:
                    names.append(mol.GetTitle())

        except ValueError:
            # Stereo is not fully defined. Use flipper with force_flip set to False
            stereo_states = _enumerate_stereoisomers(mol, force_flip=False, enum_nitrogen=True, warts=True)
            if len(stereo_states) > max_stereo_returns:
                stereo_states = stereo_states[:max_stereo_returns]

            for state in stereo_states:
                try:
                    smiles = mol_to_smiles(state, isomeric=True, mapped=False, explicit_hydrogen=explicit_h)
                except ValueError:
                    stereo_states_forced = _enumerate_stereoisomers(mol, force_flip=True, enum_nitrogen=True, warts=True)
                    if len(stereo_states_forced) > max_stereo_returns:
                        stereo_states_forced = stereo_states_forced[:max_stereo_returns]
                    for state_forced in stereo_states_forced:
                        smiles = mol_to_smiles(state_forced, isomeric=True, mapped=False, explicit_hydrogen=explicit_h)
                        if smiles not in states:
                            states.append(smiles)
                            if return_names:
                                names.append(state.GetTitle())
                if smiles not in states:
                    states.append(smiles)
                    if return_names:
                        names.append(state.GetTitle())

    if verbose:
        logger().info("{} states were generated for {}".format(len(states), oechem.OEMolToSmiles(molecule)))

    if return_names:
        return states, names

    return states

def _enumerate_tautomers(molecule, max_states=200, pka_norm=True, warts=True):
    """
    Expand reasonable tautomer states. This function generates tautomers (which might be different ionization states
    than parent) that are normalized to the predominant state at pH ~7.4
    Parameters
    ----------
    molecule : OEMol to expand states
    max_states : int
        max number of states
    pka_norm: bool, optional, default True
    warts: bool, optional default True

    Returns
    -------
    tautomers: list of oemols

    """
    from openeye import oequacpac

    tautomers = []
    tautomer_options = oequacpac.OETautomerOptions()
    tautomer_options.SetApplyWarts(warts)
    tautomer_options.SetMaxTautomersGenerated(max_states)
    i = 0
    for tautomer in oequacpac.OEGetReasonableTautomers(molecule, tautomer_options, pka_norm):
        i += 1
        tautomers.append(tautomer)
    return tautomers

def _enumerate_stereoisomers(molecule, max_states=200, force_flip=True, enum_nitrogen=True, warts=True, verbose=True):
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
    from openeye import oeomega, oechem
    stereoisomers = []
    if verbose:
        logger().debug("Enumerating stereoisomers...")
    i = 0
    for enantiomer in oeomega.OEFlipper(molecule, max_states, force_flip, enum_nitrogen, warts):
        i += 1
        enantiomer = oechem.OEMol(enantiomer)
        stereoisomers.append(enantiomer)
    return stereoisomers

def _check_nitro(molecule):
    """
    Filter out nitro that is in ([NX3](=O)=O) form. OEGetReasonableTautomers generates this form.
    Parameters
    ----------
    molecule :

    Returns
    -------

    """
    from openeye import oechem

    qmol = oechem.OEQMol()
    if not oechem.OEParseSmarts(qmol, '([NX3](=O)=O)'):
        print('OEParseSmarts failed')
    ss = oechem.OESubSearch(qmol)
    oechem.OEPrepareSearch(molecule, ss)
    matches = [m for m in ss.Match(molecule)]
    return bool(matches)
