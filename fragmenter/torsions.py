__author__ = 'Chaya D. Stern'

try:
    import openeye.oechem as oechem
except ImportError:
    pass
from fragmenter.utils import logger
import warnings
import numpy as np

warnings.simplefilter('always')


def create_mapped_smiles(molecule):
    """
    Generate an index-tagged explicit hydrogen SMILES.
    Exmaple:
    SMILES string for carbon monoxide "CO"
    With index-tagged explicit hydrogen SMILES this becomes
    '[H:3][C:1]([H:4])([H:5])[O:2][H:6]'

    Parameters
    ----------
    molecule: OEMOl

    Returns
    -------
    index-tagged explicit hydrogen SMILES str

    """
    # Check if molecule already has explicit hydrogens
    HAS_HYDROGENS = oechem.OEHasExplicitHydrogens(molecule)
    if not HAS_HYDROGENS:
        # Add explicit hydrogens
        oechem.OEAddExplicitHydrogens(molecule)
    for atom in molecule.GetAtoms():
        atom.SetMapIdx(atom.GetIdx() + 1)

    # add tag to data
    tag = oechem.OEGetTag("has_map")
    molecule.SetData(tag, bool(True))

    return oechem.OEMolToSmiles(molecule)


def find_torsions(molecule):
    """
    This function takes an OEMol (atoms must be tagged with index map) and finds the map indices for torsion that need
    to be driven.
    Parameters
    ----------
    molecule : OEMol
        The atoms in the molecule need to be tagged with map indices

    """
    # Check if molecule has map
    try:
        molecule.GetData('has_map')
    except ValueError:
        warnings.warn('Molecule does not have atom map. A new map will be generated. You might need a new tagged SMARTS if the ordering was changed')
        tagged_smiles = create_mapped_smiles(molecule)
        logger().info('If you already have a tagged SMARTS, compare it with the new one to ensure the ordering did not change')
        logger().info('The new tagged SMARTS is: {}'.format(tagged_smiles))

    mid_tors = [[tor.a, tor.b, tor.c, tor.d ] for tor in oechem.OEGetTorsions(molecule)]

    # This smarts should match terminal torsions such as -CH3, -NH2, -NH3+, -OH, and -SH
    smarts = '[*]~[*]-[X2H1,X3H2,X4H3]-[#1]'
    qmol=oechem.OEQMol()
    if not oechem.OEParseSmarts(qmol, smarts):
        warnings.warn('OEParseSmarts failed')
    ss = oechem.OESubSearch(qmol)
    mol = oechem.OEMol(molecule)
    h_tors = []
    oechem.OEPrepareSearch(mol, ss)
    unique = True
    for match in ss.Match(mol, unique):
        tor = []
        for ma in match.GetAtoms():
            tor.append(ma.target)
        h_tors.append(tor)

    # Combine middle and terminal torsions
    all_tors = mid_tors + h_tors
    # Sort all_tors so that it's grouped by central bond
    central_bonds = np.zeros((len(all_tors), 3), dtype=int)
    for i, tor in enumerate(all_tors):
        central_bonds[i][0] = i
        central_bonds[i][1] = tor[1].GetIdx()
        central_bonds[i][2] = tor[2].GetIdx()

    grouped = central_bonds[central_bonds[:, 2].argsort()]
    sorted_tors = [all_tors[i] for i in grouped[:, 0]]

    # Keep only one torsion per rotatable bond
    tors = []
    best_tor = [sorted_tors[0][0], sorted_tors[0][0], sorted_tors[0][0], sorted_tors[0][0]]
    first_pass = True
    for tor in sorted_tors:
        logger().info("Map Idxs: {} {} {} {}".format(tor[0].GetMapIdx(), tor[1].GetMapIdx(), tor[2].GetMapIdx(), tor[3].GetMapIdx()))
        logger().info("Atom Numbers: {} {} {} {}".format(tor[0].GetAtomicNum(), tor[1].GetAtomicNum(), tor[2].GetAtomicNum(), tor[3].GetAtomicNum()))
        if tor[1].GetMapIdx() != best_tor[1].GetMapIdx() or tor[2].GetMapIdx() != best_tor[2].GetMapIdx():
            new_tor = True
            if not first_pass:
                logger().info("Adding to list: {} {} {} {}".format(best_tor[0].GetMapIdx(), best_tor[1].GetMapIdx(), best_tor[2].GetMapIdx(), best_tor[3].GetMapIdx()))
                tors.append(best_tor)
            first_pass = False
            best_tor = tor
            best_tor_order = tor[0].GetAtomicNum() + tor[3].GetAtomicNum()
            logger().info("new_tor with central bond across atoms: {} {}".format(tor[1].GetMapIdx(), tor[2].GetMapIdx()))
        else:
            logger().info("Not a new_tor but now with end atoms: {} {}".format(tor[0].GetMapIdx(), tor[3].GetMapIdx()))
            tor_order = tor[0].GetAtomicNum() + tor[3].GetAtomicNum()
            if tor_order > best_tor_order:
                best_tor = tor
                best_tor_order = tor_order
    logger().info("Adding to list: {} {} {} {}".format(best_tor[0].GetMapIdx(), best_tor[1].GetMapIdx(), best_tor[2].GetMapIdx(), best_tor[3].GetMapIdx()))
    tors.append(best_tor)

    logger().info("List of torsion to drive:")
    for tor in tors:
        logger().info("Idx: {} {} {} {}".format(tor[0].GetIdx(), tor[1].GetIdx(), tor[2].GetIdx(), tor[3].GetIdx()))
        logger().info("Atom numbers: {} {} {} {}".format(tor[0].GetAtomicNum(), tor[1].GetAtomicNum(), tor[2].GetAtomicNum(), tor[3].GetAtomicNum()))

    needed_torsion_scans = dict()
    for i, tor in enumerate(tors):
        tor_name = ((tor[0].GetIdx()+1), (tor[1].GetIdx()+1), (tor[2].GetIdx()+1), (tor[3].GetIdx()+1))
        needed_torsion_scans['torsion_{}'.format(str(i))] = tor_name
    return needed_torsion_scans


# def pdb_to_psi4(starting_geom, mol_name, method, basis_set, charge=0, multiplicity=1, symmetry='C1', geom_opt=True,
#                 sp_energy=False, fixed_dih=None, mem=None, constrain='dihedral', dynamic_level=3,
#                 consecutive_backsteps=None, geom_maxiter=250, xyz_traj=True):
#     """
#
#     :param pdb: str
#         path to pdb file
#     :param method: list of str
#         QM method (see psi4 website for options)
#         If length 2, first one will be used for geom opt and second for spe.
#     :param basis_set: str
#         specification of basis set
#     :param symmetry: str
#         symmetry of molecule. Default is None.
#     :param geom_opt: bool
#         if True, will generate input file for geometry optimization
#     :param sp_energy: bool
#         if True, will run a single point energy calculation (if geom_opt also true, SPE calculation will occur after
#         geom opt
#     :param fixed_dih: str
#         string of dihedral that should be fixed at specified angle. Format: "4 7 10 14 90.00"
#         default: None - will not fix dihedral
#         Beware:
#         ------
#         Because of a bug in psi4, dihedral angle can't be exactly 0 (same would apply for 180) so use 0.001 instead
#     constrain: string. Either 'dihedral' or 'cartesian'
#     The kind of constrain to use
#
#     :param mem: int
#         memory allocation for calculation
#     :param outfile: str
#         if specified, will save file there
#     :return:
#         psi4 input string. If outfile, save file to specified path
#     """
#
#     input_string = ""
#
#     if mem is not None:
#         input_string += "\nmemory {}\n".format(mem)
#
#     input_string += "\nmolecule {}".format(mol_name)
#     input_string += " {\n"
#     input_string += "  symmetry {}\n".format(symmetry)
#     input_string += "  {} {} \n".format(charge, multiplicity)
#
#     input_string += starting_geom
#
#     input_string += "  units Angstrom\n"
#     input_string += "}\n"
#
#     if fixed_dih is not None:
#         if constrain == 'dihedral':
#             input_string += '\ndih_string = "{}"'.format(fixed_dih)
#             # ToDo add string because that's the only thing that seems to work
#             input_string += '\nset optking { fixed_dihedral = $dih_string\n'
#         elif constrain == 'cartesian':
#             input_string += '\n frozen_string = """ \n {} xyz \n {} xyz \n {} xyz \n {} xyz \n"""'.format(fixed_dih[0],
#                                                                                                           fixed_dih[2],
#                                                                                                           fixed_dih[4],
#                                                                                                           fixed_dih[6])
#             input_string += '\nset optking { opt_coordinates = cartesian\n    frozen_cartesian = $frozen_string\n'
#         else:
#             raise NameError('Only dihedral or cartesian constraints are valid')
#         if dynamic_level:
#             input_string += '    dynamic_level = {}\n'.format(dynamic_level)
#         if consecutive_backsteps:
#             input_string += '    consecutive_backsteps = {}\n'.format(consecutive_backsteps)
#         if geom_maxiter:
#             input_string += '    geom_maxiter = {}\n'.format(geom_maxiter)
#         if xyz_traj:
#             input_string += '    print_trajectory_xyz_file = True '
#         input_string += '}\n'
#
#     if geom_opt:
#         input_string += "\noptimize('{}/{}')\n".format(method[0], basis_set[0])
#
#     if sp_energy:
#         input_string += "\nenergy('{}/{}')\n".format(method[-1], basis_set[-1])
#
#     return input_string
#
#
# def generate_scan_input(root, filetype, mol_name, method, basis_set, dihedral=None, charge=0, multiplicity=1, symmetry=None,
#                         geom_opt=True, sp_energy=False, mem=None):
#     """
#     This function takes a directory and writes out psi4 input files for all files that match the filetype specified
#
#     :param root: str
#         path to files
#     :param filetype: str
#         input filetypes
#     :param mol_name: str
#         molecule name
#     :param dihedral: str
#         index of atoms that should remain fixed. format '1  2  3  4'
#     :param method: list of str
#         QM method (see psi4 website for options)
#     :param basis_set: list of str
#         see psi4 website for options
#     :param charge: int
#         default 0
#     :param multiplicity: int
#         default 1
#     :param symmetry: str
#         symmetry of molecule. default None
#     :param geom_opt: bool
#         if True, run geometry optimization
#     :param sp_energy: bool
#         if True, run a single point energy calculation after geomoetry optimization
#     :param mem: str
#         memory allocation
#
#     """
#     if not dihedral:
#         dihedral = list(filter(None, root.split('/')))[-1].split('_')
#         dihedral = dihedral[0] + ' ' + dihedral[1] + ' ' + dihedral[2] + ' ' + dihedral[3]
#     input_files = []
#     pattern = "*.{}".format(filetype)
#     for path, subdir, files in os.walk(root):
#         for name in files:
#             if fnmatch(name, pattern):
#                 input_files.append(os.path.join(path, name))
#
#     for f in input_files:
#         fixed_dih_angle = f.split('/')[-2]
#         if fixed_dih_angle == '0':
#             fixed_dih_angle = '0.001'
#         if fixed_dih_angle == '180':
#             fixed_dih_angle = '180.001'
#         if fixed_dih_angle == '360':
#             fixed_dih_angle = '360.001'
#         dihedral_string = dihedral + ' ' + fixed_dih_angle
#         mol = md.load(f)
#         starting_geom = ""
#         for i, atom in enumerate(mol.topology.atoms):
#             element = atom.element.symbol
#             # Convert to Angstroms
#             xyz = mol.xyz[0]*10
#             starting_geom += "  {}      {:05.3f}   {:05.3f}   {:05.3f}\n".format(element, xyz[i][0], xyz[i][1], xyz[i][2])
#
#         output = pdb_to_psi4(starting_geom=starting_geom, mol_name=mol_name, method=method, basis_set=basis_set, charge=charge,
#                              multiplicity=multiplicity, symmetry=symmetry, geom_opt=geom_opt, sp_energy=sp_energy,
#                              fixed_dih=dihedral_string, mem=mem)
#
#         filename = f.replace(filetype, 'dat')
#         psi4_input = open(filename, 'w')
#         psi4_input.write(output)
#         psi4_input.close()


