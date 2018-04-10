__author__ = 'Chaya D. Stern'

import mdtraj as md
import os
from fnmatch import fnmatch
import sys
from math import radians
try:
    import openeye.oechem as oechem
except ImportError:
    pass
from fragmenter.utils import logger
from openmoltools.openeye import generate_conformers
import warnings
import numpy as np
from cclib.parser import Gaussian, Psi
from cclib.parser.utils import convertor
import re
import tarfile


def generate_torsions(inp_mol, output_path, interval, base_name=None, tar=True):
    """
    This function takes a 3D molecule (pdf, mol2 or sd file) and generates structures for a torsion drive on all torsions
    in the molecule. This function uses OpenEye
    Parameters
    ----------
    mol : OEMol
        molecule to generate 1D torsion scans
    output_path: str
        path to output file directory
    interval: int
        angle (in degrees) of interval for torsion drive
    base_name: str
        base name for file. Default is None. If default, use title in OEMol for base name
    tar: bool
        If true, will compress output

    """
    if not base_name:
        base_name = inp_mol.GetTitle()

    mid_tors = [[tor.a, tor.b, tor.c, tor.d ] for tor in oechem.OEGetTorsions(inp_mol)]

    # This smarts should match terminal torsions such as -CH3, -NH2, -NH3+, -OH, and -SH
    smarts = '[*]~[*]-[X2H1,X3H2,X4H3]-[#1]'
    qmol=oechem.OEQMol()
    if not oechem.OEParseSmarts(qmol, smarts):
        warnings.warn('OEParseSmarts failed')
    ss = oechem.OESubSearch(qmol)
    mol = oechem.OEMol(inp_mol)
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
        logger().info("Idxs: {} {} {} {}".format(tor[0].GetIdx(), tor[1].GetIdx(), tor[2].GetIdx(), tor[3].GetIdx()))
        logger().info("Atom Numbers: {} {} {} {}".format(tor[0].GetAtomicNum(), tor[1].GetAtomicNum(), tor[2].GetAtomicNum(), tor[3].GetAtomicNum()))
        if tor[1].GetIdx() != best_tor[1].GetIdx() or tor[2].GetIdx() != best_tor[2].GetIdx():
            new_tor = True
            if not first_pass:
                logger().info("Adding to list: {} {} {} {}".format(best_tor[0].GetIdx(), best_tor[1].GetIdx(), best_tor[2].GetIdx(), best_tor[3].GetIdx()))
                tors.append(best_tor)
            first_pass = False
            best_tor = tor
            best_tor_order = tor[0].GetAtomicNum() + tor[3].GetAtomicNum()
            logger().info("new_tor with central bond across atoms: {} {}".format(tor[1].GetIdx(), tor[2].GetIdx()))
        else:
            logger().info("Not a new_tor but now with end atoms: {} {}".format(tor[0].GetIdx(), tor[3].GetIdx()))
            tor_order = tor[0].GetAtomicNum() + tor[3].GetAtomicNum()
            if tor_order > best_tor_order:
                best_tor = tor
                best_tor_order = tor_order
    logger().info("Adding to list: {} {} {} {}".format(best_tor[0].GetIdx(), best_tor[1].GetIdx(), best_tor[2].GetIdx(), best_tor[3].GetIdx()))
    tors.append(best_tor)

    logger().info("List of torsion to drive:")
    for tor in tors:
        logger().info("Idx: {} {} {} {}".format(tor[0].GetIdx(), tor[1].GetIdx(), tor[2].GetIdx(), tor[3].GetIdx()))
        logger().info("Atom numbers: {} {} {} {}".format(tor[0].GetAtomicNum(), tor[1].GetAtomicNum(), tor[2].GetAtomicNum(), tor[3].GetAtomicNum()))

    conf = mol.GetConfs().next()
    coords = oechem.OEFloatArray(conf.GetMaxAtomIdx() * 3)
    conf.GetCoords(coords)
    # Check if coordinates are not zero
    values = np.asarray([coords.__getitem__(i) == 0 for i in range(coords.__len__())])
    if values.all():
        # Generate new coordinates.
        mol2 = generate_conformers(mol, max_confs=1)
        conf = mol2.GetConfs().next()
        coords = oechem.OEFloatArray(conf.GetMaxAtomIdx() * 3)
        conf.GetCoords(coords)
        mol2.DeleteConfs()
    mol.DeleteConfs()

    for tor in tors:
        tor_name = str((tor[0].GetIdx())+1) + '_' + str((tor[1].GetIdx())+1) + '_' + str((tor[2].GetIdx())+1) + '_' + str((tor[3].GetIdx())+1)
        folder = os.path.join(output_path, tor_name)
        try:
            os.makedirs(folder)
        except FileExistsError:
            logger().info("Overwriting existing directory {}".format(tor_name))
        for angle in range(0, 360, interval):
            angle_folder = os.path.join(folder, str(angle))
            try:
                os.mkdir(angle_folder)
            except FileExistsError:
                logger().info("Overwriting existing directory {}".format(tor_name))
            newconf = mol.NewConf(coords)
            oechem.OESetTorsion(newconf, tor[0], tor[1], tor[2], tor[3], radians(angle))
            pdb = oechem.oemolostream('{}/{}_{}_{}.pdb'.format(angle_folder, base_name, tor_name, angle))
            oechem.OEWritePDBFile(pdb, newconf)
    if tar:
        # tar archive output
        out = tarfile.open('{}.tar.gz'.format(output_path), mode='w:gz')
        os.chdir(output_path)
        os.chdir('../')
        out.add('{}'.format(base_name))
        out.close()


def pdb_to_psi4(starting_geom, mol_name, method, basis_set, charge=0, multiplicity=1, symmetry='C1', geom_opt=True,
                sp_energy=False, fixed_dih=None, mem=None, constrain='dihedral', dynamic_level=3,
                consecutive_backsteps=None, geom_maxiter=250, xyz_traj=True):
    """

    :param pdb: str
        path to pdb file
    :param method: list of str
        QM method (see psi4 website for options)
        If length 2, first one will be used for geom opt and second for spe.
    :param basis_set: str
        specification of basis set
    :param symmetry: str
        symmetry of molecule. Default is None.
    :param geom_opt: bool
        if True, will generate input file for geometry optimization
    :param sp_energy: bool
        if True, will run a single point energy calculation (if geom_opt also true, SPE calculation will occur after
        geom opt
    :param fixed_dih: str
        string of dihedral that should be fixed at specified angle. Format: "4 7 10 14 90.00"
        default: None - will not fix dihedral
        Beware:
        ------
        Because of a bug in psi4, dihedral angle can't be exactly 0 (same would apply for 180) so use 0.001 instead
    constrain: string. Either 'dihedral' or 'cartesian'
    The kind of constrain to use

    :param mem: int
        memory allocation for calculation
    :param outfile: str
        if specified, will save file there
    :return:
        psi4 input string. If outfile, save file to specified path
    """

    input_string = ""

    if mem is not None:
        input_string += "\nmemory {}\n".format(mem)

    input_string += "\nmolecule {}".format(mol_name)
    input_string += " {\n"
    input_string += "  symmetry {}\n".format(symmetry)
    input_string += "  {} {} \n".format(charge, multiplicity)

    input_string += starting_geom

    input_string += "  units Angstrom\n"
    input_string += "}\n"

    if fixed_dih is not None:
        if constrain == 'dihedral':
            input_string += '\ndih_string = "{}"'.format(fixed_dih)
            # ToDo add string because that's the only thing that seems to work
            input_string += '\nset optking { fixed_dihedral = $dih_string\n'
        elif constrain == 'cartesian':
            input_string += '\n frozen_string = """ \n {} xyz \n {} xyz \n {} xyz \n {} xyz \n"""'.format(fixed_dih[0],
                                                                                                          fixed_dih[2],
                                                                                                          fixed_dih[4],
                                                                                                          fixed_dih[6])
            input_string += '\nset optking { opt_coordinates = cartesian\n    frozen_cartesian = $frozen_string\n'
        else:
            raise NameError('Only dihedral or cartesian constraints are valid')
        if dynamic_level:
            input_string += '    dynamic_level = {}\n'.format(dynamic_level)
        if consecutive_backsteps:
            input_string += '    consecutive_backsteps = {}\n'.format(consecutive_backsteps)
        if geom_maxiter:
            input_string += '    geom_maxiter = {}\n'.format(geom_maxiter)
        if xyz_traj:
            input_string += '    print_trajectory_xyz_file = True '
        input_string += '}\n'

    if geom_opt:
        input_string += "\noptimize('{}/{}')\n".format(method[0], basis_set[0])

    if sp_energy:
        input_string += "\nenergy('{}/{}')\n".format(method[-1], basis_set[-1])

    return input_string


def generate_scan_input(root, filetype, mol_name, method, basis_set, dihedral=None, charge=0, multiplicity=1, symmetry=None,
                        geom_opt=True, sp_energy=False, mem=None):
    """
    This function takes a directory and writes out psi4 input files for all files that match the filetype specified

    :param root: str
        path to files
    :param filetype: str
        input filetypes
    :param mol_name: str
        molecule name
    :param dihedral: str
        index of atoms that should remain fixed. format '1  2  3  4'
    :param method: list of str
        QM method (see psi4 website for options)
    :param basis_set: list of str
        see psi4 website for options
    :param charge: int
        default 0
    :param multiplicity: int
        default 1
    :param symmetry: str
        symmetry of molecule. default None
    :param geom_opt: bool
        if True, run geometry optimization
    :param sp_energy: bool
        if True, run a single point energy calculation after geomoetry optimization
    :param mem: str
        memory allocation

    """
    if not dihedral:
        dihedral = list(filter(None, root.split('/')))[-1].split('_')
        dihedral = dihedral[0] + ' ' + dihedral[1] + ' ' + dihedral[2] + ' ' + dihedral[3]
    input_files = []
    pattern = "*.{}".format(filetype)
    for path, subdir, files in os.walk(root):
        for name in files:
            if fnmatch(name, pattern):
                input_files.append(os.path.join(path, name))

    for f in input_files:
        fixed_dih_angle = f.split('/')[-2]
        if fixed_dih_angle == '0':
            fixed_dih_angle = '0.001'
        if fixed_dih_angle == '180':
            fixed_dih_angle = '180.001'
        if fixed_dih_angle == '360':
            fixed_dih_angle = '360.001'
        dihedral_string = dihedral + ' ' + fixed_dih_angle
        mol = md.load(f)
        starting_geom = ""
        for i, atom in enumerate(mol.topology.atoms):
            element = atom.element.symbol
            # Convert to Angstroms
            xyz = mol.xyz[0]*10
            starting_geom += "  {}      {:05.3f}   {:05.3f}   {:05.3f}\n".format(element, xyz[i][0], xyz[i][1], xyz[i][2])

        output = pdb_to_psi4(starting_geom=starting_geom, mol_name=mol_name, method=method, basis_set=basis_set, charge=charge,
                             multiplicity=multiplicity, symmetry=symmetry, geom_opt=geom_opt, sp_energy=sp_energy,
                             fixed_dih=dihedral_string, mem=mem)

        filename = f.replace(filetype, 'dat')
        psi4_input = open(filename, 'w')
        psi4_input.write(output)
        psi4_input.close()


