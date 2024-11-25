import os
import math 
import shutil
import json
from glob import glob
import warnings
from parmed.charmm import CharmmParameterSet
from typing import Dict, List, Tuple

def get_boundary_bonds(Q1_coords:Tuple[str,str,str], M1_coords:Tuple[str,str,str], MOL2_PATH:str, ff_type:str, params:CharmmParameterSet = None) -> Tuple[str,str,str,str]: 
    """
    Given the coordinates of the two atoms in the fronteir bond (the
    Q1 atom and M1 atom), returns the atom types in the fronteir bond
    as well as atom types for the new Q1-capping H bond

    This is necessary because the indexing of the mol2 and pdb do
    not match

    Parameters
    ----------
    Q1_coords: Tuple[str,str,str]
        coordiantes of the Q1 atom that is being capped

    M1_coords: Tuple[str, str, str]
        coordinates of the M1 atom that was attached to the Q1 atom

    MOL2_PATH: str
        path to mol2

    ff_type: str
        forcefield type, charmm or amber

    params: CharmmParameterSet or None
        parmed object containing information from the charmm topology and parameters

    Returns
    -------
    QM_bond, QM_bond_perm, QH_bond, QH_bond_perm: Tuple[str,str,str,str]
        bond information for the old bond containing the Q1 atom and M1 atom, 
        a permutation of that bond, the new Q1-capping H bond, and a permutation of that bond

    """
    out = open(glob('*.out')[0], 'a')
    # get atom types of Q1 and M1
    with open(MOL2_PATH, 'r', encoding="iso-8859-1") as mol2file:
        lines = mol2file.readlines()
        for n,l in enumerate(lines):                                        # filter lines for atoms only
            if 'ATOM' in l:
                start = n
            elif 'BOND' in l:
                end = n
        atom_lines = lines[start+1:end]
        for l in atom_lines:
            x = float("{:.3f}".format(float(l.split()[2])))
            y = float("{:.3f}".format(float(l.split()[3])))
            z = float("{:.3f}".format(float(l.split()[4])))
            if float(Q1_coords[0]) == x and float(Q1_coords[1]) == y and float(Q1_coords[2]) == z:
                Q1_at = l.split()[5]
            elif float(M1_coords[0]) == x and float(M1_coords[1]) == y and float(M1_coords[2]) == z:
                M1_at = l.split()[5]
    QM_bond = f'{Q1_at:<2}' + '-' + f'{M1_at:<2}'
    QM_bond_perm = f'{M1_at:<2}' + '-' + f'{Q1_at:<2}' #X-Y bond is the same as Y-X bond
    if ff_type == 'amber':
        QH_bond = f'XC' + '-' + f"{'H1':<2}"
        QH_bond_perm = f"{'H1':<2}" + '-' + f'XC'
    else:
        try:
            QH_bond = f'CT2-HB2'
            QH_bond_perm = f"HB2-CT2"
            atom1 = QH_bond.split('-')[0].strip()
            atom2 = QH_bond.split('-')[1].strip()
            test = params.bond_types[(atom1, atom2)].req
        except:
            QH_bond = f'CT2-HB'
            QH_bond_perm = f"HB-CT2"
            atom1 = QH_bond.split('-')[0].strip()
            atom2 = QH_bond.split('-')[1].strip()
            test = params.bond_types[(atom1, atom2)].req
    out.close()
    return QM_bond, QM_bond_perm, QH_bond, QH_bond_perm


def calc_C_HL(QM_bond:str, QM_bond_perm:str, QH_bond:str, QH_bond_perm:str, ff_type:str, params:CharmmParameterSet = None, PARM_PATH:str= None) -> float:
    """
    This function returns the scaling factor for a given cut bond based
    on the forcefield bond sretch parameters of the Q1-M1 bond and the
    new Q1-capping H bond

    Parameters
    ----------
    QM_bond:str
        atoms in the old bond that was cut with one atom in the QM region 
        and one atom in the MM region

    QM_bond_perm:str
        permutation of QM_bond

    QH_bond:str
        atoms in the new bond that is being capped with one atom in QM region 
        and one atom is the new capping hydrogen 

    QH_bond_perm:str
        permutation of QH_bond

    ff_type: str
        forcefield type, charmm or amber

    params: CharmmParameterSet or None
        parmed object containing information from the charmm topology and parameters

    PARM_PATH: str or None
        path to parameters for amber forcefield

    Returns
    -------
    R0_Q1H1/R0_Q1M1: float
        length of scaling factor
    """
    out = open(glob('*.out')[0], 'a')
    # get bond parameters and divide them
    if ff_type =='amber':
        with open(PARM_PATH, 'r') as parmfile:
            lines = parmfile.readlines()
            for l in lines:
                if l.startswith(QM_bond+'  ') or l.startswith(QM_bond_perm+'  '):
                    line = l.split('  ')
                    R0_Q1M1 = float(line[3])
                if l.startswith(QH_bond+'  ') or l.startswith(QH_bond_perm+'  '):
                    line = l.split('  ')
                    R0_Q1HL = float(line[3])
    else:
        atom1 = QM_bond.split('-')[0].strip()
        atom2 = QM_bond.split('-')[1].strip()
        R0_Q1M1 = params.bond_types[(atom1, atom2)].req

        atom1 = QH_bond.split('-')[0].strip()
        atom2 = QH_bond.split('-')[1].strip()
        R0_Q1HL = params.bond_types[(atom1, atom2)].req
    out.write(f'MM bond distance parameter (R0) for {QH_bond} stretch: {R0_Q1HL}\n')
    out.write(f'MM bond distance parameter (R0) for {QM_bond} stretch: {R0_Q1M1}\n')
    out.close() 
    return R0_Q1HL/R0_Q1M1

def cap(no_HL:Dict[str,List[int]], num_broken_bonds:int, PDB_PATH:str, MOL2_PATH:str, CAPPED_PDB_PATH:str, ff_type:str, params:CharmmParameterSet = None, PARM_PATH:str=None) -> None:
    """
    caps Q1 atoms with hydrogen, scaled along Q1 M1 bond; produces new dictionary 'with_HL'

    Parameters
    ----------
    no_HL: Dict[str, List[int]]
        Dictionary containing keys for each region and values of a list of each atom in that region

    num_broken_bonds: int
        number of bonds that were broken in cut_protein and need to be capped with hydrogens on the QM side

    PDB_PATH: str
        path to cx PDB

    MOL2_PATH: str
        path to mol2

    CAPPED_PDB_PATH: str
        path to the pdb containing the new capping atoms

    ff_type: str
        forcefield type, charmm or amber

    params: CharmmParameterSet or None
        parmed object containing information from the charmm topology and parameters

    PARM_PATH: str or None
        path to parameters for amber forcefield

    Returns
    -------
    None
    """
    out = open(glob('*.out')[0], 'a')
    with open(PDB_PATH, 'r') as pdbfile:
        lines = pdbfile.readlines()
        lines = [l for l in lines if l.startswith('ATOM')]
        lines = ['0']+lines                           # so that python indexing matches atom indexing
    for bond in range(1,num_broken_bonds+1):                           # loop through broken bonds/frontier regions
        out.write(f'Capping Q1 atom of bond {bond}:\n')
        Q1_atom = no_HL[f'Q1_{bond}']
        out.write(f'Q1 atom ID: {Q1_atom}\n')
        M1_atom = no_HL[f'M1_{bond}']
        out.write(f'M1 atom ID: {M1_atom}\n')
        for line in lines:
            if len(list(line.split())) > 5:
                x_coord = line[29:38].strip()
                y_coord = line[38:46].strip()
                z_coord = line[46:54].strip()
                if line[6:11].strip() == str(Q1_atom[0]):
                    Q1_coords = [x_coord, y_coord, z_coord]
                    Q1_pdb_atom_type = line[11:16].strip()
                if line[6:11].strip() == str(M1_atom[0]):
                    M1_coords = [x_coord, y_coord, z_coord]
                    M1_pdb_atom_type = line[11:16].strip()
        out.write(f'Q1 atom type: {Q1_pdb_atom_type}\n')
        out.write(f'M1 atom type: {M1_pdb_atom_type}\n')
        Q1_coords = [float(x) for x in Q1_coords]
        M1_coords = [float(x) for x in M1_coords]
        R_Q1M1 = math.dist(Q1_coords, M1_coords)
        a,b,c,d = get_boundary_bonds(Q1_coords, M1_coords, MOL2_PATH, ff_type = ff_type, params = params)
        out.close()
        out = open(glob('*.out')[0], 'a')
        if ff_type =='amber':
            C_HL = calc_C_HL(a,b,c,d, ff_type, PARM_PATH = PARM_PATH) #doi:10.1021/jp0446332 eqn (4)
        else:
            C_HL = calc_C_HL(a,b,c,d, ff_type, params = params)                             # doi:10.1021/jp0446332 eqn (4)
        out.write(f'C_HL: {C_HL:.4f}\n')
        R_Q1HL = C_HL * R_Q1M1                                          # doi:10.1021/jp0446332 eqn (3)
        out.write(f'R_Q1HL(Å)={R_Q1HL:.4f}\n\n')
        R_ratio = (R_Q1HL / R_Q1M1)
        HL_coords = []
        for n,x in enumerate(Q1_coords):                                # calc coordinates for HL
            HL_coords.append("{:.3f}".format((1-R_ratio)*x+R_ratio*M1_coords[n]))
        if os.path.exists(CAPPED_PDB_PATH) == False:
            shutil.copyfile(PDB_PATH, CAPPED_PDB_PATH)
        with open(CAPPED_PDB_PATH, 'r') as cappedfile:    
            cap_lines = cappedfile.readlines()
            CONECT = []
            for n, l in enumerate(cap_lines):
                if l.startswith('ATOM'):
                    idx_last_ATOM_line = n
                elif l.startswith('CONECT'):
                    CONECT.append(int(n))
            idx_last_HETATM = min(CONECT)
            HL_atom_num = int(cap_lines[int(idx_last_HETATM)-1][6:11].strip())
            HL_atom_num = HL_atom_num + bond
            HL_bond =f'HL{bond}'
            HL_line =f"{'ATOM':<6}"+f"{HL_atom_num:>5}"+f"{HL_bond:>5} "+'LIN K '+ f"{bond:03d}"
            HL_line += f"{HL_coords[0]:>12}"+f"{HL_coords[1]:>8}"+f"{HL_coords[2]:>8}"
            HL_line += f"{'1.00':>6}"+f"{'0.00':>6}"+f"{'H':>12}"+'\n'
            cap_lines.insert(idx_last_ATOM_line+1, HL_line)
        with open(CAPPED_PDB_PATH, 'w') as cappedfile:                 # write new PDB to add link H
            cappedfile.writelines(cap_lines)
        if 'with_HL' in globals() or 'with_HL' in locals():
            with_HL[f'HL_{bond}'] = [HL_atom_num]
        else:
            with_HL = no_HL.copy()
            with_HL[f'HL_{bond}'] = [HL_atom_num]
    with open('with_HL.dat', 'w+') as wfile:
        json.dump(with_HL, wfile)
    out.close()
    return 

def run_cap(ff_type:str, rtf:str = None, prm:str = None) -> None:    
    """
    caps Q1 atoms with hydrogen, scaled along Q1 M1 bond; produces new dictionary 'with_HL'

    Parameters
    ----------
    ff_type: str
        forcefield type (charmm or amber)

    rtf: str or None
        path to topology for charmm forcefield. required if ff_type is charmm to look up ff parameters

    prm: str or None
        path to paramteres for charmm forcefield. required if ff_type is charmm to look up ff parameters

    Returns
    -------
    None
    """
    out = open(glob('*.out')[0], 'a')
    out.write('----------------------------------------------------------------------------------------------------\n')
    out.write('cap'.center(100)+'\n')
    out.write('----------------------------------------------------------------------------------------------------\n')
    out.write('Atom IDs refer to cx_autocap_fixed.pdb\n')
    out.close()
        #out.write('C_HL should be about 0.71. (doi:https://doi.org/10.1021/jp0446332\n')
    # load json dictionary
    with open('dictionary.dat', 'r') as dictfile:
        no_HL = json.load(dictfile)
    
    num_broken_bonds = int(list(no_HL.keys())[-1].split('_')[-1])
    
    PDB_PATH = 'cx_autocap_fixed.pdb' 
    MOL2_PATH = 'prot_autocap_fixed.mol2'
    CAPPED_PDB_PATH = 'CAPPED-prot_autocap_fixed.pdb'
        
    if ff_type =='amber':
        path_to_env = os.environ['AMBERHOME']
        PARM_PATH = f'{path_to_env}/dat/leap/parm/parm19.dat' 
        cap(no_HL, num_broken_bonds, PDB_PATH, MOL2_PATH, CAPPED_PDB_PATH, ff_type ='amber', PARM_PATH = PARM_PATH)
    else:
        if os.path.isabs(rtf):
            rtf_file = rtf
        else:
            rtf_file = '../' + rtf
        
        if os.path.isabs(prm):
            prm_file = prm
        else:
            prm_file = '../' + prm
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            params = CharmmParameterSet(prm_file,rtf_file)
        cap(no_HL, num_broken_bonds, PDB_PATH, MOL2_PATH, CAPPED_PDB_PATH, 'charmm', params = params)
    out = open(glob('*.out')[0], 'a')
    out.write(f'Hydrogen link atoms to cap the QM region have been added to {CAPPED_PDB_PATH}\n')
    out.close()
