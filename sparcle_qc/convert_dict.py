import os
import json
import sys

from pymol.cgo import *
from pymol import cmd
from glob import glob
import os
from typing import List
from typing import Dict

def match_resi_neighborhood(Me_PDB_lines: List[List[str]], Cl_PDB_lines: List[List[str]]) -> Dict[str,str]:
    """
    This function will use the 'neighborhood' of a given residue
    in one PDB (called Me) to match it to the residue in another
    PDB (called Cl).  

    Inputs are the PDB lines from both.  

    Because the resnum may not be the same between the pdbs, we consider
    the protein sequence to match the neighborhoods of the residues.
    This is necessary if the QM region is desired to be the same between
    two similar protein:ligand complexes, but there are small differences
    between protonation states or residue numbers are different.

    Output will be a dictionary that has the cl residue as the key and
    the corresponding me residue as the value, where 'residue' is the
    3 letter code + residue number (eg LEU251).  

    Parameters
    ----------
    Me_PDB_Lines: List(List(str))
        List containing the lines of the reference PDB where each line is another list containing each category from the PDB
    Cl_PDB_lines: List(List(str))
        List containing the lines of the PDB to map where each line is another list containing each category from the PDB
    Returns
    -------
    mapping_dict: Dict[str, str]
        the key contains the residue in the reference PDB and the value is the corresponding residue in the PDB to map
    """
    me_resis = []
    cl_resis = []
    mapping_dict = {}
    for i in range(len(Me_PDB_lines)):
        if [Me_PDB_lines[i][3][:-1], Me_PDB_lines[i][5]] not in me_resis:
            me_resis.append([Me_PDB_lines[i][3][:-1], Me_PDB_lines[i][5]])
    for i in range(len(Cl_PDB_lines)):
        if [Cl_PDB_lines[i][3][:-1], Cl_PDB_lines[i][5]] not in cl_resis:
            cl_resis.append([Cl_PDB_lines[i][3][:-1], Cl_PDB_lines[i][5]])
    for me_num, me_resi in enumerate(me_resis):
        for cl_num,cl_resi in enumerate(cl_resis):
            true_count = 0
            for i in range(0,5):
                try:
                    if cl_resis[cl_num+i][0]==me_resis[me_num+i][0]:
                        true_count +=1
                except:
                    if cl_resis[cl_num+i-5][0]==me_resis[me_num+i-5][0]:
                            true_count +=1
            if true_count==5:
                mapping_dict[cl_resi[0]+cl_resi[1]] = me_resi[0] + me_resi[1]
                break
    return mapping_dict


def atoms_in_resi(Me_d: Dict[str,str],Me_PDB_lines: List[List[str]]) -> None:
    """
    This function checks to see how much of a given residue is in the QM
    region in the reference PDB to take the same amount of that residue
    later from the current PDB.  

    Will return f for full residue, c for only the carbonyl or xc for
    everything except the carbonyl or n for none of the residue.

    Parameters ---------- Me_D: Dict[str,str]
        Dictionary containing the atoms from the reference PDB and
        which region they belong to
    Me_PDB_Lines: List(List(str))
        List containing the lines of the reference PDB where each line
        is another list containing each category from the PDB

    Returns
    -------
    return_dict: Dict[str, str]
        the key contains the residue in the reference PDB and the value corresponds to how much of that residue is in the QM region
    """
    QM_list = []
    for key in Me_d:
        if 'Q' in key:
            for entry in Me_d[key]:
                QM_list.append(entry)
    resi_dict = {}
    return_dict = {}
    for i in range(len(Me_PDB_lines)):
        if Me_PDB_lines[i][3][:-1]+Me_PDB_lines[i][5] not in resi_dict:
            resi_dict[Me_PDB_lines[i][3][:-1]+Me_PDB_lines[i][5]] = [Me_PDB_lines[i][1]]
        else:
            resi_dict[Me_PDB_lines[i][3][:-1]+Me_PDB_lines[i][5]].append(Me_PDB_lines[i][1])
    for resi in resi_dict:
        atom_present_counter=0
        for atom in resi_dict[resi]:
            if int(atom) in QM_list:
                atom_present_counter+=1
        if atom_present_counter == len(resi_dict[resi]):
            return_dict[resi] = 'f'
        elif atom_present_counter ==2:
            return_dict[resi] = 'c'
        elif atom_present_counter == len(resi_dict[resi])-2:
            return_dict[resi] = 'xc'
        elif atom_present_counter ==0:
            return_dict[resi] = 'n'
        else:
            print(f'ERROR: Residue {resi} in Me_PDB does not match the three cases')
            sys.exit()
    return return_dict



def convert_dictionary(cutoff: str,template_path:str) -> None:
    """ 
    This function will take the dictionary for the template and find
    the corresponding residues in the current pdb and put the same amount
    of each residue in the QM region.  If there are residues, waters,
    or anything else that is in the current PDB but not in the template
    PDB, whether it is QM or MM will be determined by the cutoff radius.

    Parameters
    ----------
    cutoff: str
        cutoff in angstroms to determine the radius of the QM region for residues 
        that are in the current PDB but not in the reference PDB and for waters
    template_path: str
         path to the reference PDB

    Returns
    -------
    """

    #reading in the reference dict, reference pdb, and current pdb and parsing it into the correct data structure
    #path to reference pdb
    if os.path.isabs(template_path):
        Me_PDB_PATH = templte_path
    else:
        Me_PDB_PATH = '../' + template_path 
    Me_DICT_PATH = glob(f'{os.path.dirname(Me_PDB_PATH)}/dictionary.dat')[0]
    with open(Me_DICT_PATH, 'r') as dictfile:
        Me_d = json.load(dictfile)
    #path for the current pdb 
    Cl_PDB_PATH = 'cx_autocap_fixed.pdb'
    Me_PDB_lines = []
    with open(Me_PDB_PATH, 'r') as Me_PDB_file:
        all_Me_PDB_lines = Me_PDB_file.readlines()
        for line in all_Me_PDB_lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                Me_PDB_lines.append([line[0:6].strip(),line[6:11].strip(), line[11:16].strip(), line[16:20].strip(), line[20:22].strip(), line[22:26].strip(), line[26:38].strip(), line[38:46].strip(), line[46:54].strip(), line[54:60].strip(), line[60:66].strip(), line[66:79].strip()])
    Cl_PDB_lines = []
    with open(Cl_PDB_PATH, 'r') as Cl_PDB_file:
        all_Cl_PDB_lines = Cl_PDB_file.readlines()
        for line in all_Cl_PDB_lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                Cl_PDB_lines.append([line[0:6].strip(),line[6:11].strip(), line[11:16].strip(), line[16:20].strip(), line[20:22].strip(), line[22:26].strip(), line[26:38].strip(), line[38:46].strip(), line[46:54].strip(), line[54:60].strip(), line[60:66].strip(), line[66:79].strip()])
    
    #finds how much of the reference residues are in the QM region and maps each residue between the two PDBS
    grab_dict = atoms_in_resi(Me_d,Me_PDB_lines)
    mapping = match_resi_neighborhood(Me_PDB_lines, Cl_PDB_lines)
   
    #directly maps the link regions
    me_atom_dict = {}
    for line in Me_PDB_lines:
        me_atom_dict[line[1]]= line[3][:-1]+line[5]
    Cl_d = {'QM':[],'MM':[]}
    boundary_atoms = {}
    unique_keys = [key for key in Me_d.keys() if key!='QM' and key!='MM']
    for key in unique_keys:
        for atom in Me_d[key]:
            res = me_atom_dict[str(atom)]
            identifier = res+Me_PDB_lines[int(atom)-1][2]
            boundary_atoms[identifier] = key
    old_res = None
    with open('ligand.pdb') as lig:
        lig_lines = lig.readlines()
    lig_res = lig_lines[0][16:20].strip()
    for line in Cl_PDB_lines:
        if line[3]==lig_res:
            pass
        else:
        #look up Me resi
            if line[3][:-1]+line[5] in mapping and 'HOH' not in line[3]+line[5]:
                mapped_resi = mapping[line[3][:-1]+line[5]] 
                grab = grab_dict[mapped_resi]
                identifier = mapped_resi+line[2]
                if identifier in boundary_atoms:
                    key = boundary_atoms[identifier]
                    if key in Cl_d:
                        Cl_d[key].append(int(line[1]))
                    else:
                        Cl_d[key] = [int(line[1])]
                elif grab == 'f':
                    Cl_d['QM'].append(int(line[1]))
                elif grab == 'n':
                    Cl_d['MM'].append(int(line[1]))
                elif grab == 'c':
                    if line[2]=='C' or line[2]=='O':
                        Cl_d['QM'].append(int(line[1]))
                    else:
                        Cl_d['MM'].append(int(line[1]))
                elif grab == 'xc':
                    if line[2]=='C' or line[2]=='O':
                        Cl_d['MM'].append(int(line[1]))
                    else:
                        Cl_d['QM'].append(int(line[1]))
                else:
                    pass
            else:
                #if the current residue is not in the reference PDB, go based on the cutoff
                if line[3]+line[5]!=old_res:
                    old_res = line[3]+line[5]
                    cmd.reinitialize()
                    # Load PDB
                    cmd.load(Cl_PDB_PATH,"pdb")
                    cmd.show("sticks", "all")
                    cmd.label("all", "name")
                    cmd.select('close', f'organic and not solvent and not resname NME and not resname NMA and not resname ACE around {cutoff}')
                    cmd.select('QM',f'close and id {line[1]}')
                    inqm = cmd.count_atoms('QM')
                    #check distance
                    if inqm == 1:
                        Cl_d['QM'].append(int(line[1]))
                        QM = True
                    else:
                        Cl_d['MM'].append(int(line[1]))
                        QM = False
                else:
                    if QM:
                        Cl_d['QM'].append(int(line[1]))

                    else:
                        Cl_d['MM'].append(int(line[1]))
    with open('dictionary.dat', 'w+') as wfile:
        json.dump(Cl_d, wfile)
