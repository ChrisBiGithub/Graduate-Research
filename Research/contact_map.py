from Bio.PDB import *
import numpy as np
import warnings
warnings.filterwarnings('ignore')

AMINO_ACIDS = np.array(['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG'
               , 'SER','THR', 'VAL', 'TRP', 'TYR'])

"""Get contact map of the structure, will be use as adjacency matrix"""
def getContactMap(p, path):
    parser = PDBParser()
    structure = parser.get_structure(p, path)[0]
 
    p_chains = []
    for chain in structure:
        for residue in chain:
            if residue.get_resname() in AMINO_ACIDS:
                p_chains.append(chain)
                break
    
    
    c_map = []
    for chain in p_chains:
        for residue in chain:
            if residue.get_resname() not in AMINO_ACIDS: 
                continue
            temp = []
            
            for c in p_chains:
                for r in c:
                    if r.get_resname() not in AMINO_ACIDS: 
                        continue
                    print(residue.get_resname(), r.get_resname())
                    dis = np.linalg.norm(residue['CA'].get_coord() - r['CA'].get_coord())
                    if dis <= 7:
                        temp.append(1)
                    else:
                        temp.append(0)

            c_map.append(temp)
    
    
    print("Finished calculation")
    return c_map

def main():
    f = open('2XZN_contact.txt', 'w+')
    p = '2XZN'
    path = '../datasets/Tongji/benchmark/pdb/2XZN.pdb'
    
    contact_map = getContactMap(p, path)
    for i in contact_map:
        f.write('\t'.join(str(a) for a in i))
        f.write('\n')
    
    f.close()
    return


if __name__ == '__main__':
    main()