from Bio.PDB import *
import os, os.path

import __main__
__main__.pymol_argv = ['pymol', '-qc']
import pymol
from pymol import cmd, stored
pymol.finish_launching()


in_dir = "./pdb/"
out_dir = "~/Desktop/Test code/Research/pdb_protein/"
parser = MMCIFParser()

rna_aa = ['A', 'G', 'U', 'T', 'C']
rna_chain = []
def main():

    counter = 5
    for file in os.listdir(in_dir):
        if counter == 0: break
        structure = parser.get_structure(file[:4], in_dir+file)
        for chain in structure[0]:
            #print(chain.get_full_id())
            for residue in chain:
                #print(residue.get_full_id())
                if residue.get_resname() in rna_aa:
                    #print(residue.get_resname())
                    rna_chain.append(chain)
                    break
        
        cmd.load(in_dir + file)
        for chain in rna_chain:
            cmd.remove('chain ' + chain.get_full_id()[2] )
            cmd.save(out_dir + file[:4]+'_protein.cif')
        cmd.reinitialize()
        counter -= 1

if __name__ == '__main__':
    main()
