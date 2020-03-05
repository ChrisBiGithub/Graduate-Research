from Bio.PDB import *
import numpy as np
import warnings
import os
from PCM import PCM

warnings.filterwarnings('ignore')

AMINO_ACIDS = np.array(['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER','THR', 'VAL', 'TRP', 'TYR'])

def outputFile(cm, cm_no, protein_name, out_dir):
    f_name = protein_name + '_cm' + str(cm_no) + '.txt'
    f = open(out_dir + f_name, 'w+')
    for i in cm:
        f.write('\t'.join(str(a) for a in i))
        f.write('\n')

    f.close()
    print("Finished {}".format(f_name))
    return

def main():
    directory = '../datasets/Tongji/benchmark/pdb/'
    out_dir = './contact_maps/'

    pdb_list = []
    for file in os.listdir(os.getcwd() + '/trial_data'):
        if not file.endswith('.txt'):
            continue
        pdb_list.append(file[:4])

    tester = 2
    for file in os.listdir(directory):
        if tester == 0:
            break

        if not file.endswith('.pdb'):
            continue

        name = file[:4]
        if name in pdb_list:
            pcm = PCM(name, directory + file)
            cm1, _ = pcm.getContactMap1()
            outputFile(cm1, 1, name, out_dir)

            cm2, _ = pcm.getContactMap2()
            outputFile(cm2, 2, name, out_dir)

            cm3, _ = pcm.getContactMap3()
            outputFile(cm3, 3, name, out_dir)

            cm4, _ = pcm.getContactMap4()
            outputFile(cm4, 4, name, out_dir)

            tester -= 1

    return


if __name__ == '__main__':
    main()
