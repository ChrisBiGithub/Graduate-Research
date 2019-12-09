import sys
import os, os.path

import __main__
__main__.pymol_argv = ['pymol', '-qc']
import pymol
from pymol import cmd, stored
pymol.finish_launching()

cmd.set('dot_solvent', 1)
cmd.set('dot_density', 3)

in_dir1 = './pdb/'
out_dir1 = './SASAs/complexes/'

in_dir2 = './pdb_protein/'
out_dir2 = './SASAs/proteins/'

def main():
    in_dir, out_dir = '', ''
    if len(sys.argv) != 2 and (sys.argv[1] != 'protein' or sys.argv[1] != 'complex'):
        print('Insufficient argument. (Can only be either protein or complex)')
        print(sys.argv[1])
        return

    if sys.argv[1] == 'protein':
        in_dir, out_dir = in_dir2, out_dir2
    else:
        in_dir, out_dir = in_dir1, out_dir1

    counter = 2
    for file in os.listdir(in_dir):
        if counter == 0:
            break
        if file[-4:] != '.cif':
            continue

        print(file)
        outfile = open(out_dir + file + '_SASA.txt', 'w+')
        cmd.load(in_dir + file)
        stored.residues = []
        cmd.iterate('name ca', 'stored.residues.append(resi)')

        sasa_per_residue = []
        for i in stored.residues:
            #sasa_per_residue.append(cmd.get_area('resi %s' % i))
            outfile.write('resi %s' % i + ' ' + str(cmd.get_area('resi %s' %i)))
            outfile.write('\n')
        outfile.close()
        print('Finished calculating protein {}'.format(file))
        cmd.reinitialize()
        counter -= 1

    return

if __name__ == '__main__':
    main()
