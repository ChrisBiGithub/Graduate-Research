import numpy as np
import scipy.sparse as sp
from Bio.PDB import *
from matplotlib import pyplot as plt
from mendeleev import element as ele

AMINO_ACIDS = np.array(['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG'
               , 'SER','THR', 'VAL', 'TRP', 'TYR'])

class PCM:
    def __init__(self, p, path):
        self.p = p
        self.path = path

    #Find center of mass for side chains
    @staticmethod
    def getSideComCoord(residue):
        if residue.get_resname() == 'GLY':
            return residue['CA'].get_coord()
        else:
            temp = np.array([0., 0., 0.])
            sc = 0

            sc_flag = False
            for atom in residue:
                if atom.get_id() != 'CB' and not sc_flag:
                    continue

                sc_flag = True
                temp += atom.get_coord()
                sc += 1

        return temp/sc


    """First version of contact map, unweight and undirected"""
    def getContactMap1(self):
        parser = PDBParser()
        structure = parser.get_structure(self.p, self.path)[0]

        res_list = []

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

                res_list.append(residue.get_resname())

                temp = []
                counter = 0
                for c in p_chains:
                    for r in c:
                        if r.get_resname() not in AMINO_ACIDS:
                            continue
                        counter += 1
                        dis = np.linalg.norm(residue['CA'].get_coord() - r['CA'].get_coord())
                        if dis <= 7:
                            temp.append(1)
                        else:
                            temp.append(0)
                c_map.append(temp)


        return c_map, res_list

    #Three different criterions, unweighted and undirected
    def getContactMap2(self):
        parser = PDBParser()
        structure = parser.get_structure(self.p, self.path)[0]

        res_list = []

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

                res_list.append(residue.get_resname())

                temp = []
                counter = 0
                for c in p_chains:
                    for r in c:
                        if r.get_resname() not in AMINO_ACIDS:
                            continue
                        counter += 1
                        #Carbon alpha distance
                        dis_alpha = np.linalg.norm(residue['CA'].get_coord() - r['CA'].get_coord())

                        #Carbon beta distance, if the residue is Glycine, carbon alpha is used
                        if residue.get_resname() == 'GLY' and r.get_resname() != 'GLY':
                            dis_beta = np.linalg.norm(residue['CA'].get_coord() - r['CB'].get_coord())
                        elif residue.get_resname() != 'GLY' and r.get_resname() == 'GLY':
                            dis_beta = np.linalg.norm(residue['CB'].get_coord() - r['CA'].get_coord())
                        elif residue.get_resname() == 'GLY' and r.get_resname() == 'GLY':
                            dis_beta = np.linalg.norm(residue['CA'].get_coord() - r['CA'].get_coord())
                        elif residue.get_resname() != 'GLY' and r.get_resname() != 'GLY':
                            dis_beta = np.linalg.norm(residue['CB'].get_coord() - r['CB'].get_coord())

                        #Side chain distance
                        #dis_sc = 10
                        dis_sc = np.linalg.norm(self.getSideComCoord(residue) - self.getSideComCoord(r))

                        if dis_alpha <= 7 or dis_beta <= 7 or dis_sc <= 7:
                            temp.append(1)
                        else:
                            temp.append(0)
                c_map.append(temp)

        return c_map, res_list

    #Three different criterions, weighted and undirected
    def getContactMap3(self):
        parser = PDBParser()
        structure = parser.get_structure(self.p, self.path)[0]

        res_list = []

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

                res_list.append(residue.get_resname())

                temp = []
                counter = 0
                for c in p_chains:
                    for r in c:
                        if r.get_resname() not in AMINO_ACIDS:
                            continue
                        counter += 1
                        #Carbon alpha distance
                        dis_alpha = np.linalg.norm(residue['CA'].get_coord() - r['CA'].get_coord())

                        #Carbon beta distance, if the residue is Glycine, carbon alpha is used
                        if residue.get_resname() == 'GLY' and r.get_resname() != 'GLY':
                            dis_beta = np.linalg.norm(residue['CA'].get_coord() - r['CB'].get_coord())
                        elif residue.get_resname() != 'GLY' and r.get_resname() == 'GLY':
                            dis_beta = np.linalg.norm(residue['CB'].get_coord() - r['CA'].get_coord())
                        elif residue.get_resname() == 'GLY' and r.get_resname() == 'GLY':
                            dis_beta = np.linalg.norm(residue['CA'].get_coord() - r['CA'].get_coord())
                        elif residue.get_resname() != 'GLY' and r.get_resname() != 'GLY':
                            dis_beta = np.linalg.norm(residue['CB'].get_coord() - r['CB'].get_coord())

                        #Side chain distance
                        #dis_sc = 10
                        dis_sc = np.linalg.norm(self.getSideComCoord(residue) - self.getSideComCoord(r))

                        """Set the weight to the closest link, could be changed later"""
                        if dis_alpha <= 7 or dis_beta <= 7 or dis_sc <= 7:
                            temp.append(7 - min(dis_alpha, dis_beta, dis_sc))
                        else:
                            temp.append(0)
                c_map.append(temp)

        return c_map, res_list

    #Three different criterions, weighted and undirected
    def getContactMap4(self):
        parser = PDBParser()
        structure = parser.get_structure(self.p, self.path)[0]

        res_list = []

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

                res_list.append(residue.get_resname())

                temp = []
                counter = 0
                for c in p_chains:
                    for r in c:
                        if r.get_resname() not in AMINO_ACIDS:
                            continue
                        counter += 1
                        #Carbon alpha distance
                        dis_alpha = np.linalg.norm(residue['CA'].get_coord() - r['CA'].get_coord())

                        #Carbon beta distance, if the residue is Glycine, carbon alpha is used
                        if residue.get_resname() == 'GLY' and r.get_resname() != 'GLY':
                            dis_beta = np.linalg.norm(residue['CA'].get_coord() - r['CB'].get_coord())
                        elif residue.get_resname() != 'GLY' and r.get_resname() == 'GLY':
                            dis_beta = np.linalg.norm(residue['CB'].get_coord() - r['CA'].get_coord())
                        elif residue.get_resname() == 'GLY' and r.get_resname() == 'GLY':
                            dis_beta = np.linalg.norm(residue['CA'].get_coord() - r['CA'].get_coord())
                        elif residue.get_resname() != 'GLY' and r.get_resname() != 'GLY':
                            dis_beta = np.linalg.norm(residue['CB'].get_coord() - r['CB'].get_coord())

                        #Side chain distance
                        #dis_sc = 10
                        dis_sc = np.linalg.norm(self.getSideComCoord(residue) - self.getSideComCoord(r))

                        """Set the weight to the closest link, could be changed later"""
                        if dis_alpha <= 7 or dis_beta <= 7 or dis_sc <= 7:
                            temp.append(7 - min(dis_alpha, dis_beta, dis_sc))
                        else:
                            temp.append(0)
                        """If the bond is polypeptide bond, double the score, this can be change"""
                        if chain == c and abs(residue.get_id()[1] - r.get_id()[1]) == 1:
#                             print(True)
                            temp[-1] = temp[-1] * 2

                c_map.append(temp)

        return c_map, res_list

