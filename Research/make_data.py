#!/usr/bin/env python
# coding: utf-8

# ### 

# ##### Quick notes on the precomputed .mat files:  it also contains label for RNA, have to remove those labels before using them

# In[1]:


import warnings
warnings.filterwarnings('ignore')


# In[2]:


import Bio
import scipy.io
import numpy as np
import os, os.path
import math
from Bio.PDB import *
from sklearn.preprocessing import OneHotEncoder
from mendeleev import element as ele


# In[3]:


AMINO_ACIDS = np.array(['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG'
               , 'SER','THR', 'VAL', 'TRP', 'TYR'])


# In[4]:


"""
Extract labels from tongji dataset
"""
lab_path = "../datasets/Tongji/benchmark/precomputed/"
out_file = open("./label/labels.txt", "w+")

odditties = ['1FJG', '1JJ2', '1N32']

for file in os.listdir(lab_path):
#     print(file)
    mat = scipy.io.loadmat(lab_path + file)
    name = mat['File'][0][0][2][0]
#     if name in odditties:
#         continue #skip the wierd ones for now
    labels = mat['File'][0][0][3]

    out_file.write((str(name) + '\t').rstrip('\n'))
    for i, label in enumerate(labels):
        out_file.write((str(label[0]) + '\t').rstrip('\n'))

    out_file.write('\n')

out_file.close()


# In[5]:


in_path = "../datasets/Tongji/benchmark/pdb/"
out_path = "./trial_data/"


# In[6]:


temp = []
for file in os.listdir(lab_path):
    temp.append(file[:4])

print(len(set(temp)))
print(set(temp))


# In[7]:


"""
Load protein structures based on the labels
"""
labels_file_path = './label/labels.txt'

proteins = []
missing_protiens = []
protein_structures = {}
parser = PDBParser()

labels_file = open(labels_file_path, 'r')

for label in labels_file:
    if label.split('\t')[0] in odditties:
        continue #exclude the one with weird pdbs
    proteins.append(label.split('\t')[0])

for p in proteins:
    file_name = p + '.pdb'
    if file_name not in os.listdir(in_path):
        missing_protiens.append(p)

if len(missing_protiens) == 0:
    print("There are no missing pdb")
else:
    print("Missing proteins pdb:", missing_protiens)

for p in proteins:
    protein_structures[p] = parser.get_structure(p, in_path + p + '.pdb')

print(protein_structures)
labels_file.close()


# In[8]:


print(set(temp) - set(proteins))


# In[52]:


"""
Clean up the ASA files to remove non-amino acid residue
"""
asa_path = './SASAs/proteins/'
for p in proteins:
    with open(asa_path + p + '_protein.pdb_SASA.txt', 'r') as f:
        lines = f.readlines()
    with open(asa_path + p + '_SASA.txt', 'w') as f:
        for line in lines:
            if line.rstrip('\n').split(' ')[2] in AMINO_ACIDS:
                f.write(line)
        


# In[10]:


"""
Calculating protein COM base on coordinates and mass (equation from CALCOM)
"""
def proteinCOM(p_chains):
    atoms = ['N', 'C', 'O', 'H', 'S']
    com = np.array([0., 0., 0.])
    total_mass = 0
    
    for chain in p_chains:
        for residue in chain:
            for atom in residue:
                name = atom.get_id()[0]
#                 print(atom.get_id())
                if name not in atoms:
#                     print(atom.get_id())
                    break
                a = ele(name)
                
                com += atom.get_coord() * a.mass
                total_mass += a.mass
    
    result = com/total_mass
    return result


# In[ ]:


"""
Make files with features
"""

'''Load encoder'''
encoder = OneHotEncoder(sparse = False)
encoder.fit(AMINO_ACIDS.reshape(-1, 1))

'''Get file path ready'''
write_path = './trial_data/'

'''Get labels'''
label_file = open('./label/labels.txt', 'r')
labels = {}
for line in label_file:
    content = line.rstrip('\n').split('\t')
    labels[content[0]] = content[1:-1]

i = 1
for p in proteins:
    if os.path.isfile(write_path + p + '.txt'):
        print("Trial data for protein {} aleady exists.".format(p))
        continue
        
    print("Processing protein {}.".format(p))
    
    structure = protein_structures[p]
    res_label = labels[p]

    asa_file = open(asa_path + p + '_SASA.txt', 'r')
    out_file = open(write_path + p + '.txt', 'w+')


    #Separate the proteins and RNA from each other in the complex
    p_chains = []
    n_chains = []
    for chain in structure[0]:
        for residue in chain:
            #No need to continue if the chain is not an amino acid chain
            if residue.get_resname() not in AMINO_ACIDS:
                if chain not in n_chains:
                    n_chains.append(chain.get_id)
                break
            if chain not in p_chains:
                p_chains.append(chain)
                break

    #Calculte center of mass of the protein structure
    #com = proteinCOM(p_chains)
    com = proteinCOM(p_chains)

    '''Temporarily loading all residues'''
    res_list = []
    for chain in p_chains:
        for residue in chain:
            if residue.get_resname() not in AMINO_ACIDS:
                continue
            res_list.append(residue)


    '''Getting labels for only amino acids'''
    ares_label = res_label[-len(res_list):]
    
    for residue, lab in zip(res_list, ares_label):
        #Get the one hot representation for the residue 
        if residue.get_resname() not in AMINO_ACIDS:
            continue
        
        output = {} #temporarily saving the output

        r = np.array([residue.get_resname()])

        one_hot = encoder.transform(r.reshape(-1, 1))[0]
        output['one_hot'] = one_hot

        '''Get index'''
        '''maybe try the sinusoidal '''
        index = residue.get_id()[1]
        output['index'] = str(index)

        '''Get distance from center of mass'''
        #Some residue's CA is not recorded in pdb files
        dis = 0 
        try:
            dis = np.linalg.norm(residue['CA'].get_coord() - com)
        except:
            print('No CA in residue {}.'.format(residue.get_resname()))
        output['dis'] = str(dis)

        '''Get hydrophobicity'''
        hydro = 0
        if residue.get_resname() == 'GLY':
            output['hydrophobicity'] = str(hydro)
        else:
            temp = np.array([0., 0., 0.])
            sc_mass = 0

            sc_flag = False
            for atom in residue:
                if atom.get_id() != 'CB' and not sc_flag:
                    continue

                sc_flag = True
                mass = ele(atom.get_id()[0]).mass
                temp += atom.get_coord() * mass
                sc_mass += mass


            sc_com = temp/sc_mass
            dot = np.dot(sc_com, com)
            norm_sc = np.linalg.norm(sc_com)
            norm_protein = np.linalg.norm(com)
            hydro = dot / (norm_sc * norm_protein)
            output['hydrophobicity'] = str(hydro)

############################################################################################################            
        '''Get ASA'''
        '''ASA is calculated through another python script and saved in certain location'''
        chain_id = residue.get_parent().get_id()
        res_id = str(residue.get_id()[1])
        res_name = residue.get_resname()
        
        for line in asa_file:
            content = line.rstrip('\n').split(' ')
            if chain_id == content[0] and res_id == content[1] and res_name == content[2]:
                output['ASA'] = content[-1]
                break

        '''Get propensity'''


        '''Get label'''
        output['label'] = lab
        
        out_file.write(residue.get_resname().rstrip('\n') + '\t')
        for key in output:
            if key == 'one_hot':
                out_file.write('\t'.join(str(t) for t in output[key]).rstrip('\n') + '\t')
            else:
                out_file.write(output[key] + '\t')
        out_file.write('\n')


    out_file.close()
    asa_file.close()
    print("Protein {} is done.".format(p))    


label_file.close()



# In[ ]:


"""Temporarily disabled code"""
# """Downloading data from RCSB PDB"""
# if(not os.path.isdir('./pdb')):
#     pdbl = PDBList()
#     download_proteins_set = set(proteins)

#     for i in download_proteins_set:
#         pdbl.retrieve_pdb_file(pdb_code = i,pdir='../datasets/Tongji/pdb')

###############################################################################################################

# """Information about the downloaded files"""

# print("Number of unique proteins the dataset {} is {}.".format("Tongji", len(set(proteins))))

# pdb_dir = '../datasets/Tongji/benchmark/pdb'
# print("Number of protein structure have been downloaded is %2d." %len([name for name in os.listdir(pdb_dir) if os.path.isfile(os.path.join(pdb_dir, name))]))

# """Files that are obselete (for whatever reasons)"""
# pdbs = [name[:4].upper() for name in os.listdir(pdb_dir)]
# obsolete_pdbs = list(set(proteins) - set(pdbs))
# print("Number of obsolete pdbs of RPI2241 datasets in RCSB protein databank is %2d, list as follow:" %(len(obsolete_pdbs)))
# i = 0 
# for name in obsolete_pdbs:
#     if i < 12:
#         print(name, end = '\t')
#     else:
#         i = 0
#         print()
#         print(name, end = '\t')
    
#     i += 1

# """Empty header?"""
# structure = protein_structures[0]
# print(structure.header)

# print(structure[0]['A'].get_residues())

###############################################################################################################

# """Get the protein chain"""
# protein_chains = []
# for chain in structure[0]:
#     for residue in chain:
#         if residue.get_resname() in AMINO_ACIDS and chain not in protein_chains:
#             protein_chains.append(chain)
#         break
#     continue

# for chain in protein_chains: 
#     print(chain)

###############################################################################################################

# """Make a dictionary that contain all the positive chain and residue sequence (only amino acid)"""

# file.seek(0) #move back to top of file
# positive_dic = {}

# for line in file:
#     content = line.split('\t')
    
#     if content[0][:4] in obsolete_pdbs:
#         continue
        
#     if content[0][:4] in positive_dic:
#         positive_dic[content[0][:4]].update({content[0][-1]: content[2]})
#     else:
#         positive_dic[content[0][:4]] = {content[0][-1]: content[2]}
        
# # print(positive_dic)
# print("Number of protein in dictionary is %2d." %len(positive_dic))


# In[ ]:




