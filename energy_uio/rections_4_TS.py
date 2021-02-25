# 'The code is supplemntary for the work 'Studies of mechanical properties of MOFs with Reactive Force fields 
# using automatic parameter optimization and training set generation.'
# Any further work that uses the setting must therefore include the reference to the work.
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License.
# __author__ = "Anna Shchygol, shchygol@scm.com, Ganna.Shchygol@UGent.be"
# __license__ = "GPL"


# Script to generate reactions for the training data, based on energies extracted from ADF calculation; 

# In: molecules, participating in the reactions; 
# Out: Energy block of the trainset.in; 
#1. Extract molecules from the list of molecules: []
#2. Loop over to get the reactions IDs: 2 loops;
#3. Function to get reactions coefficients [INT];
import os
from scm.plams import *

init()

import numpy as np
from roots import find_reac_coeff
from calc_sigma import get_sigma
from bonds_guessing import occurrences, guess_bonds_alt


@add_to_class(Molecule)
def get_formula_dict(self): # Will use it while extracting molecules using PLAMS; 
    ret = {}
    for atom in self:
        if atom.symbol not in ret:
            ret[atom.symbol] = 0
        ret[atom.symbol] += 1
    return ret


def get_reaction_matrix(substrates, products, atom_order):
    #A = numpy.matrix([[8,0,0,-1], [24,1,-2,-2], [8,2,0,-2], [6,0,-1,0]])
    tmp = []
    for item in substrates:
        tmp.append([item[atom] for atom in atom_order])
    for item in products:
        tmp.append([-1*item[atom] for atom in atom_order])
    return np.array(tmp).T


def mol_reac_info(path_to_job_dir, atom_order = ['C', 'O', 'H', 'Zr']):  # <-default atom order; could be done smarter; 

    print("This path:", path_to_job_dir)
    mol_id = path_to_job_dir.split('/')[-1]

    s = Settings()
    s.input.basis.type='None'

    asd=ADFJob.load_external(path=path_to_job_dir, settings=s)
    res = asd.results

    energy = res.readkf('Energy', 'Bond Energy') # Up to which order of magnitude ? 
    hartree_to_kcalpmol = 627.509474 # Conversion coefficient 
    energy = energy * hartree_to_kcalpmol
    mol = asd.results.get_main_molecule()
    mol_dict = mol.get_formula_dict()

    upd_mol = guess_bonds_alt(mol)

    d = {}
    d['id'] = mol_id
    d['energy']  = energy
    d['dict'] = add_keys_to_dict(mol_dict, atom_order)
    d['bonds_occ'] = occurrences(upd_mol)

    #print("D:", d)
    return d 


def add_keys_to_dict(mol_dict, atom_order):
    #print("Start dict:", mol_dict)
    for atom in atom_order:
        if atom not in mol_dict.keys():
            mol_dict[atom] = 0
            #print("This key ->", atom)

    #print("End dict:", mol_dict, '\n')
    return mol_dict


def get_reaction_str(substrate_mols, product_mols, coeff_lst): # N(coeffs in the list) =  N(substr_mol) + N(prod_mol)
    #0. TODO: calculate the sigma(weight), using calc_sigma.py 


    # These extra dictionaries are created to keep compatiability with the calc_sigma code (by Toon)
    bonds_counts = {}
    reaction = {}
    coeff_iter = iter(coeff_lst)
    
    for mol in substrate_mols:
        bonds_counts[mol['id']] = mol['bonds_occ']
        reaction[mol['id']] = next(coeff_iter)
        
    for mol in product_mols:
        bonds_counts[mol['id']] = mol['bonds_occ']
        reaction[mol['id']] = -next(coeff_iter)

    print(" >> BONDS count:", bonds_counts)
    print(" >> corresp keys", bonds_counts.keys())
    print(" >> COEFFICIENTS:", coeff_lst)
    print(" >> REACTION:", reaction)

    #get reaction: #TODO: SIMPLIFY!

    ##############################################    
    # The rest of the string:
    reaction_str = '' #'1.0 ' # starting default weight: 
    reaction_energy = 0

    substr_num = len(substrate_mols)

    for i in range(0, substr_num):
        cur_coeff = coeff_lst[i]
        cur_id = substrate_mols[i]['id']
        cur_energy = substrate_mols[i]['energy']
        reaction_str += ' + ' + '/'.join([cur_id,  str(cur_coeff)])
        reaction_energy += cur_energy / cur_coeff
        reaction[cur_id] = cur_coeff


    for i in range(0, len(product_mols)):
        cur_coeff = coeff_lst[i + substr_num]
        cur_id = product_mols[i]['id']
        cur_energy = product_mols[i]['energy']
        reaction_str += ' - ' + '/'.join([cur_id,  str(cur_coeff)])
        reaction_energy -= cur_energy / cur_coeff
        reaction[cur_id] = (-1) * cur_coeff


    weight = get_sigma(bonds_counts, reaction)

    reaction_str += '  ' + str(reaction_energy)
    print ("Reaction:  ", reaction_str)

    weighted_reaction_str = str(weight) + reaction_str
    print ("Weighted react:  ", reaction_str)
    
    return weighted_reaction_str #reaction_str 

########################################
# UIO + ZR block: 

#1. Extract all we need from the folder:
start_path = '/Users/anna/Documents/Training_sets/newZrMOF/structures-MBIS+smearq/'
struct_lst = os.listdir(start_path)

#2. Make 2 lists: UIO & ZR;  
uio_struct_lst = [s for s in struct_lst if s.startswith('uio')]
zr_struct_lst  = [s for s in struct_lst if s.startswith('zr')]

uio_path_lst = [os.path.join(start_path, i) for i in uio_struct_lst] 
zr_path_lst  = [os.path.join(start_path, i) for i in zr_struct_lst]

h2o_path = [os.path.join(start_path, s) for s in struct_lst if s.startswith('h2o')][0]
hformic_path = [os.path.join(start_path, s) for s in struct_lst if s.startswith('hformic')][0]

print("Struct lst:", struct_lst)
print("UIO lst:", uio_struct_lst)
print("ZR  lst:", zr_struct_lst)
print("UIO path lst:", uio_path_lst)
print("ZR path lst:", zr_path_lst)

uio_items = [mol_reac_info(p) for p in uio_path_lst] # TODO: should be a path to job! 
zr_items  = [mol_reac_info(p) for p in zr_path_lst]

h2o = mol_reac_info(h2o_path)
hformic = mol_reac_info(hformic_path)

atom_order = ['C', 'O', 'H', 'Zr']
h2o['dict']     = add_keys_to_dict(h2o['dict'], atom_order)
hformic['dict'] = add_keys_to_dict(hformic['dict'], atom_order)

reactions = [get_reaction_matrix([uio['dict'], h2o['dict']], [zr['dict'], hformic['dict']], atom_order)
             for uio in uio_items for zr in zr_items]

lst_of_coeffs_lst =  [find_reac_coeff(A) for A in reactions]

react_items = [(uio, zr) for uio in uio_items for zr in zr_items]

react_str_lst = [get_reaction_str([r[0], h2o], [r[1], hformic], c)
                 for r,c in zip(react_items, lst_of_coeffs_lst)]

print('####### Reactions #########')
print('LEN of reactions:', len(reactions))
print('LEN of coefficients:', len(lst_of_coeffs_lst))

########################################


out_str = ""
out_str += "ENERGY\n"
out_str += '\n'.join([i for i in react_str_lst])
out_str += "\nENDENERGY\n"


########################################

print (" >> FILE:\n",  out_str)
with open("trainset.in-2", "w") as text_file:
    text_file.write(out_str)

finish()
