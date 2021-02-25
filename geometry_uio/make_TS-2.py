# 'The code is supplemntary for the work 'Studies of mechanical properties of MOFs with Reactive Force fields 
# using automatic parameter optimization and training set generation.'
# Any further work that uses the setting must therefore include the reference to the work.
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License.
# __author__ = "Anna Shchygol, shchygol@scm.com, Ganna.Shchygol@UGent.be"
# __license__ = "GPL"



import os
from scm.plams import *

init()

import numpy as np
import h5py as h5
from itertools import combinations, product
from bonds_guessing import guess_bonds_alt

def tors (atomI, atomJ, atomK, atomL):
    xyzI = np.array(atomI.coords)
    xyzJ = np.array(atomJ.coords)
    xyzK = np.array(atomK.coords)
    xyzL = np.array(atomL.coords)
    # Bond vectors
    f = xyzI - xyzJ
    g = xyzJ - xyzK
    h = xyzL - xyzK

    # Normal to i-j-k plane
    a  = np.cross(f,g)
    aa = np.dot(a,a)

    # Normal to j-k-l plane
    b  = np.cross(h,g)
    bb = np.dot(b,b)

    if aa > 1.0e-10 and bb > 1.0e-10:
        cosphi = np.dot(a,b)/np.sqrt(aa*bb)
        if cosphi > 1.0:
            cosphi = 1.0
        if cosphi < -1.0:
            cosphi = -1.0
        tIJKL = np.arccos(cosphi)
    else:
        tIJKL = 0.0

    rchk = np.cross(a,b)
    dum = np.dot(rchk,g)
    if dum > 0.0:
        tIJKL = - tIJKL
    return Units.convert(tIJKL, 'radian', 'degree') # radiands -> degrees 


def makeUpdateFile(path, fileName, inputStr = None):
    filePath = os.path.join(path, fileName)
    logging.info("FilePath " + filePath)

    if os.path.exists(filePath):
        os.remove(filePath)

    file = open(fileName, 'a')
    if inputStr:
        file.write(inputStr)
    file.close()
    return file 


@add_to_class(ADFResults)
def get_hirshfeld_charges(self):
    mapping = self._int2inp()
    charges_int = self.readkf('Properties', 'FragmentCharge Hirshfeld')
    charges_inp = [charges_int[mapping[i]-1] for i in range(len(charges_int))]
    return charges_inp


def get_charges_MBIS(mbis_path): 
    charges_str = ""
    w_charge = 0.05
    mol_id = mbis_path.split('/')[-2]

    with h5.File(mbis_path) as f:
        charges = f["charges"][:]
    #print(" >>> MBIS CHARGES:", charges)

    for i in range(0, len(charges)):
        charges_str += '{} {} {} {:.6}\n'.format(mol_id, w_charge, i+1, charges[i]) 
        #print("MBIS > ", mol_id, w_charge, i+1, charges[i])

    #print(" <<< MBIS CHARGES STR: >>> \n", charges_str, "<<< >>>") 

    return charges_str


def geo_data(abs_path):
    
    mol_id = abs_path.split('/')[-1]
    print("MOL ID:", mol_id)

    #asd = load(dill_path)
    s = Settings()
    s.input.basis.type='None'

    asd=ADFJob.load_external(path=abs_path, settings=s)
    res = asd.results

    out_mol = res.get_main_molecule() # -> Geometries -> Process it; 
    energy = res.readkf('Energy', 'Bond Energy')
    #charges = res.get_hirshfeld_charges() #res.readkf('Properties', 'FragmentCharge Hirshfeld') # 'AtomCharge Mulliken' Hirshfeld

    print("Out mol:", out_mol)
    print("Energy:", energy) # [U]: Hartree -> kcal / mol 
    #print("Hirshfield Charge:", charges) # TODO: MBIS charges instead, or separately;

    bonds = out_mol.guess_bonds() #unreliable
    new_mol = guess_bonds_alt(out_mol)

    print("Bonds:", out_mol.bonds)
    out_mol = new_mol # updating the molecule, with the other bond_guessing algorithm.
    bonds = new_mol.bonds
    print("Bonds:", out_mol.bonds)
    #for i in out_mol.bonds: print(i)
    #for i in out_mol.neighbors(out_mol[1]): print(i)

    mol = out_mol
    mol.set_atoms_id()

    # print("Neighbors:")
    # for i in mol.neighbors(mol[1]): print(i)  # <- list of atoms, connected to one another; 

    #out_str = ""
    charges_str = ""
    bonds_str = ""
    angles_str = ""
    tors_str = ""

    # WEIGHTS"
    w_charge = 0.05
    w_bonds = 0.01 * 2 
    w_angles = 5.00
    w_tors = 10.0
    w_energy = 1.0 / 4
     
#    TODO: comment back this block! Testing Hirshfeld charges; 
#    print("CHARGES:")
#    charges = get_charges_MBIS()
#    for i in range(0,len(charges)):
#        #print(mol_id, w_charge, i+1, format(charges[i], '.6f'))
#        charges_str += '{} {} {} {:.6}\n'.format(mol_id, w_charge, i+1, charges[i])

    #BONDS
    print("BONDS:")
    for bond in mol.bonds:
        #print(mol_id, w_bonds, bond.atom1.id, bond.atom2.id, format(bond.length(), '.6f'))
        bonds_str += '{} {} {} {} {:.6}\n'.format(mol_id, w_bonds, bond.atom1.id, bond.atom2.id, bond.length())
    #print(">> MY OUT_STR:\n%s" %(out_str))

    # ANGLES
    print("ANGLES:")
    for atom in mol: 
        #print("\n>>", atom.id, atom) #,  "\nCorresponding Angles:")
        for at1, at2 in combinations(mol.neighbors(atom), 2):
            angle = atom.angle(at1, at2)
            angle = Units.convert(angle, 'radian', 'degree') # <- Fixed conversion to degrees.
            #print(at1.id, atom.id, at2.id, "AT1:", at1, "AT2:", at2, "Angle:", angle)
            #print(mol_id, w_angles, at1.id, atom.id, at2.id, format(angle, '.6f')) # >> OUT
            angles_str += '{} {} {} {} {} {:.6}\n'.format(mol_id, w_angles, at1.id, atom.id, at2.id, angle)
    #print(">> MY OUT_STR:\n%s" %(out_str))

    # TORSIONS;
    print("TORSIONS:")
    for bond in mol.bonds:
        nl1 = mol.neighbors(bond.atom1)
        nl2 = mol.neighbors(bond.atom2)
        nl1.remove(bond.atom2)
        nl2.remove(bond.atom1)
        if (len(nl1) > 0 and len(nl2) > 0):
           for k, l in  product(nl1, nl2):

                # Extra condition for MOFs: to exclude torsions formed by 2 pairs of angles(3 atoms) if both > 90 degrees:
                angle1 = bond.atom1.angle(k, bond.atom2)
                angle2 = bond.atom2.angle(bond.atom1, l)
                angle1 = Units.convert(angle1, 'radian', 'degree')
                angle2 = Units.convert(angle2, 'radian', 'degree')
                if angle1 <= 150 and angle2 <= 150:

                    t = tors(k, bond.atom1, bond.atom2, l)         
                    if not np.isnan(t):
                        #print(mol_id, w_tors, k.id, bond.atom1.id, bond.atom2.id, l.id, format(t, '.6f'))
                        tors_str += '{} {} {} {} {} {} {:.6}\n'.format(mol_id, w_tors, k.id, bond.atom1.id, bond.atom2.id, l.id, t)

    #print(">> MY OUT_STR:\n%s" %(out_str))
    return (charges_str, bonds_str, angles_str, tors_str) # Now charges are empty


d = os.getcwd()
str_dir = os.path.join(d, 'structures-MBIS+smearq') # 'structures-MBIS+smearq'
str_list = os.listdir(str_dir)
print("LIST of STRUCTURES:", str_list)
abs_str_list = [os.path.join(str_dir, i) for i in str_list]
print("ABS LIST:", abs_str_list)

#Add path to the mbis.h5 file;
mbis_path_lst = []
for p in abs_str_list:
    mbis_path = os.path.join(p, 'mbis.h5')
    if os.path.exists(mbis_path):
        mbis_path_lst.append(mbis_path)
print("MBIS LIST:", mbis_path_lst)


# MBIS charges:
charges = [get_charges_MBIS(mbis_path) for mbis_path in mbis_path_lst]
out_str = ""
out_str += "CHARGE\n"
print("MBIS charges I've get >>>>>>>>", charges)
out_str += ''.join([i for i in charges]) # Replace it by MBIS charges:
out_str += "ENDCHARGE\n"


geo = [geo_data(path) for path in abs_str_list] # dill_lst


out_str += "GEOMETRY\n"
out_str += ''.join([i[1] for i in geo])
out_str += ''.join([i[2] for i in geo])
out_str += ''.join([i[3] for i in geo])
out_str += "ENDGEOMETRY\n"

#######################################
#Block with reactions:


# THIS is done in a separate script: reactions_4_TS.py        
#1. Return corresponding molecules, by ID 
#2. Get their chemical formula: get_formula; 
#3. Make vectors out of it
#4. Construct A matrix out of it
#5. Call roots.py -> Get needed coefficients ; 


#######################################

with open("trainset.in", "w") as text_file:
    text_file.write(out_str)

