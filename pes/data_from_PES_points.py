
# 1. Loop through the folders, containing *.rkf names 
#  --> Something similar to TS-2, since we have to extract all the files from #      there. 
# 2. From each PES#.rkf point extract: geo, gradients, Energy. 
# 3. Generate a unique identifier.  -> geo & trainset.in 

import argparse
import os

import numpy as np
import collections
import matplotlib.pyplot as plt

from scm.plams.tools.kftools import KFReader
from scm.plams.tools.units import Units
from scm.plams.core.errors import FileError
from scm.plams import Molecule, Atom, Bond

from build_TS.ts_from_mol import geo_data # Check!
from bonds_guessing import guess_bonds_alt

import pandas as pd

import kfhistory


ap = argparse.ArgumentParser()
ap.add_argument("resultdir", type=str)
args = ap.parse_args()

cur_dir = os.getcwd()
dirs = [x[0] for x in os.walk(cur_dir)]

PES_points = [(fname, dirname) for dirname in dirs
                          for fname in os.listdir(dirname)
                          if fname.endswith('.rkf')]
data = []

for point in PES_points:

    path = point[1].split('/')
    PESnum = str(point[0]).split('.')[0]
    direction = path[-1]

    trainset_id = 'scan_' + path[-2] + '_' + path[-1] + '_' + PESnum
    
    pes_atoms = path[-2].split('-')
    if len(pes_atoms) > 2:  # <------[skipping dihedral data] TODO: change this!
        continue

    pes_atom1_id = int(pes_atoms[0]) 
    pes_atom2_id = int(pes_atoms[1])

    kf = KFReader(os.path.join(point[1], point[0]))
    NumAtoms         = kf.read('Molecule','nAtoms') # Same for all the points;
    AtomNamesList    = kf.read('Molecule','AtomicNumbers')
    AtomSymbolList   = kf.read('Molecule','AtomSymbols')
    Coords           = Units.convert(kf.read('Molecule','Coords'), 'Bohr', 'Angstrom')
    energy           = Units.convert(kf.read('AMSResults', 'Energy'), 'au', 'kcal/mol') # 1 Hartree = 1 a.u
    
    npCoords = np.array(Coords)
    npCoords = npCoords.reshape(int(npCoords.size/3),3)
    
    mol = Molecule()
    AtomSymbolList = AtomSymbolList.strip()
    SymbolList = list(AtomSymbolList.split()) 
    SymbolList = [s.strip() for s in SymbolList]
    
    for at, coord in zip(SymbolList, npCoords):
        mol.add_atom(Atom(symbol=at, coords=coord))

    # Keep in mind, that coordinate indexing starts from 0;
    pes_at1 = Atom(symbol = SymbolList[pes_atom1_id - 1], coords = npCoords[pes_atom1_id - 1])
    pes_at2 = Atom(symbol = SymbolList[pes_atom2_id - 1], coords = npCoords[pes_atom2_id - 1])
    pes_bond = Bond(pes_at1, pes_at2)
    line = [(pes_atom1_id, pes_atom2_id), direction, PESnum, pes_bond.length(), energy]
    data.append(line)

column_names = ['id', 'str/sq', 'PES', 'Bond length [A]', 'Energy [kcal / mol]']
table_data = pd.DataFrame(data, columns = column_names)
table = table_data.set_index('id', 'PES')
table.sort_values(['id', 'str/sq', 'Bond length [A]']).to_csv("PES_bonds.csv") 

print("SUCCESS :)")
