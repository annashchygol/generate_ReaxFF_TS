
from itertools import combinations
from scm.plams import load, Molecule

def bond_type(bond):
    if bond.atom1.atnum < bond.atom2.atnum:
        return bond.atom1.symbol, bond.atom2.symbol 
    else:
        return bond.atom2.symbol, bond.atom1.symbol

def sort_atom(atom1, atom2):
    if atom1.atnum < atom2.atnum:
        return atom1, atom2 
    else:
        return atom2, atom1


def occurrences(molecule):
    bonds = [bond_type(b) for b in molecule.bonds]  
    types = set(bonds)
    return {key: bonds.count(key) for key in types}


def get_bonds_occurances(mol):
    symbols = [(i.atom1.symbol, i.atom2.symbol) for i in mol.bonds]
    print("SYMBOLS:", set(symbols))
    occurance_dict = {}
    for item in set(symbols):
        occ = symbols.count(item)
        occurance_dict[item] = occ
        #print(item, ":", occ)
    print("OCCURANCES:", occurance_dict)
    return occurance_dict

def neighbours_sorted_by_dist(mol, centre):
    neighbours = mol.neighbors(centre)
    neighbours.sort(key=lambda x: centre.distance_to(x), reverse=True)
    return neighbours


#@add_to_class(Molecule)
def guess_bonds_alt(mol, radius_coeff = 1.3, angle_threshold = 60, dist_threshold = 0.8):
    print (" >>> ALTERNATIVE BONDS !!!")
    mol.delete_all_bonds()
    print(mol)

    # Step 1: 
    for atom_1, atom_2 in combinations(mol, 2):
        if atom_1.distance_to(atom_2) < (atom_1.radius + atom_2.radius) * radius_coeff:
            mol.add_bond(atom_1, atom_2)

    # Step 2:
    for centre in mol:

        neighbours = mol.neighbors(centre)
        neighbours_sort_by_dist = neighbours_sorted_by_dist(mol, centre)


        for n1, n2 in combinations(neighbours_sort_by_dist, 2):
            #print ("n1, n2 = ", n1, n2)
            angle = centre.angle(n1,n2, result_unit='degree')
            #print ("ANGLE: ", angle)

            if angle < angle_threshold: #result_unit='degree'
                print (" TO remove ?")

                dist_1 = centre.distance_to(n1)
                dist_2 = centre.distance_to(n2)

                if dist_1 < dist_2 * dist_threshold:
                    mol.delete_bond(centre, n2)

                elif dist_2 < dist_1 * dist_threshold:
                    mol.delete_bond(centre, n1)
    return  mol


#def guess_bonds_new(mol, coeff=None, ):
def guess_bonds_new(mol, ij_radii_factor=1.3, ijk_angle_treshold=60, ijk_dist_treshold=0.8):
    symb = [a.symbol for a in mol.atoms]
    print("Atoms:", symb)
    print( "First atom", mol.atoms[1]) 

    mol.delete_all_bonds()

    # loop over pairs
    for i in range(1, len(mol.atoms)+1):
        atom_i = mol[i]
        for j in range(1, i):
            #print(i, j)
            atom_j = mol[j]
            distance_ij = atom_i.distance_to(atom_j)
            radiisum_ij = atom_i.radius + atom_j.radius
            if distance_ij < radiisum_ij * ij_radii_factor:
                # TODO: swap atom_i & atom_j if needed
                [atom_I, atom_J] = sort_atom(atom_i, atom_j)
                mol.add_bond(atom_I, atom_J)

    for atom_i in mol.atoms:
        print (mol.atoms.index(atom_i), '<<<<<<<<<')

        for atom_j in mol.neighbors(atom_i):
            for atom_k in mol.neighbors(atom_i):
                #print((mol.atoms.index(atom_j) , mol.atoms.index(atom_k)))
                if (mol.atoms.index(atom_j) <= mol.atoms.index(atom_k)):
                    break
                print((mol.atoms.index(atom_j) , mol.atoms.index(atom_k)))


                if (atom_i.angle(atom_j, atom_k, result_unit='degree') < ijk_angle_treshold):
                    print("Candidate to remove!", atom_i.angle(atom_j, atom_k, result_unit='degree'))
                    dist_j = atom_i.distance_to(atom_j)
                    dist_k = atom_i.distance_to(atom_k)

                    print ("Distances:", dist_j, dist_k)
                    if dist_j < dist_k * ijk_dist_treshold: # shorter one 
                        mol.delete_bond(atom_i, atom_j) #atom_i.delete_bond(atom_j)
                        print ("DELETING bond J !")
                    elif dist_k < dist_j * ijk_dist_treshold:
                        mol.delete_bond(atom_i, atom_k)
                        print ("DELETING bond K !")

    return mol
