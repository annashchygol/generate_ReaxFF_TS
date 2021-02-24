# Eaxmple reaction
#1.0  + uio66zr_C08/24 + h2o/6 - zr2o4A/8 - hformic/3  -25.218958938920622

# Get python 3 division, even in python 2.
from __future__ import division
# Get print function from python 3, even in python 2.
from __future__ import print_function

# You need to derive the following dictionary from the geometry
#bond_counts = {
#    # Make sure the keys are sorted, i.e. lowest atomic number first (or last).
#    "uio66zr_C08": {(8, 40): 40, (6, 8): 16, (1, 6): 8},
#    "h20": {(1, 8): 2},
#    "zr2o4A": {(8, 40): 6},
#    "hformic": {(6, 8): 2, (1, 8): 1, (1, 6): 1},
#}

bond_counts = {
    # Make sure the keys are sorted, i.e. lowest atomic number first (or last).
    "uio66zr_C08": {('H', 'C'): 8, ('C', 'O'): 16, ('O', 'Zr'): 40},
    "h20": {('H', 'O'): 2},
    "zr2o4A": {('O', 'Zr'): 6},
    "hformic": {('H', 'C'): 1, ('H', 'O'): 1, ('C', 'O'): 2},
}


# This comes from your stochiometry calculation
# TODO: Be careful with the floating point operations!
reaction = {
    "uio66zr_C08": 24, #  1/24
    "h20": 6,          #  1/6
    "zr2o4A": -8,      # -1/8
    "hformic": -3,     # -1/3
}

def get_sigma(bond_counts, reaction): 
    # Compute the net number of bonds formed (+) or broken (-)
    bonds_balance = {}
    for name, bonds in bond_counts.items():
        coeff =  1 / reaction[name]
        for element_pair, nbond in bonds.items():
            bonds_balance[element_pair] = \
                bonds_balance.get(element_pair, 0.0) \
                + coeff*nbond

    print(bonds_balance)
    sigma_per_bond = 1.0  # in kcal/mol
    # Absolute value is important! in the next line. Bond broken or formed
    # is in practice the same thing. This makes the algorithm insensitive to
    # reaction reversal.
    sigma = sigma_per_bond*sum(abs(nbond) for nbond in bonds_balance.values())
    print("sigma:", sigma/4 ) 
    return sigma / 4  # <- Feb 24. 2021: Decision to change sigma, based on D-optimality;


get_sigma(bond_counts, reaction)

# My output:
#{(1, 8): 0.0, (8, 40): 0.9166666666666665, (1, 6): 0.0,6,8):0.0}                                           
#0.916666666667  - 