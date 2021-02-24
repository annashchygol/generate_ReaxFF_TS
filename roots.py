import numpy
from pyomo.environ               import *
from pyomo.opt                   import SolverStatus, TerminationCondition, SolverFactory


def gcd(a,b):
    """Compute the greatest common divisor of a and b"""
    while b > 0:
        a, b = b, a % b
    return a

def lcm(a, b):
    """Compute the lowest common multiple of a and b"""
    return a * b / gcd(a, b)

def lcm_lst(l):
    res = l[0]
    for n in l[1:]:
        res = lcm(res, n)
        #print("RES:", res)
    return res

#parameters
def init_A ( model, i, j ):
    return A[i-1,j-1]

# A*x = 0
def equality_rule ( model, i ):
    return sum([ model.A[i,j] * model.x[j] for j in model.j ]) == 0 

#make sure that there is at least 1 x that is nonzero -- this prevents the trivial x=0 solution
def gte1_rule ( model ):
    return sum([ model.x[j] for j in model.j ]) >= 1 


def find_reac_coeff(A):
    #initialize A with given vector: 
    #uio_dict     = {'C': 8, 'O': 24, 'H': 8, 'Zr': 6}
    #h2o_dict     = {'C': 0, 'O': 1,  'H': 2, 'Zr': 0}
    #zro2_dict    = {'C': 0, 'O': 2,  'H': 0, 'Zr': 1} # -> negative
    #hformic_dict = {'C': 1, 'O': 2,  'H': 2, 'Zr': 0} # -> negative

    #A = numpy.matrix([[8,24,8,6], [0,1,2,0], [0,-2,0,-1], [-1,-2,-2,0]])
    #A = numpy.matrix([[8,0,0,-1], [24,1,-2,-2], [8,2,0,-2], [6,0,-1,0]])



    #initialize length of vector x
    m = A.shape[0]
    n = A.shape[1]

    model = ConcreteModel()

    #parameters
    def init_A ( model, i, j ):
        return A[i-1,j-1]

    # A*x = 0
    def equality_rule ( model, i ):
        return sum([ model.A[i,j] * model.x[j] for j in model.j ]) == 0 

    #sets 
    model.i      = Set( initialize = [x+1 for x in range(m)] )
    model.j      = Set( initialize = [x+1 for x in range(n)] )

    model.A      = Param( model.i, model.j , initialize = init_A )

    #variables
    # change NonNegativeIntegers to PositiveIntegers if all x must be > 0
    model.x      = Var( model.j , domain = NonNegativeIntegers)
    #model.x      = Var( model.j , domain = PositiveIntegers)

    #objective function 
    # -- option 1 : no objective function
    #model.obj = Objective( expr =  0, sense = minimize )

    #-- option 2 : minimize the sum of the x's
    model.obj = Objective( expr =  sum([ model.x[j] for j in model.j ]), sense = minimize )

    # A*x = 0
    model.con1 = Constraint( model.i, rule=equality_rule ) 

    #make sure that there is at least 1 x that is nonzero -- this prevents the trivial x=0 solution
    model.con2 = Constraint( rule=gte1_rule ) 

    #solve model
    opt = SolverFactory("couenne") #"glpk" #"couenne"
    results = opt.solve( model, tee=True, keepfiles=False )

    if results.solver.termination_condition == TerminationCondition.optimal:
        #print( "\nCONVERGED!!")
        #TODO: uncomment
        #model.obj.display()
        #model.x.display()
        #print("Variables -> ", model.x.values())
        x_val = [ model.x[j].value for j in model.j ]
        x_int = [int(round(x,0)) for x in x_val ]
        #print("VAL:", x_val, type(x_val))
        #print("INT:", x_int, type(x_val))
        coeff = lcm_lst(x_int)
        #print("COEFF:", coeff)
        coeff_lst = [int(coeff/x) for x in x_int] 
        #print("COEFF_LST:", coeff_lst)
        return coeff_lst
    else:
        #print ("\nINFEASIBLE!")
        return coeff_lst

################################
#initialize A with given vector: 
#uio_dict     = {'C': 8, 'O': 24, 'H': 8, 'Zr': 6}
#h2o_dict     = {'C': 0, 'O': 1,  'H': 2, 'Zr': 0}
#zro2_dict    = {'C': 0, 'O': 2,  'H': 0, 'Zr': 1} # -> negative
#hformic_dict = {'C': 1, 'O': 2,  'H': 2, 'Zr': 0} # -> negative
#A = numpy.matrix([[8,0,0,-1], [24,1,-2,-2], [8,2,0,-2], [6,0,-1,0]])
#A = numpy.matrix([[8,0,0,-1], [24,1,-8,-2], [8,2,0,-2], [6,0,-4,0]])
#find_reac_coeff(A) # <- Used to work! 

