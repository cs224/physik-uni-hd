
from sympy import *
from sympy.solvers import solve
# init_printing()

dotV = Symbol("\dot V", positive=True)
eta = Symbol("\eta", positive=True)
R = Symbol("R", positive=True)
p_B = Symbol("p_B", positive=True)
p_0 = Symbol("p_0", positive=True)
l = Symbol("l", positive=True)
e1 = Eq(dotV, (pi/8*eta)*R**4*((p_B-p_0)/l))

eR = Eq(R, solve(e1, R)[1])

# ----------------------

dotV_0 = Symbol("\dot V_0", positive=True)

dotV_2 = Symbol("\dot V_2", positive=True)
R_2 = Symbol("R_2", positive=True)
l_2 = Symbol("l_2", positive=True)

dotV_3 = Symbol("\dot V_3", positive=True)
R_3 = Symbol("R_3", positive=True)
l_3 = Symbol("l_3", positive=True)

e2 = Eq(dotV_2, dotV_0)
e3 = Eq(dotV_3, dotV_0)


eR_2_over_R_3 = simplify(eR.rhs.subs(dotV, dotV_0).subs(l, l_2) / eR.rhs.subs(dotV, dotV_0).subs(l, l_3))

# ----------------------

Re = Symbol("Re", positive=True)
Re_K = Symbol("Re_K", positive=True)
rho = Symbol("\\rho", positive=True)
barv = Symbol("\\bar v", positive=True)
R = Symbol("R", positive=True)
e4 = Eq(Re,2*R*rho*barv/eta)

e5 = Eq(dotV, pi*R**2*barv)

solve([e4, e5], [barv, R], dict=True)

e4.subs(barv, solve(e5, barv)[0])
