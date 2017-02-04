
from sympy import *
from sympy.solvers import solve
# init_printing()

# V˙2=π8ηR42pB−p0l2

dotV_2 = Symbol("\dot V_2", positive=True)
eta = Symbol("\eta", positive=True)
R_2 = Symbol("R_2", positive=True)
p_B = Symbol("p_B", positive=True)
p_0 = Symbol("p_0", positive=True)
l_2 = Symbol("l_2", positive=True)
e1 = Eq(dotV_2, (pi/8*eta)*R_2**4*((p_B-p_0)/l_2))
solve(e1, R_2)

dotV_3 = Symbol("\dot V_3", positive=True)
R_3 = Symbol("R_3", positive=True)
l_3 = Symbol("l_3", positive=True)
e2 = Eq(dotV_3, (pi/8*eta)*R_3**4*((p_B-p_0)/l_3))

e3 = Eq(dotV_2, dotV_3)
solve([e1, e2, e3], R_3, dict=True)

r1 = solve(Eq(simplify(e1.rhs / e2.rhs),1), R_3, dict=True)
r1[0][R_3]

Re = Symbol("Re", positive=True)
Re_K = Symbol("Re_K", positive=True)
rho = Symbol("\\rho", positive=True)
barv = Symbol("\\bar v", positive=True)
R = Symbol("R", positive=True)
e4 = Eq(Re,2*R*rho*barv/eta)

dotV = Symbol("\dot V", positive=True)
e5 = Eq(dotV, pi*R**2*barv)
solve(e5, barv)[0]
e4_1 = e4.subs(barv, solve(e5, barv)[0]).subs(Re, Re_K)

dotV_1 = Symbol("\dot V_1", positive=True)
R_1 = Symbol("R_1", positive=True)
e6 = Eq(dotV_1, 2*dotV_2)

solve(e4_1.subs(R, R_1).subs(dotV, dotV_1),R_1)[0]
solve(e4_1.subs(R, R_2).subs(dotV, dotV_2),R_2)[0]

#e4.subs([(R, R_1)])

#solve(e1, R_2, dict=True)
#solve(e1 / e2, R_3, dict=True)

x = Symbol('x', real=True)
a, b, c = symbols('a b c')
Cx0=Symbol('{C_{x_{0}}')
latex(Cx0)
