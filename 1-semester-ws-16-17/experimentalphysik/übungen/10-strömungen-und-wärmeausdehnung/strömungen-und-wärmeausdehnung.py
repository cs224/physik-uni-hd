
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


#solve(e1, R_2, dict=True)
#solve(e1 / e2, R_3, dict=True)

x = Symbol('x', real=True)
a, b, c = symbols('a b c')
Cx0=Symbol('{C_{x_{0}}')
latex(Cx0)
