from sympy import symbols, Matrix, Eq, solve, nsolve, linsolve, nonlinsolve
from pprint import pprint

# Define the variables
x, y, a = symbols('x y a')

# Define the matrices A and H
A = Matrix([[1.1163353, 0.1196577], [-0.05982886, 0.99667761]])
H = Matrix([[558.84381, -18.29773], [-18.29773, 624.81490]])
pprint(A)
pprint(H)

# Define x0
x0 = Matrix([0.1, 0.2])

Nmin = 9
h = [A.T * H * A]  # Initial value of h

# Calculate the values of h
for i in range(1, Nmin):
    h.append(A.T * h[i - 1] * A)

psi0, psi1 = symbols('psi0 psi1')
psi = Matrix([psi0, psi1])
pprint(psi)
pprint(h)
print(":FKSJD:FKJSDF")
len(h)

fun = Matrix([0, 0])
for i in range(Nmin):
    bottom = (psi.T * h[i].inv() * psi)**0.5
    fun += ((h[i].inv() * psi) / bottom[0])

# bottom = (psi.T * h[0].inv() * psi)**0.5
# fun = (h[0].inv() * psi) / bottom[0]
pprint(fun)
# Solve the system of equations
eq1 = Eq(-x0[0]/a, fun[0])
eq2 = Eq(-x0[1]/a, fun[1])
eq3 = Eq(psi[0]**2 + psi[1]**2, 1)

sol = solve((eq1, eq2, eq3), (a, psi0, psi1))

pprint(sol)
