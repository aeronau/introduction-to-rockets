from scipy import *
from scipy.optimize import fsolve

F = 1e3
R0 = 8.314
Tc = 1500
M = 14.4e-3
k = 1.4
pa = 1e5
pe = 3e5
pc = 100e5

def A2diam(A):
    return 2 * sqrt(A/pi)

epsilon = sqrt((k - 1)/2) * (2 / (k+1))**((k+1) / 2 / (k-1)) / (pe/pc)**(1/k) / sqrt(1 - (pe/pc)**((k - 1)/k))
Ve = sqrt( 2*k*R0*Tc / (k-1) / M * (1 - (pe/pc)**((k-1)/k)))
K = sqrt(k) * (2 / (k+1))**((k+1) / 2 / (k-1))
Cstar = sqrt(R0 / M * Tc) / K;

def eqs(x):
    eq1 = -epsilon**2 + 1/x[1]**2 * (2 / (k+1) * (1 + (k-1) / 2 * x[1]**2 ))**((k+1) / (k-1))
    eq2 = -F + x[0] * Ve + x[2] * epsilon * (pe - pa)
    eq3 = -x[0] + pc * x[2] / Cstar
    return [eq1, eq2, eq3]

x = fsolve(eqs, array([0.4838, 2.9355, 9.2e-3]))

print("epsilon =", epsilon)
print("Ve =", Ve)
print("Cstar =", Cstar)
print('mp =', x[0]) 
print('Me =', x[1])
print('dt =', A2diam(x[2]))
print('de =', A2diam(epsilon * x[2]))
