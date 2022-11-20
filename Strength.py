import numpy as np
import matplotlib.pyplot as plt
import math
import sympy
import scipy
import mpmath

N = 100
lu = 1
R = 5
pi = math.pi

phi = lu * pi * R**2

value = -sympy.EulerGamma.n(10) + scipy.special.expi(phi) - math.log(phi, math.e)
# value = phi* (1 + phi/4 + phi**2/18 + phi**3/96 + phi**4/600)

som = 0
result = [som]

for n in range(1, N):
    som += (phi)**n/ (n * math.factorial(n))
    result.append(som)

fig, ax = plt.subplots()
plt.plot(range(N), result)
plt.axhline(value)
plt.show()

print(som - value)