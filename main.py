import numpy as np
import matplotlib.pyplot as plt
import math
import sympy

def distance_distribution(x, lambda1, lambda2, p):
    return 1 - math.exp(-lambda1 * math.pi * x**2) * math.exp(-(1-p)*lambda2 * math.pi*x**2)

def expected_capacity(lambda1, lambda2, p):
    lambda_tilde = lambda1 + (1-p)*lambda2
    alpha = 2
    c = 10**3
    phi = lambda_tilde * math.pi * c**(2/alpha)
    return lambda_tilde/2 * alpha/2 * (math.log2(phi) - sympy.EulerGamma - lambda_tilde * math.pi)

xs = np.arange(0, 2, 0.1)

lambda1, lambda2 = 0.05, 0.05

fig, ax = plt.subplots()
plt.plot(xs, [distance_distribution(x, lambda1, lambda2, 0) for x in xs], label = 'p = 0')
plt.plot(xs, [distance_distribution(x, lambda1, lambda2, 0.5) for x in xs], label = 'p = 0.5')
plt.plot(xs, [distance_distribution(x, lambda1, lambda2, 1) for x in xs], label = 'p = 1')
plt.legend()
plt.xlabel('Distance to closest BS')
plt.ylabel('CDF')
plt.show()


fig, ax = plt.subplots()
ps = np.arange(0, 1, 0.05)
plt.plot(ps, [expected_capacity(lambda1, lambda2, p) for p in ps], label = '$\lambda_1 = \lambda_2$')
plt.plot(ps, [expected_capacity(2* lambda2, lambda2, p) for p in ps], label = '$\lambda_1 = 2\lambda_2$')
plt.plot(ps, [expected_capacity(lambda2, 2 * lambda2, p) for p in ps], label = '$2\lambda_1 = \lambda_2$')
# plt.plot(ps, [expected_capacity(4 * lambda2, lambda2, p) for p in ps], label = '$\lambda_1 = 4\lambda_2$')
# plt.plot(ps, [expected_capacity(5 * lambda2, lambda2, p) for p in ps], label = '$\lambda_1 = 5\lambda_2$')
# plt.plot(ps, [expected_capacity(4*lambda2, lambda2, 0) for p in ps], ':', label = '$4 \cdot \lambda_2$')
# plt.plot(ps, [expected_capacity(5*lambda2, lambda2, 0) for p in ps], ':', label = '$5 \cdot \lambda_2$')
plt.legend()
plt.xlabel('$p$')
plt.ylabel('$E(C)$')
plt.show()

lambdas = np.arange(0.001, 1/math.pi, 0.1)
fig, ax = plt.subplots()
plt.plot(lambdas, [expected_capacity(l, 0, 0) for l in lambdas])
plt.show()