import math

import matplotlib.pyplot as plt
import numpy as np

colors = ['blueviolet', 'dodgerblue', 'mediumseagreen', 'deeppink', 'coral', 'royalblue', 'midnightblue',
          'yellowgreen', 'darkgreen', 'mediumblue', 'DarkOrange', 'green', 'red', 'MediumVioletRed',
          'darkcyan', 'orangered', 'purple', 'cornflowerblue', 'saddlebrown', 'indianred', 'fuchsia', 'DarkViolet',
          'salmon', 'sandybrown', 'seagreen', 'seashell', 'sienna', 'silver', 'skyblue', 'slateblue',
          'slategray', 'snow', 'springgreen', 'steelblue', 'tan', 'teal', 'thistle', 'black', 'grey',
          'tomato', 'turquoise', 'violet', 'wheat', 'whitesmoke', 'yellow']

pi = math.pi
lu = 0.01
lbs = 0.01
w = 1

phi = lu * pi


def expectation_S(r, N, p):
    return w / (phi * r ** 3) * 1 / N * ((N - 1) * p + 1) * (1 - math.exp(-N * phi * r ** 2))


def check_max(n):
    phi = lu * pi
    return (1 - math.exp(n * phi) + n - math.exp(-phi) * n) / (1 - math.exp(n * phi) - n ** 2 + math.exp(-phi) * n ** 2)


def weighted_expectation_S(r, N, p):
    return (1 + (N - 1) * (1 - p)) / (N * (1 + (N - 1) * p)) * (lbs * w) / (lu * r) * (1 - math.exp(-N * phi * r ** 2))


def coverage_probability(r, N, p, lbs):
    lbs = (1 + (N - 1) * (1 - p)) * lbs
    xmax, ymax = 1000, 1000
    area = xmax * ymax
    nbs = lbs * area
    return 1 - (1 - pi * r ** 2 / area) ** nbs


R = np.arange(0.01, 10, 0.0001)
delta = 0.05
ps = np.arange(0, 1 + delta, delta)

N = 2
print(weighted_expectation_S(10, 1, 0), weighted_expectation_S(10, 5, 0))
fig, ax = plt.subplots()
plt.plot(ps, [weighted_expectation_S(1, 1, 0) for r in ps], label='N=1', color=colors[0])
plt.plot(ps, [weighted_expectation_S(1, 2, p) for p in ps], label=f'N={N}, p = 0', color=colors[1])
plt.plot(ps, [weighted_expectation_S(1, 3, p) for p in ps], label=f'N={N}, p = 0.5', color=colors[2])
plt.plot(ps, [weighted_expectation_S(1, 4, p) for p in ps], label=f'N={N}, p = 1', color=colors[3])
plt.legend()
# plt.yscale('log')
plt.ylabel('$E(S)$')
plt.xlabel('Rmax')
plt.show()


ns = np.arange(1, 20, 1)

fig, ax = plt.subplots()
plt.scatter(ns, [weighted_expectation_S(3, n, 0) for n in ns], label = 'p=0', color = colors[0])
plt.scatter(ns, [weighted_expectation_S(3, n, 0.1) for n in ns], label = 'p=0.1', color = colors[1])
plt.scatter(ns, [weighted_expectation_S(3, n, 0.2) for n in ns], label = 'p=0.2', color = colors[2])
plt.scatter(ns, [weighted_expectation_S(3, n, 0.5) for n in ns], label = 'p=0.5', color = colors[3])
plt.scatter(ns, [weighted_expectation_S(3, n, 1) for n in ns], label = 'p=1', color = colors[4])

plt.xlabel('N')
plt.legend()
plt.show()

print( [weighted_expectation_S(3, n, 0) for n in ns])