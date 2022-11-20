import math
import sympy
import scipy

import matplotlib.pyplot as plt
import mpmath
import numpy as np

xmin, xmax = 0, 100
ymin, ymax = 0, 100
xDelta = xmax - xmin
yDelta = ymax - ymin

def distance(coord1, coord2):
    x = np.minimum((coord1[0] - coord2[0]) % xDelta, (coord2[0] - coord1[0]) % xDelta)
    y = np.minimum((coord1[1] - coord2[1]) % yDelta, (coord2[1] - coord1[1]) % yDelta)
    # x = coord1[0] - coord2[0]
    # y = coord1[1] - coord2[1]
    return math.sqrt(x ** 2 + y ** 2)

number_of_bs = 100
number_of_users = 5000
R = 1
pi = math.pi
bandwidth = 1
w = bandwidth
k = 2
p = 0.5

bs_points = np.random.uniform(xmin, xmax, (number_of_bs, 2))
user_points = np.random.uniform(xmin, xmax, (number_of_users, 2))

lambda_BS = number_of_bs / (xDelta * yDelta)
lambda_U = number_of_users / (xDelta * yDelta)

phi = lambda_U * pi * R ** 2

def fR(r):
    if r <= R:
        return 2 * r / R ** 2
    else:
        return 0

def FR(r):
    if 0 <= r <= R:
        return r ** 2 / R ** 2
    else:
        return 1

def fD(d):
    return (lambda_U * pi * R ** 2) ** d / math.factorial(d) * math.exp(-phi)

def FD(d):
    return mpmath.gammainc(1 + d, lambda_U * pi * R ** 2) / mpmath.gamma(1 + d)

def fC(j):
    j = j - 1
    return math.comb(k - 1, j) * p ** j * (1 - p) ** (k - 1 - j)

def FC(v):
    som = 0
    for j in range(1, v + 1):
        som += fC(j)
    return som

def FDR(s):
    if 0 < s <= R+ scipy.special.expi(phi) - math.log(phi, math.e):
        return math.exp(-phi) * (1 + lambda_U * pi * s ** 2 * mpmath.hyper((1, 1, 1), (2, 2, 2), phi))
    elif s > R:
        y = math.ceil(s / R)
        return mpmath.gammainc(y, phi) / mpmath.gamma(y) + 1 / (y ** 2 * math.factorial(y)) * (
                math.exp(-phi) * phi ** y * R ** (-2) * s ** 2 * mpmath.hyper((1, y, y), (1 + y, 1 + y, 1 + y),
                                                                              phi))
    else:
        return math.exp(-phi)

def fDR(x):
    return ((2 * x/R**2)*(-sympy.EulerGamma + scipy.special.expi(phi) - math.log(phi, math.e))+1)*math.exp(-phi)

def FS_definition(s):
    som = 0
    if s > 0:
        for j in np.arange(1, k + 1, 1):
            som += (1 - FDR(w * j / s)) * fC(j)
        return som
    else:
        return 0

def weighted_expectation_S(r, N, p):
    lbs = number_of_bs / (xmax * ymax)
    return coverage_probability(r, N, p, lbs) * w / (phi * r ** 3) * 1 / N * ((N - 1) * p + 1) * (
                1 - math.exp(-N * phi * r ** 2))

def expectation_S(r, N, p):
    return w / (phi * r ** 3) * 1 / N * ((N - 1) * p + 1) * (1 - math.exp(-N * phi * r ** 2))


def coverage_probability(r, N, p, lbs):
    lbs = (1 + (N - 1) * (1 - p)) * lbs
    xmax, ymax = 1000, 1000
    area = xmax * ymax
    nbs = lbs * area
    return (1 - pi * r ** 2 / area) ** nbs

BS_degrees = np.zeros(number_of_bs)
user_degrees = np.zeros(number_of_users)
distances = []

for i in range(number_of_users):
    for j in range(number_of_bs):
        dist = distance(user_points[i], bs_points[j])
        if 0 < dist <= R:
            distances.append(dist)
            BS_degrees[j] += 1
            user_degrees[i] += 1

alpha = 1
strengths = []
Ws = []
Cs = []
DR = []
Dinv = []

for j in range(number_of_bs):
    n = np.random.binomial(k - 1, p)
    Ws.append(bandwidth * (n + 1))
    Cs.append(n + 1)
    if BS_degrees[j] > 0:
        Dinv.append(1/BS_degrees[j])

for i in range(number_of_users):
    for j in range(number_of_bs):
        dist = distance(user_points[i], bs_points[j])
        if dist <= R:
            DR.append(BS_degrees[j] * dist)
            strengths.append(Ws[j] / (BS_degrees[j] * dist))


# xs = np.arange(0, max(DR), 1)
# # plt.plot(xs, [FS_definition(r) for r in xs])
# plt.hist(strengths,  density=True, cumulative=False)
# plt.xlabel('strengths')
# plt.show()

print(sum(strengths)/len(strengths))
print(weighted_expectation_S(R, k, p))
print(expectation_S(R, k, p))

print('bandwidth:', sum(Ws)/len(Ws), w + (k-1)*p*w)
print('1/degree:', sum(Dinv)/len(Dinv), 1/(k * lambda_U * math.pi * R**2) * (1 - math.exp(-k*lambda_U * math.pi * R**2)))