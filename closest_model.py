import math

import numpy as np

xmin, xmax = 0, 1
ymin, ymax = 0, 1
xDelta = xmax - xmin
yDelta = ymax - ymin


def distance(coord1, coord2):
    x = np.minimum((coord1[0] - coord2[0]) % xDelta, (coord2[0] - coord1[0]) % xDelta)
    y = np.minimum((coord1[1] - coord2[1]) % yDelta, (coord2[1] - coord1[1]) % yDelta)
    # x = coord1[0] - coord2[0]
    # y = coord1[1] - coord2[1]
    return math.sqrt(x ** 2 + y ** 2)


def find_squared_distance(x, y, xbs, ybs):
    x = np.minimum((x - np.array(xbs)) % xDelta, (np.array(xbs) - x) % xDelta)
    y = np.minimum((y - np.array(ybs)) % yDelta, (np.array(ybs) - y) % yDelta)
    return (x ** 2 + y ** 2)


def find_bs(xpop, ypop, xbs, ybs):
    indices = find_squared_distance(xpop, ypop, xbs, ybs).argsort()
    return indices[0]


pi = math.pi
bandwidth = 100
w = bandwidth
k = 3
p = 1

number_of_bs = 100
number_of_users = k * 1000

lambda_BS_t = number_of_bs / (xDelta * yDelta)
lambda_U_t = number_of_users/k / (xDelta * yDelta)
lambda_U = number_of_users / (xDelta * yDelta)

lambda_BS = lambda_BS_t + (k-1)*(1-p)*lambda_BS_t
number_of_bs = int(lambda_BS * (xDelta * yDelta))

bs_points = np.random.uniform(xmin, xmax, (number_of_bs, 2))
user_points = np.random.uniform(xmin, xmax, (number_of_users, 2))

xbs = [i for i, j in bs_points]
ybs = [j for i, j in bs_points]



def fR(r):
    return 2 * (lambda_BS * pi * r ** 2) / r * math.exp(-lambda_BS * pi * r ** 2)


def FR(r):
    return 1 - math.exp(-lambda_BS * pi * r ** 2)


BS_degrees = np.zeros(number_of_bs)
user_degrees = np.zeros(number_of_users)
distances = []
Rinv = []

for i in range(number_of_users):
    j = find_bs(user_points[i][0], user_points[i][1], xbs, ybs)
    dist = distance(user_points[i], bs_points[j])
    distances.append(dist)
    Rinv.append(1 / dist)
    BS_degrees[j] += 1

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
        Dinv.append(1 / BS_degrees[j])

for i in range(number_of_users):
    j = find_bs(user_points[i][0], user_points[i][1], xbs, ybs)
    dist = distance(user_points[i], bs_points[j])
    DR.append(BS_degrees[j] * dist)
    strengths.append(Ws[j] / (BS_degrees[j] * dist))

print('E(W) =', (k - 1) * p * w + w, 'real =', sum(Ws) / len(Ws))
print('E(1/D) =', lambda_BS / lambda_U, 'real =', sum(Dinv) / len(Dinv))
print('E(1/R) =', math.sqrt(lambda_BS) * pi, 'real =', sum(Rinv) / len(Rinv))
print('E(S) =', ((k - 1) * p + 1) * w * lambda_BS / lambda_U * math.sqrt(lambda_BS) * pi, 'real=',
      sum(strengths) / len(strengths))

print('K=', w * lambda_BS_t ** (3 / 2) * pi / lambda_U_t)
print('K root(N)=', w * lambda_BS_t ** (3 / 2) * pi / lambda_U_t * math.sqrt(k))
