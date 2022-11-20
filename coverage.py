import math

import matplotlib.pyplot as plt
import numpy as np

xmin, xmax = 0, 10
ymin, ymax = 0, 10
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
bandwidth = 0.1
w = 0.1
k = 3
p = 0.5

bs_points = np.random.uniform(xmin, xmax, (number_of_bs, 2))
axis = np.linspace(xmin, xmax, number_of_users)
# x, y = np.meshgrid(axis, axis)
# user_points = np.stack([x, y]).transpose(1,2,0).reshape(-1,2)
user_points = np.random.uniform(xmin, xmax, (number_of_users, 2))

x_bs = [bs[0] for bs in bs_points]
y_bs = [bs[1] for bs in bs_points]

lambda_BS = number_of_bs / (xDelta * yDelta)
lambda_U = number_of_users**2 / (xDelta * yDelta)

phi = lambda_U * pi * R ** 2

user_degrees = np.zeros(number_of_users)

for i in range(number_of_users):
    for j in range(number_of_bs):
        dist = distance(user_points[i], bs_points[j])
        if 0 < dist <= R:
            user_degrees[i] = 1

print((number_of_users - sum(user_degrees)) / number_of_users)
print('Theoretical:', (1 - pi * R**2/(xmax * ymax))**number_of_bs)

fig, ax = plt.subplots()
for bs in range(len(x_bs)):
    circle = plt.Circle(bs_points[bs], R, color='r', alpha=0.5)
    ax.add_patch(circle)
ax.set_box_aspect(1)
plt.scatter(x_bs, y_bs, marker='2', color='k')
# plt.scatter([u[0] for u in user_points], [u[1] for u in user_points])
plt.show()
