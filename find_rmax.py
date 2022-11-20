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
lu = 0.1
lbs = 0.01
w = 1

phi = lu * pi

def dR(R, N, p):
    term1 = (1 + (N-1)*(1-p))/(N * (1 + (N-1)*p))
    term2 = 2 * lbs/lu
    term3 = 2*N*phi*math.exp(-N*phi*R**2) - 1/R**2 * (1 - math.exp(-N*phi*R**2))
    return term1 * term2 * term3


def function(phi):
    return 2*pi* 1/math.sqrt(phi)
phis = [2, 3, 4, 5, 6, 7, 8]
vals = [4.47177, 3.65118, 3.16202,2.8282, 2.58178, 2.39026, 2.23588]

plt.scatter(phis, vals)
plt.scatter(phis, [function(phi) for phi in phis])
plt.show()