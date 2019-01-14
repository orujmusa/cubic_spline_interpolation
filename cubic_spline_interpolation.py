import numpy as np
from math import sqrt

def cubic_spline_interpolation(x0, x, y):

x = np.asfarray(x)
y = np.asfarray(y)
size = len(x)
xdiff = np.diff(x)
ydiff = np.diff(y)
Li = np.empty(size)
Li_1 = np.empty(size-1)
z = np.empty(size)
Li[0] = sqrt(2*xdiff[0])
Li_1[0] = 0
B0 = 0
z[0] = B0 / Li[0]

for i in range(1, size-1, 1):
Li_1[i] = xdiff[i-1] / Li[i-1]
Li[i] = sqrt(2*(xdiff[i-1]+xdiff[i]) - Li_1[i-1] * Li_1[i-1])
Bi = 6*(ydiff[i]/xdiff[i] - ydiff[i-1]/xdiff[i-1])
z[i] = (Bi - Li_1[i-1]*z[i-1])/Li[i]
i = size - 1
Li_1[i-1] = xdiff[-1] / Li[i-1]
Li[i] = sqrt(2*xdiff[-1] - Li_1[i-1] * Li_1[i-1])
Bi = 0
z[i] = (Bi - Li_1[i-1]*z[i-1])/Li[i]
i = size-1
z[i] = z[i] / Li[i]

for i in range(size-2, -1, -1):
z[i] = (z[i] - Li_1[i-1]*z[i+1])/Li[i]
index = x.searchsorted(x0)
np.clip(index, 1, size-1, index)
xi1, xi0 = x[index], x[index-1]
yi1, yi0 = y[index], y[index-1]
zi1, zi0 = z[index], z[index-1]
hi1 = xi1 - xi0
f0 = zi0/(6*hi1)*(xi1-x0)**3 + \
zi1/(6*hi1)*(x0-xi0)**3 + \
(yi1/hi1 - zi1*hi1/6)*(x0-xi0) + \
(yi0/hi1 - zi0*hi1/6)*(xi1-x0)
return f0

if __name__ == '__main__':
import matplotlib.pyplot as plt
x = np.linspace(0.9, 13.3, num=21)
y = np.array([1.3, 1.5, 1.85, 2.1, 2.6, 2.7, 2.4, 2.15, 2.05, 2.1, 2.25, 2.3, 2.25, 1.95,
1.4, 0.9, 0.7, 0.6, 0.5, 0.4, 0.25])
plt.scatter(x, y, color = 'red')
x_new = np.linspace(0.9, 13.3, num=50)
plt.plot(x_new, cubic_spline_interpolation(x_new, x, y), 'k--')
plt.ylim(-1,5)
plt.title('Cubic Spline Interpolation')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.ylim(-1,5)
plt.show()
