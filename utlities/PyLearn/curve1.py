import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
 
def func(x, a, b, c):
    return a * np.exp(-b * x) + c

xdata =[]
ydata =[]
with open("a",'r') as fh:
    lines = fh.readlines()
    for line in lines:
        a = line.split(',')
        xdata.append(int(a[0]))
        ydata.append(int(a[1]))


plt.plot(xdata, ydata, 'b-', label='data')
 
# Fit for the parameters a, b, c of the function `func`
popt, pcov = curve_fit(func, xdata, ydata)
plt.plot(xdata, func(xdata, *popt), 'r-', label='fit')
 
# Constrain the optimization to the region of ``0 < a < 3``, ``0 < b < 2``
# and ``0 < c < 1``:
popt, pcov = curve_fit(func, xdata, ydata, bounds=(0, [3., 2., 1.]))
plt.plot(xdata, func(xdata, *popt), 'g--', label='fit-with-bounds')
 
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()
