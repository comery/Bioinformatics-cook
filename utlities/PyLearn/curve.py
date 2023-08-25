import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
'''#定义x、y散点坐标
x = np.arange(1, 16, 1)
num = [4.00, 5.20, 5.900, 6.80, 7.34,
       8.57, 9.86, 10.12, 12.56, 14.32,
       15.42, 16.50, 18.92, 19.58, 20.00]
'''
def func(x, a, b):
    return a * np.exp(b/x)
x =[]
y =[]
with open("a",'r') as fh:
    lines = fh.readlines()
    for line in lines:
        a = line.split(',')
        x.append(int(a[0]))
        y.append(int(a[1]))
    
#非线性最小二乘法拟合
popt, pcov = curve_fit(func, x, y)
#获取popt里面是拟合系数
a = popt[0] 
b = popt[1]
yvals = func(x,a,b) #拟合y值
print("u'系数a:',"+ str(a) )
print("u'系数b:',"+ str(b))
 
#绘图
plot1 = plt.plot(x, y, 's',label='original values')
plot2 = plt.plot(x, yvals, 'r',label='polyfit values')
plt.xlabel('x')
plt.ylabel('y')
plt.legend(loc=4) #指定legend的位置右下角
plt.title('curve_fit')
plt.show()
plt.savefig('test2.png')
