import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 3:
    print("usage: python3 {} <data.txt> <column>".format(sys.argv[0]))
    exit()

lie = sys.argv[2]
#if type(lie) is not 'int':
#    print("the second parameter must be int")
#    exit()

x = []
with open(sys.argv[1], "r") as f:
    for i in f:
        tmp = i.strip().split("\t")
        x.append(tmp[int(lie)-1])


plt.title("Depth of Hyena genome mapping")
plt.hist(x, 50)
plt.xlabel("depth")
plt.ylabel("count")
plt.show()
