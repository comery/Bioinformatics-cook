import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Data for plotting
#t = np.arange(0.0, 2.0, 0.01)
#s = 1 + np.sin(2 * np.pi * t)
t = np.arange(0, 208, 1)
s = []
with open(sys.argv[1], 'r') as fh:
    for i in fh:
        tmp = i.strip().split("\t")
        #t.append(tmp[0])
        s.append(int(tmp[1]))

fig, ax = plt.subplots()
ax.plot(t, s)

ax.set(xlabel='coverage (X)', ylabel='genome site number',
       title='The coverage distribution of spotted hyena genome')
ax.grid()

fig.savefig("hyena_depth.png", dpi=300)
plt.show()
