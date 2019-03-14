#/usr/bin/env python3
# __**__coding:utf-8 __*__

import os,sys
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)

for i in np.arange(10):
	data = np.loadtxt("output0%02d.dat" % i)
	x = data[:,0]
	y = data[:,1]

	ax.plot(x,y)

	ax.set_xlim([0,1])
	ax.set_ylim([-0.5,0.5])

	ax.set_xlabel('X')
	ax.set_ylabel('Y', rotation=0)
major_ticks = np.arange(0,1,0.5)
minor_ticks = np.arange(0,1,0.02)
ax.set_xticks(minor_ticks, minor=True)

plt.title('8th order')
plt.savefig('diff.png', format='png', dpi=200, bbox_inches='tight')
plt.close(fig)
