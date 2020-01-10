#/usr/bin/env python3
# __**__coding:utf-8 __*__

import os,sys
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

data = np.loadtxt('output000.dat')

x = data[:,0]
y = data[:,1]
z = data[:,2]
nx = np.int(np.sqrt(np.shape(x)[0]))
ny = nx
X = np.reshape(x,(nx,ny))
Y = np.reshape(y,(nx,ny))
Z = np.reshape(z,(nx,ny))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(X, Y, Z)

ax.set_xlim([0,1])
ax.set_ylim([0,1])
ax.set_zlim([-0.5,0.5])

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.title('8th order')
plt.savefig('diff.png', format='png', dpi=200, bbox_inches='tight')
plt.close(fig)
