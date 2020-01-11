#/usr/bin/env python3
# __**__coding:utf-8 __*__

import os,sys
import numpy as np
from netCDF4 import Dataset
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
  print('python input_file')
  sys.exit()
else:
	fn = sys.argv[1]

f= Dataset(fn, 'r')
x = f.variables['x'][:]
y = f.variables['y'][:]
h = f.variables['rho'][:,:][0]

X, Y = np.meshgrid(x,y)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(X, Y, h)

ax.set_xlim([0,1])
ax.set_ylim([0,1])
ax.set_zlim([-0.5,0.5])

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.title('8th order')
plt.savefig('res.png', format='png', dpi=300, bbox_inches='tight')
plt.close(fig)
