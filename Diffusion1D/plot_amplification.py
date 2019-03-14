#/usr/bin/env python3
# __**__coding:utf-8 __*__

import os,sys
import numpy as np
import matplotlib.pyplot as plt


fig = plt.figure()

ax = fig.add_subplot(111)

dx = np.pi/100
x = np.linspace(0,np.pi,100)
a2 = dx**2/(2-2*np.cos(np.pi))

y2 = 1 - a2*(2-2*np.cos(x))/dx**2

a4 = dx**4/(6-8*np.cos(np.pi)+2*np.cos(2*np.pi))
y4 =  1 - a4*(6-8*np.cos(x)+2*np.cos(2*x))/dx**4

a6 = dx**6/(20-30*np.cos(np.pi)+12*np.cos(2*np.pi)-2*np.cos(3*np.pi))
y6 = 1 - a6*(20-30*np.cos(x)+12*np.cos(2*x)-2*np.cos(3*x))/dx**6

ax.plot(x, y2, label = '2nd')
ax.plot(x, y4, label = '4th')
ax.plot(x, y6, label = '6th')


ax.set_xlim(0,np.pi)
ax.set_ylim(0,1)

label_xmajor = ['0', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$']

major_ticks = np.arange(0, np.pi+np.pi/4, np.pi/4)
minor_ticks = np.arange(0, np.pi+np.pi/20, np.pi/20)

ax.set_xticks(minor_ticks, minor=True)
ax.set_xticks(major_ticks)
ax.set_xticklabels(label_xmajor)

ax.grid(which='major', linestyle='dotted')
ax.set_xlabel(r'$k\Delta x$')
ax.set_ylabel(r'$|G|^2$', rotation=0)
ax.legend(frameon=True)

plt.savefig('amplification.png', format='png', dpi=200, bbox_inches='tight')
plt.close(fig)

