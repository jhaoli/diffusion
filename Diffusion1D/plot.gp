
# Copyright(C) 2016 Western University
#
# Plot the numerical solution u(x,t) of the heat equation
#
#     u  = au  ,  x in [-L, L], u(x,0) = f(x), u(-L,t) = u(L,t) = 0.
#      t     xx
#
# using gnuplot.
set xrang[0:1]
set yrang[-0.5:0.5]
set nokey
set term X11

plot "output000.dat" using 1:2 w line
pause -1 "Hit return to exit."
plot "output001.dat" using 1:2 w line
pause -1 "Hit return to exit."
plot "output002.dat" using 1:2 w line
pause -1 "Hit return to exit."
plot "output003.dat" using 1:2 w line
pause -1 "Hit return to exit."
plot "output004.dat" using 1:2 w line
pause -1 "Hit return to exit."
plot "output005.dat" using 1:2 w line
pause -1 "Hit return to exit."
plot "output006.dat" using 1:2 w line
pause -1 "Hit return to exit."
plot "output007.dat" using 1:2 w line
pause -1 "Hit return to exit."
plot "output008.dat" using 1:2 w line
pause -1 "Hit return to exit."
plot "output009.dat" using 1:2 w line
pause -1 "Hit return to exit."
plot "output010.dat" using 1:2 w line
pause -1 "Hit return to exit."
