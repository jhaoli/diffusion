set term X11
set xrang[0:1]
set yrang[-0.5:0.5]

plot "output000.dat" using 1:2 w line, "output001.dat" using 1:2 w line, "output002.dat" using 1:2 w line, "output003.dat" using 1:2 w line, "output004.dat" using 1:2 w line,"output005.dat" using 1:2 w line,"output006.dat" using 1:2 w line,"output007.dat" using 1:2 w line,"output008.dat" using 1:2 w line,"output009.dat" using 1:2 w line,"output010.dat" using 1:2 w line
