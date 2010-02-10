# !/usr/local/bin/gnuplot

set terminal gif
set xlabel "Simulation Time (cycles)"

set ylabel "Variance"
set output 'conns.gif' 
plot 'variance_sqrt.txt' using 1:2 with lines ti "sqrt", \
     'variance_fixed.txt' using 1:2 with lines ti "fixed", \
     'variance_random.txt' using 1:2 with lines ti "random", \
     'variance_full.txt' using 1:2 with lines ti "full"

set ylabel "Number of Connections"
set output 'variance.gif' 
plot 'variance_sqrt.txt' using 1:3 with lines ti "sqrt", \
     'variance_fixed.txt' using 1:3 with lines ti "fixed", \
     'variance_random.txt' using 1:3 with lines ti "random", \
     'variance_full.txt' using 1:3 with lines ti "full"

