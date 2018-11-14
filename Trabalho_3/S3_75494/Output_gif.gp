set terminal gif animate delay 1
set output 'Output_gif.gif'
stats 'gif.dat' nooutput
set xrange [-1:1]
set yrange [-1.5:1.5]

do for [i=1:int(STATS_blocks)] {

t = system(sprintf("awk 'NR==%d{print $1}' '%s'", i-1, 'time.dat'))
set label 1 sprintf("t=%.3f", t+0) at screen 0.2,0.9
    plot 'gif.dat' index (i-1) using 1:2 with lines notitle
}