set terminal postscript eps size 11cm,15cm enhanced color font 'Helvetica,30' #--> put large font size!
set output 'epsilon_105.eps'
set border linewidth 0.8

set style line 11 linecolor rgb '#D22D56' linewidth 3 linetype 1 #--> you can define your own line style
unset key

#Set large ticks scale!
set tics scale 2.5
set tics front

set multiplot layout 2,1
unset title
set ylabel '{/Symbol e}_1 (Kubo)' offset 1.5,0.0
set mytics 5
set xlabel 'Energy (eV)' offset 0.0,0.2
set mxtics 5
unset key
plot 'tube.eps1.xyy.105' using 1:2 w l ls 11
#
set xlabel 'Energy (eV)' offset 0.0,0.2
set mxtics 5
set ylabel '{/Symbol e}_2 (Kubo)' offset 1.5,0.0
set mytics 5
unset key
plot 'tube.eps2.xyy.105' using 1:2 w l ls 11

unset multiplot
