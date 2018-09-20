set terminal postscript eps size 10cm,15cm enhanced color font 'Helvetica,30' #--> put large font size!
set output 'energy_160.eps'
set border linewidth 0.8

set style line 11 linecolor rgb '#D22D56' linewidth 3 linetype 1 #--> you can define your own line style
set style line 22 linecolor rgb '#666666' linewidth 2 linetype 0
unset key

#Set large ticks scale!
set tics scale 2.5
set tics front

#set xrange [-0.3:0.3]
#set yrange [-1.5:1.5]
#set ytics 0,0.05,0.2

#Don't forget axis properties!
set xlabel 'k_{||}' offset 0.0,0.2
set mxtics 5
#set ylabel '{/Symbol e}_2' offset 1.5,0.0
set ylabel 'E (eV)' offset 1.5,0.0
set mytics 5

set xzeroaxis ls 22

#plot "tube.elecDOS.xyy.1611" u 2:1 w l ls 11, \
#"tube.elecDOS.xyy.1611" u 2:($2 > 0.88 ? $1 : 0):1 w filledcurves lc rgb '#cccccc'

N = system("awk 'NR==1{print NF}' tube.Enk.xyy.160")
plot for [col=2:N] "tube.Enk.xyy.160" u 1:col w l ls 11
