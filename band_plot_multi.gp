set terminal postscript eps size 19cm,11cm enhanced color font 'Helvetica,30' #--> put large font size!
set output 'ph_en_dos.eps'
set border linewidth 0.8

set style line 11 linecolor rgb '#D22D56' linewidth 3 linetype 1 #--> you can define your own line style
unset key

#Set large ticks scale!
set tics scale 2.5
set tics front

#Don't forget axis properties!
set xlabel 'DOS (states/atom/eV)' offset 0.0,0.2
set mxtics 5
set ylabel 'Energy (eV)' offset 1.5,0.0
set mytics 5

set multiplot layout 1,2
#set bmargin 5
#
set title "E(q)"
unset key
set xlabel 'q ({/Symbol p}/T)' offset 0.0,0.2
set mxtics 5
set ylabel 'Energy (eV)' offset 1.5,0.0
set mytics 5
set yrange [-0.01:0.201]

N = system("awk 'NR==1{print NF}' tube.Emq.xyy.default")
plot for [col=2:N] "tube.Emq.xyy.default" u 1:col w l ls 11

#
set title "DOS"
unset key
set xlabel 'DOS (states/atom/eV)' offset 0.0,0.2
set xtics 0,20,90
set mxtics 5
set ylabel 'Energy (eV)' offset 1.5,0.0
set mytics 5
set yrange [-0.01:0.201]
plot "tube.phDOS.xyy.default" u 1:2 w l ls 11

unset multiplot
