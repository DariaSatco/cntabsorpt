set terminal postscript eps size 15cm,10cm enhanced color font 'Helvetica,30' #--> put large font size!
set output 'im_sigm_perp_1010.eps'
set border linewidth 0.8

set style line 11 linecolor rgb '#D22D56' linewidth 3 linetype 1 #--> you can define your own line style
set style line 22 linecolor rgb '#666666' linewidth 2 linetype 0
set style line 33 linecolor rgb '#000080' linewidth 7 linetype 0
set style line 44 linecolor rgb '#90ee90' linewidth 3 linetype 1
set style line 55 linecolor rgb '#006400' linewidth 7 linetype 0
set key

#Set large ticks scale!
set tics scale 2.5
set tics front

#set yrange [-5:10]
#set ytics 0,0.05,0.2

#Don't forget axis properties!

set xzeroaxis ls 22

set xlabel 'Energy (eV)' offset 0.0,0.2
set mxtics 5
set ylabel 'Im {/Symbol s}_{/Symbol \136} (e_2/h)' offset 1.5,0.0
set mytics 5


#plot "tube.elecDOS.xyy.1611" u 2:1 w l ls 11, \
#"tube.elecDOS.xyy.1611" u 2:($2 > 0.88 ? $1 : 0):1 w filledcurves lc rgb '#cccccc'

plot "tube.sigm2_intra.xyy.101012" u 1:2 title 'intraband, E_F = 1.22 eV' w l ls 11, \
"tube.sigm2_inter.xyy.101012" u 1:2 title 'interband, E_F = 1.22 eV' w l ls 44
