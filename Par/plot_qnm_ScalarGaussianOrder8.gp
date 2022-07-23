reset

set logscale y

set format y "%.1e"

set xrange [0:400]
set yrange [1e-21:100]

set key top right

set xlabel 'T(M)'
set ylabel '|{/Symbol Y}|

p 'psi.1.extract.asc' u 1:(abs($2)) w l  t 'l=0'
rep 'psi.2.extract.asc' u 1:(0.01*abs($2)) w l  t 'l=1'
rep 'psi.3.extract.asc' u 1:(0.0001*abs($2)) w l  t 'l=2'

rep 12*abs(exp(-0.1048957201*x)*sin(0.1104549399*x+1.82)) t 'expected l=0'
rep 0.06*abs(exp(-0.09765998891*x)*sin(0.2929361333*x+pi/2)) t 'expected l=1'
rep 4.5e-4*abs(exp(-.09675877598*x)*sin(0.4836438722*x+1.1)) lt 7 t 'expected l=2'

