reset

set logscale xy

set format y "%.1e"

set xrange [80:*]

set key bottom left

set xlabel 'T(M)'
set ylabel '|{/Symbol Y}|

p 'psi.1.extract.asc' u 1:(abs($2)) w l  t 'l=0, Horizon'
rep 'psi.1.extract.asc' u 1:(abs($4)) w l  t 'l=0, Scri+'
rep 'psi.2.extract.asc' u 1:(abs($2)) w l  t 'l=1, Horizon'
rep 'psi.2.extract.asc' u 1:(abs($4)) w l  t 'l=1, Scri+'
rep 'psi.3.extract.asc' u 1:(abs($2)) w l  t 'l=2, Horizon'
rep 'psi.3.extract.asc' u 1:(abs($4)) w l  t 'l=2, Scri+'

rep 12*x**(-2) w l dt 2 t 't^{-2}'
rep 100*x**(-3) w l dt 2 t 't^{-3}'
rep 13*x**(-3) w l dt 2 lt 8 notitle
rep 12*x**(-4) w l dt 2 lt 9 t 't^{-4}'
rep 200*x**(-5) w l dt 2 lt 10 t 't^{-5}'
