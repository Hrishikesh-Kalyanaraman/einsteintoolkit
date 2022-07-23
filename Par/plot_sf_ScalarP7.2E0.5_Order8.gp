reset

set logscale y

set format y "%.1e"

set xlabel 'T (M)'
set ylabel '|F_r|

p 'Schwarzschild_wave.frl.asc' u 1:(abs($2)) w l t 'l=0'
rep '' u 1:(abs($3)) w l t 'l=1'
rep '' u 1:(abs($4)) w l t 'l=2'
rep '' u 1:(abs($5)) w l t 'l=3'
rep '' u 1:(abs($6)) w l t 'l=4'
rep '' u 1:(abs($7)) w l t 'l=5'
rep '' u 1:(abs($8)) w l t 'l=6'
rep '' u 1:(abs($9)) w l t 'l=7'
rep '' u 1:(abs($10)) w l dt 2  t 'l=8'
rep '' u 1:(abs($11)) w l dt 2  t 'l=9'
rep '' u 1:(abs($12)) w l dt 2  t 'l=10'
