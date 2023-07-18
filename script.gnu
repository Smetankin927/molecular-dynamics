set term gif animate
set terminal gif animate delay 0.00005
set output "simulation.gif"
stats "data.txt" name "A" nooutput
set xrange [0:10]
set yrange [0:10]
set zrange [0:10]
do for [i=0 : int(A_blocks - 1)]{splot "data.txt" index i w p pt 7}
pause -1 "Hit any key to continue"
