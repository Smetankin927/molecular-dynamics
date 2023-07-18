set term gif animate
set terminal gif animate delay 0.00005
set output "maxwell.gif"
stats "Maxwell.txt" name "A" nooutput
set xrange [A_min_x:A_max_x]
set yrange [A_min_y - 1:A_max_y + 1]
do for [i=0 : int(A_blocks - 1)]{plot "Maxwell.txt" index i w p pt 7}
pause -1 "Hit any key to continue"
