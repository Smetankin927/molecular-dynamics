set xrange [0:50]
set yrange [-150:300]
set nokey
set terminal png size 600,500 enhanced font "Helvetica,10"
set output "Energies.png"
plot "Energies.txt" using 4:1 w p pt 5, "Energies.txt" using 4:2 w p pt 3, "Energies.txt" using 4:3 w p pt 2 
