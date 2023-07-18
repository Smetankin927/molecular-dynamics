
set terminal png size 600,500 enhanced font "Helvetica,10"
set output "MSD.png"
stats "MSD.txt" name "A" nooutput

plot "MSD.txt" using 2:1 w l
