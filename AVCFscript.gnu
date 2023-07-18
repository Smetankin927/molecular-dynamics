set xrange [0:50]
set yrange [-1:3.5]
set terminal png size 600,500 enhanced font "Helvetica,10"
set output "AVCF.png"
plot "AVCF.txt" using 2:1 w l
