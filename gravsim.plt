#!/opt/local/bin/gnuplot -persist
set term aqua 0
#splot "test_grav.txt" using 1:10:11 with lines, "test_grav.txt" using 1:4:5 with lines, "test_grav.txt" using 1:16:17 with lines

set term aqua 1
#set xrange [0:1e06]
#plot "test_grav_dt1.txt" using 1:(log10($3/$2)) with lines,"test_grav_dt.1.txt" using 1:(log10($3/$2)) with lines,"test_grav_dt10.txt" using 1:(log10($3/$2)) with lines, -6 with lines
plot "rk1test.txt" using 1:(log10($3/$2)) with lines,"rk4test.txt" using 1:(log10($3/$2)) with lines,"sprk6test.txt" using 1:(log10($3/$2)) with lines,"sprk6test2.txt" using 1:(log10($3/$2)) with lines

set term aqua 2
#set xrange [0:1e06]
#plot "test_grav_dt1.txt" using 1:($3/$2) with lines,"test_grav_dt.1.txt" using 1:($3/$2) with lines,"test_grav_dt10.txt" using 1:($3/$2) with lines, 10**-6 with lines

plot "rk1test.txt" using 1:($3/$2) with lines,"rk4test.txt" using 1:($3/$2) with lines,"sprk6test.txt" using 1:($3/$2) with lines,"sprk6test2.txt" using 1:($3/$2) with lines
