#!/opt/local/bin/gnuplot -persist
set term aqua 0
plot "test_ngrad.txt" using 1:2:4:5 with vectors
set term aqua 1
plot "test_ngrad.txt" using 1:2:7:8 with vectors
set term aqua 2
plot "test_ngrad.txt" using 1:2:10:11 with vectors
