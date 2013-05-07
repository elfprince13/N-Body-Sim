#!/opt/local/bin/gnuplot -persist
set term aqua 0
#splot "test_grav.txt" using 1:10:11 with lines, "test_grav.txt" using 1:4:5 with lines, "test_grav.txt" using 1:16:17 with lines
set xrange [0:1e06]
plot "data/precision/gdata-0-1000.0-10000000.0-1-0(ForwardEuler)-0(Newtonian).txt" using 1:(log10($4/$2)) with lines, "data/precision/gdata-0-1000.0-10000000.0-1-0(ForwardEuler)-1(BarnesHut).txt" using 1:(log10($4/$2)) with lines, "data/precision/gdata-0-1000.0-10000000.0-1-1(RK4)-0(Newtonian).txt" using 1:(log10($4/$2)) with lines, "data/precision/gdata-0-1000.0-10000000.0-1-2(SPRK6Ruth)-0(Newtonian).txt" using 1:(log10($4/$2)) with lines, "data/precision/gdata-0-1000.0-10000000.0-1-3(SPRK6SuzukiTrotter)-0(Newtonian).txt" using 1:(log10($4/$2)) with lines
unset xrange

set term aqua 1
set xrange [0:1e03]
#plot "test_grav_dt1.txt" using 1:(log10($3/$2)) with lines,"test_grav_dt.1.txt" using 1:(log10($3/$2)) with
plot "rk1test.txt" using ($1/1000):(log10($3/$2)) with lines,"rk4test.txt" using ($1/1000):(log10($3/$2)) with lines,"sprk6testRuth.txt" using ($1/1000):(log10($3/$2)) with lines,"sprk6testSuzukiTrotter.txt" using ($1/1000):(log10($3/$2)) with lines
#plot "data/precision/gdata-0-10000.0-100000000.0-1-0(ForwardEuler)-0(Newtonian).txt" using 1:(log10($4/$2)) with lines, "data/precision/gdata-0-10000.0-100000000.0-1-0(ForwardEuler)-1(BarnesHut).txt" using 1:(log10($4/$2)) with lines, "data/precision/gdata-0-10000.0-100000000.0-1-1(RK4)-0(Newtonian).txt" using 1:(log10($4/$2)) with lines, "data/precision/gdata-0-10000.0-100000000.0-1-2(SPRK6Ruth)-0(Newtonian).txt" using 1:(log10($4/$2)) with lines, "data/precision/gdata-0-10000.0-100000000.0-1-3(SPRK6SuzukiTrotter)-0(Newtonian).txt" using 1:(log10($4/$2)) with lines
unset xrange
set term aqua 2
#set xrange [0:1e06]
#plot "test_grav_dt1.txt" using 1:($3/$2) with lines,"test_grav_dt.1.txt" using 1:($3/$2) with lines,"test_grav_dt10.txt" using 1:($3/$2) with lines, 10**-6 with lines

plot "data/precision/gdata-0-100.0-1000000.0-1-0(ForwardEuler)-0(Newtonian).txt" using ($1/100):(log10($4/$2)) with lines, "data/precision/gdata-0-100.0-1000000.0-1-0(ForwardEuler)-1(BarnesHut).txt" using ($1/100):(log10($4/$2)) with lines, "data/precision/gdata-0-1000.0-10000000.0-1-0(ForwardEuler)-0(Newtonian).txt" using ($1/1000):(log10($4/$2)) with lines, "data/precision/gdata-0-1000.0-10000000.0-1-0(ForwardEuler)-1(BarnesHut).txt" using ($1/1000):(log10($4/$2)) with lines, "data/precision/gdata-0-10000.0-100000000.0-1-0(ForwardEuler)-0(Newtonian).txt" using ($1/10000):(log10($4/$2)) with lines, "data/precision/gdata-0-10000.0-100000000.0-1-0(ForwardEuler)-1(BarnesHut).txt" using ($1/10000):(log10($4/$2)) with lines

set term aqua 3
plot "data/softspeed_newtonian.csv" using (log10($1)):(log10($2)):(log10($2-$3)):(log10($2+$3)) with yerrorbars,"data/softspeed_barneshut.csv" using (log10($1)):(log10($2)):(log10($2-$3)):(log10($2+$3)) with yerrorbars

set term aqua 4

f1(x) = a1*x**2 + b1*x + c1
fit f1(x) "data/softspeed_newtonian.csv" using 1:(1/$2):(abs(1/($2-$3))) via a1,b1,c1

f2(x) = a2 + b2 * (x * log(x)/log(4))
fit f2(x) "data/softspeed_barneshut.csv" using 1:(1/$2):(abs(1/($2-$3))) via a2,b2

plot "data/softspeed_newtonian.csv" using 1:(1/$2):(1/($2-$3)):(1/($2+$3)) with yerrorbars,"data/softspeed_barneshut.csv" using 1:(1/$2):(1/($2-$3)):(1/($2+$3)) with yerrorbars,f1(x) title "O(N**2)",f2(x) title "O(N log4(N))"
