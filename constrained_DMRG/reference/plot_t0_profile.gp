# Reproduce the bundled PNG with: gnuplot plot_t0_profile.gp
datafile = "one_kink_t0_profile.csv"
set terminal pngcairo size 1500,950 enhanced font "Arial,16"
set output "one_kink_t0_profile.png"
set datafile separator ","
set multiplot layout 2,1 title "Accepted t=0 one-kink state: N=100, Delta=3, chi_max=64" font ",20"
set grid ytics lc rgb "#dddddd"
set border lw 1.1
set key opaque box
set xrange [1:100]
set arrow 1 from 50.5, graph 0 to 50.5, graph 1 nohead dt 2 lw 1.5 lc rgb "#666666"

set title "Target and measured local polarization (stage 0, after 45 sweep pairs)"
set ylabel "<Sz_i>"
set xlabel ""
set yrange [-0.52:0.52]
plot datafile using 1:2 with lines dt 2 lw 2 lc rgb "#777777" title "target <Sz_i>", \
     datafile using 1:3 with linespoints pt 7 ps 0.45 lw 1.7 lc rgb "#d95f02" title "DMRG <Sz_i>"

unset arrow 1
set arrow 2 from 50.5, graph 0 to 50.5, graph 1 nohead dt 2 lw 1.5 lc rgb "#666666"
set title "Local Lagrange multiplier from the same checkpoint"
set ylabel "lambda_i"
set xlabel "site i"
set yrange [*:*]
plot datafile using 1:4 with linespoints pt 7 ps 0.45 lw 1.7 lc rgb "#2878b5" title "lambda_i"
unset multiplot
unset output
