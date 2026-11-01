# Render the version-controlled reference figures from locally generated outputs:
#   gnuplot kinks-in-real-space/render_reference_profiles.gp

set terminal pngcairo size 1800,820 font "Arial,22"
set datafile separator comma
set key outside top center horizontal
set grid ytics lc rgb "#d8d8d8"
set xrange [1:100]
set yrange [-0.52:0.52]
set xlabel "site i"
set ylabel "<S_i^z>"
set zeroaxis lw 1 lc rgb "#555555"

set output "kinks-in-real-space/figures/one_kink_reference.png"
set title "One kink: restricted core optimization, N = 100, {/Symbol D} = 3"
set object 1 rect from 45,-0.52 to 56,0.52 fc rgb "#f2c14e" fs transparent solid 0.22 behind
set arrow 1 from 50.5,-0.52 to 50.5,0.52 nohead dt 2 lw 2 lc rgb "#444444"
plot "kinks-in-real-space/outputs/one-kink-accurate/profiles.csv" \
     using ($1==20 ? $2 : 1/0):3 with linespoints pt 7 ps 0.55 lw 2.2 lc rgb "#df6c00" title "optimized A|B state", \
     NaN with boxes fc rgb "#f2c14e" fs transparent solid 0.22 title "optimized core: sites 45-56"
unset object 1
unset arrow 1

set output "kinks-in-real-space/figures/two_kinks_reference.png"
set title "Two kinks: restricted defect-block optimization, N = 100, {/Symbol D} = 3"
set arrow 1 from 36.5,-0.52 to 36.5,0.52 nohead dt 2 lw 2 lc rgb "#444444"
set arrow 2 from 68.5,-0.52 to 68.5,0.52 nohead dt 2 lw 2 lc rgb "#444444"
plot "kinks-in-real-space/outputs/two-kink-verify/profiles.csv" \
     using ($1==7 ? $2 : 1/0):3 with linespoints pt 7 ps 0.55 lw 2.2 lc rgb "#377eb8" title "optimized A|B|A state", \
     NaN with lines dt 2 lw 2 lc rgb "#444444" title "detected walls: 36|37, 68|69"
unset arrow 1
unset arrow 2
