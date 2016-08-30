set format '$%g$'
set terminal epslatex standalone color
set xlabel "$\\Wavenumber (cm^{-1})$"
set ylabel "$|K|^2$"
set style line 1 lw 3
set output "shift.tex"
plot [0:2000] "shift.dat" u 1:2 w i
