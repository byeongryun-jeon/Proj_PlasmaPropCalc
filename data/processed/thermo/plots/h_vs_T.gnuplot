set datafile separator ','
set terminal pngcairo size 1400,900 enhanced
set output '/home/brjeon/Proj_PlasmaPropCalc/data/processed/thermo/plots/h_vs_T.png'
set grid
set key outside right top
set xlabel 'Temperature [K]'
set ylabel 'h [J/kg]'
set title 'Specific Enthalpy'
set format y '%.3e'
plot \
'/home/brjeon/Proj_PlasmaPropCalc/data/processed/thermo/argon_thermo_0p1_1_4atm.csv' u 1:(abs($2-0.1)<1e-12 ? $4 : 1/0) w l lw 2 lc rgb '#1f77b4' title '0.1 atm', \
'/home/brjeon/Proj_PlasmaPropCalc/data/processed/thermo/argon_thermo_0p1_1_4atm.csv' u 1:(abs($2-1)<1e-12 ? $4 : 1/0) w l lw 2 lc rgb '#d62728' title '1 atm', \
'/home/brjeon/Proj_PlasmaPropCalc/data/processed/thermo/argon_thermo_0p1_1_4atm.csv' u 1:(abs($2-4)<1e-12 ? $4 : 1/0) w l lw 2 lc rgb '#2ca02c' title '4 atm'
