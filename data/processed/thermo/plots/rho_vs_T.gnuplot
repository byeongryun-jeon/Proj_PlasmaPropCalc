set datafile separator ','
set terminal pngcairo size 1400,900 enhanced
set output '/home/brjeon/Proj_PlasmaPropCalc/data/processed/thermo/plots/rho_vs_T.png'
set grid
set key outside right top
set xlabel 'Temperature [K]'
set ylabel 'rho [kg/m^3]'
set title 'Density'
set logscale y
set format y '%.3e'
plot \
'/home/brjeon/Proj_PlasmaPropCalc/data/processed/thermo/argon_thermo_0p1_1_4atm.csv' u 1:(abs($2-0.1)<1e-12 ? $3 : 1/0) w l lw 2 lc rgb '#1f77b4' title 'This work 0.1 atm', \
'/home/brjeon/Proj_PlasmaPropCalc/data/processed/thermo/argon_matf_reference_0p1_1_4atm.csv' u 1:(abs($2-0.1)<1e-12 ? $6 : 1/0) w l lw 2 dt 2 lc rgb '#1f77b4' title 'MATF 0.1 atm', \
'/home/brjeon/Proj_PlasmaPropCalc/data/processed/thermo/argon_thermo_0p1_1_4atm.csv' u 1:(abs($2-1)<1e-12 ? $3 : 1/0) w l lw 2 lc rgb '#d62728' title 'This work 1 atm', \
'/home/brjeon/Proj_PlasmaPropCalc/data/processed/thermo/argon_matf_reference_0p1_1_4atm.csv' u 1:(abs($2-1)<1e-12 ? $6 : 1/0) w l lw 2 dt 2 lc rgb '#d62728' title 'MATF 1 atm', \
'/home/brjeon/Proj_PlasmaPropCalc/data/processed/thermo/argon_thermo_0p1_1_4atm.csv' u 1:(abs($2-4)<1e-12 ? $3 : 1/0) w l lw 2 lc rgb '#2ca02c' title 'This work 4 atm', \
'/home/brjeon/Proj_PlasmaPropCalc/data/processed/thermo/argon_matf_reference_0p1_1_4atm.csv' u 1:(abs($2-4)<1e-12 ? $6 : 1/0) w l lw 2 dt 2 lc rgb '#2ca02c' title 'MATF 4 atm'
