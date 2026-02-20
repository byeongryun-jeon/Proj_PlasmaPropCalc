set datafile separator ','
set terminal pngcairo size 1500,950 enhanced
set output '/home/brjeon/Proj_PlasmaPropCalc/data/processed/transport_properties/plots/kappa_components_4atm.png'
set grid
set key outside right top
set xlabel 'Temperature [K]'
set ylabel 'kappa component [W/(m K)]'
set title 'Thermal Conductivity Components (4 atm)'
set format y '%.3e'
plot \
'/home/brjeon/Proj_PlasmaPropCalc/data/processed/transport_properties/argon_transport_0p1_1_4atm.csv' u 1:(abs($2-4)<1e-12 ? $4 : 1/0) w l lw 2 lc rgb '#000000' title 'kappa total', \
'/home/brjeon/Proj_PlasmaPropCalc/data/processed/transport_properties/argon_transport_0p1_1_4atm.csv' u 1:(abs($2-4)<1e-12 ? $5 : 1/0) w l lw 2 lc rgb '#1f77b4' title 'k_H', \
'/home/brjeon/Proj_PlasmaPropCalc/data/processed/transport_properties/argon_transport_0p1_1_4atm.csv' u 1:(abs($2-4)<1e-12 ? $6 : 1/0) w l lw 2 lc rgb '#2ca02c' title 'k_e', \
'/home/brjeon/Proj_PlasmaPropCalc/data/processed/transport_properties/argon_transport_0p1_1_4atm.csv' u 1:(abs($2-4)<1e-12 ? $7 : 1/0) w l lw 2 lc rgb '#ff7f0e' title 'k_int', \
'/home/brjeon/Proj_PlasmaPropCalc/data/processed/transport_properties/argon_transport_0p1_1_4atm.csv' u 1:(abs($2-4)<1e-12 ? $8 : 1/0) w l lw 2 lc rgb '#d62728' title 'k_reac'
