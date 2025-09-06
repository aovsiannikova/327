#!/bin/sh

for i in 20 50 100 200 500 1000 2000 4000 
do
./OpNovice2 ./gammas/gamma"$i".mac 
cp ./results_h2_absXY.csv ./xy_q512/q512_abs_XY_h_"$i"keV.csv
cp ./results_h2_X_ev.csv ./xy_q512/q512_X_ev_"$i"keV.csv
cp ./results_nt_absorption.csv ./xy_q512/q512_abs_XY_"$i"keV.csv
cp ./results_nt_phot_count.csv ./xy_q512/q512_counts_"$i"keV.csv

./OpNovice2 ./gammas/2gamma"$i".mac 
cp ./results_h2_absXY.csv ./xy_q512_1/q512_1_abs_XY_h_"$i"keV.csv
cp ./results_h2_X_ev.csv ./xy_q512_1/q512_1_X_ev_"$i"keV.csv
cp ./results_nt_absorption.csv ./xy_q512_1/q512_1_abs_XY_"$i"keV.csv
cp ./results_nt_phot_count.csv ./xy_q512_1/q512_1_counts_"$i"keV.csv

./OpNovice2 ./gammas/3gamma"$i".mac 
cp ./results_h2_absXY.csv ./xy_q512_2/q512_2_abs_XY_h_"$i"keV.csv
cp ./results_h2_X_ev.csv ./xy_q512_2/q512_2_X_ev_"$i"keV.csv
cp ./results_nt_absorption.csv ./xy_q512_2/q512_2_abs_XY_"$i"keV.csv
cp ./results_nt_phot_count.csv ./xy_q512_2/q512_2_counts_"$i"keV.csv

#cp ./results_nt_phot_count.csv ./25/q511_25_counts_"$i"keV.csv
#cp ./results_nt_status.csv ./25/q511_25_abs_"$i"keV.csv
done


