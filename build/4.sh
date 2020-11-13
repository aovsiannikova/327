#!/bin/sh

for i in 20 50 100 200 500 1000 2000 4000 
do
./OpNovice2 .gamma"$i".mac 
cp ./results_h2_absXY.csv ./xy_q512/q512_abs_XY_h_"$i"keV.csv
done


