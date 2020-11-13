#!/bin/sh

for i in 20 50 100 200 500 1000 2000 4000
do
cp ./1gamma.mac ./gamma"$i".mac
echo "/gun/direction 0 0 1" >> ./gamma"$i".mac
echo "/gun/position 0 0 -10.1 mm" >> ./gamma"$i".mac
echo "/gun/energy "$i" keV" >> ./gamma"$i".mac
echo "/run/beamOn 10000" >> ./gamma"$i".mac
done
