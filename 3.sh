#!/bin/sh

for i in 20 50 100 200 500 1000 2000 4000
do
cp ./1gamma.mac ./gammas/gamma"$i".mac
echo "/gun/direction 0 0 1" >> ./gammas/gamma"$i".mac
echo "/gun/position 0 0 -10.1 mm" >> ./gammas/gamma"$i".mac
echo "/gun/energy "$i" keV" >> ./gammas/gamma"$i".mac
echo "/run/beamOn 10000" >> ./gammas/gamma"$i".mac

cp ./1gamma.mac ./gammas/2gamma"$i".mac
echo "/gun/direction 0 0 1" >> ./gammas/2gamma"$i".mac
echo "/gun/position 0.71 0.71 -10.1 mm" >> ./gammas/2gamma"$i".mac
echo "/gun/energy "$i" keV" >> ./gammas/2gamma"$i".mac
echo "/run/beamOn 10000" >> ./gammas/2gamma"$i".mac

cp ./1gamma.mac ./gammas/3gamma"$i".mac
echo "/gun/direction 0 0 1" >> ./gammas/3gamma"$i".mac
echo "/gun/position 1.41 1.41 -10.1 mm" >> ./gammas/3gamma"$i".mac
echo "/gun/energy "$i" keV" >> ./gammas/3gamma"$i".mac
echo "/run/beamOn 10000" >> ./gammas/3gamma"$i".mac

done
