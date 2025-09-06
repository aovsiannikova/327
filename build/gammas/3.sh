#!/bin/sh

for i in 20 50 100 200 500 1000 2000 4000
do
cp ./1gamma.mac ./gamma"$i".mac
echo "/gun/direction .46631 0 1" >> ./gamma"$i".mac
echo "/gun/position 0 0 -5.1 mm" >> ./gamma"$i".mac
echo "/gun/energy "$i" keV" >> ./gamma"$i".mac
echo "/run/beamOn 10000" >> ./gamma"$i".mac

cp ./1gamma.mac ./2gamma"$i".mac
echo "/gun/direction 1.191754 0 1" >> ./2gamma"$i".mac
echo "/gun/position 0 0 -5.1 mm" >> ./2gamma"$i".mac
echo "/gun/energy "$i" keV" >> ./2gamma"$i".mac
echo "/run/beamOn 10000" >> ./2gamma"$i".mac

cp ./1gamma.mac ./3gamma"$i".mac
echo "/gun/direction 0 0 1" >> ./3gamma"$i".mac
echo "/gun/position 0.71 0.71 -5.1 mm" >> ./3gamma"$i".mac
echo "/gun/energy "$i" keV" >> ./3gamma"$i".mac
echo "/run/beamOn 10000" >> ./3gamma"$i".mac

cp ./1gamma.mac ./4gamma"$i".mac
echo "/gun/direction 0 0 1" >> ./4gamma"$i".mac
echo "/gun/position 1.41 1.41 -5.1 mm" >> ./4gamma"$i".mac
echo "/gun/energy "$i" keV" >> ./4gamma"$i".mac
echo "/run/beamOn 10000" >> ./4gamma"$i".mac

done
