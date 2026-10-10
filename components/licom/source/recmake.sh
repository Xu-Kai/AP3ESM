#!/usr/bin/bash


#LICOM_RES="1km"
LICOM_RES="2km"
LICOM_RES="10km-fcst"
LICOM_RES="10km"
cp def-undef.h-${LICOM_RES} $PWD/sunway/def-undef.h
cp def-undef.h-${LICOM_RES} $PWD/host-sunway/def-undef.h
cp def-undef.h-${LICOM_RES} def-undef.h
cp host-sunway/paramsw.h sunway/paramsw.h
rm -rf build/
mkdir build
cd build
cmake ..
make -j8
pwd
cp licom3 /home/export/online1/mdt00/shisuan/swiap/bjh/run-10km-obp 
