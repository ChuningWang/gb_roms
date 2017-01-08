#!/bin/bash

fyear=1995
lyear=2014

for year in $(seq $fyear $lyear ) ; do

    python make_clim_file.py $year

    cd /Volumes/R1/scratch/chuning/data/gb_roms/clim
    ncrcat GlacierBay_clim_${year}_??_SODA3.3.1.nc -o GlacierBay_clim_SODA3.3.1_y${year}.nc
    rm GlacierBay_clim_${year}_??_SODA3.3.1.nc
    cd ~/git/gb_roms/make_clm

done

