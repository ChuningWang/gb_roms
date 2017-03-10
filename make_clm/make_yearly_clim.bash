#!/bin/bash

fyear=2000
lyear=2000

dst_dir=/Volumes/R1/scratch/chuning/gb_roms/data/clim/

for year in $(seq $fyear $lyear ) ; do

    python make_clim_file.py $year $dst_dir

    cd $dst_dir 
    ncrcat GlacierBay_clim_${year}_??_SODA3.3.1.nc -o GlacierBay_clim_SODA3.3.1_y${year}.nc
    rm GlacierBay_clim_${year}_??_SODA3.3.1.nc
    cd ~/git/gb_roms/make_clm

done

