#!/bin/bash

fyear=2000
lyear=2000

grd_name=GlacierBay_usgs
dst_dir=/Volumes/R1/scratch/chuning/gb_roms/clm/

for year in $(seq $fyear $lyear ) ; do

    python make_clim_file.py $year $dst_dir

    cd $dst_dir 
    ncrcat -O ${grd_name}_clim_${year}_??_SODA3.3.1.nc -o ${grd_name}_clim_${year}_SODA3.3.1.nc
    rm ${grd_name}_clim_${year}_??_SODA3.3.1.nc
    cd ~/git/gb_roms/make_clm

done

