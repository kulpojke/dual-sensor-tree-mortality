#!/bin/sh
# usage:
# ./start_crowns_rocker.sh /path/to/laz/dir /path/to/output_dir

laz_volume=$1
out_volume=$2

#docker build docker -t wrtc_rocker_geo && \
docker run --rm -v $PWD:/code -v $out_volume:/work -v $laz_volume:/data -w /work  wrtc_rocker_geo \
    Rscript /code/crown_delineation.R /work /data
