#!/bin/sh
# usage:

BUCKET=$1
UTM_ZONE=$2

# make a list of the laz files in the bucket
gcloud storage ls gs://$BUCKET/*.laz >>  bucket_laz.list

# make a list of files with ttops.gpkg in the bucket
gcloud storage ls gs://$BUCKET/*_ttops.gpkg >>  bucket_ttops.list

# make a list of files with crowns.gpkg in the bucket
gcloud storage ls gs://$BUCKET/*_crowns.gpkg >>  bucket_crowns.list

# run the python script to compare lists
python compare_file_lists.py \
    --laz_list=bucket_laz.list \
    --ttops_list=bucket_ttops.list \
    --crowns_list=bucket_crowns.list

# above step writes a new file, todo.list
cat todo.list | parallel \
    --bar \
    --joblog gnu_parallel.log \
    -j 4 \
    gcloud storage cp gs://$BUCKET/{} {} && \
    Rscript crown_delineation_restart.R {} $UTM_ZONE


BASENAME="${INPUT_FILE%.*}"
CROWNS="${BASENAME}_crowns.gpkg"
TTOPS="${BASENAME}_ttops.gpkg"


# TODO:change R script to work on only 1 file, adjust outfile names
# TODO: add log
Rscript crown_delineation_restart.R /work /data /cold $utm_zone 

# TODO: tries and error handlingfor this

# just chuck em' in the bucket, do directories
gcloud storage cp $CROWNS gs://$BUCKET
gcloud storage cp $TTOPS gs://$BUCKET


