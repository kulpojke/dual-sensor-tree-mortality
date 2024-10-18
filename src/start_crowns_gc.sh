#!/bin/sh
# install R
#sudo apt-get update
#sudo apt -y install gdal-bin libgdal-dev libudunits2-dev
#sudo apt -y install r-base r-base-dev
#install R packages
#mkdir -p ~/R/library
#export R_LIBS_USER=~/R/library
#Rscript -e 'install.packages(c("plyr", "dplyr", "terra", "sf", "lidR", "future"), lib="~/R/library")'


# usage: ./start_crowns_gc.sh chd-trintiy-bucket 10

BUCKET="$1"
UTM_ZONE="$2"

x=0
function do_stuff {
    # make a list of the laz files in the bucket
    gcloud storage ls gs://$BUCKET/*.laz >>  bucket_laz.list

    # make a list of files with ttops.gpkg in the bucket
    gcloud storage ls gs://$BUCKET/*_ttops.gpkg >>  bucket_ttops.list

    # make a list of files with crowns.gpkg in the bucket
    gcloud storage ls gs://$BUCKET/*_crowns.gpkg >>  bucket_crowns.list

    # run the python script to compare lists
    python3 compare_file_lists.py \
        --laz_list=bucket_laz.list \
        --ttops_list=bucket_ttops.list \
        --crowns_list=bucket_crowns.list

    # above step writes a new file, todo.list
    export UTM_ZONE=$UTM_ZONE
    export BUCKET=$BUCKET

    cat todo.list | parallel \
        --bar \
        --joblog crowns_log.txt \
        -j 2 \
        --env UTM_ZONE --env BUCKET \
        'gcloud storage cp {} {/} && \
        gcloud storage cp {.}.tiff {/.}.tiff && \
        Rscript crown_delineation_gc.R {/} {/.}.tiff $UTM_ZONE && \
        gcloud storage cp crowns/{/.}_crowns.gpkg gs://$BUCKET && \
        gcloud storage cp ttops/{/.}_ttops.gpkg gs://$BUCKET && \
        gcloud storage cp {/.}_log.txt gs://$BUCKET && \
        rm {/.}_log.txt ttops/{/.}_ttops.gpkg crowns/{/.}_crowns.gpkg {/} {/.}.tiff'

    # when done relist the files
    rm bucket_laz.list
    rm bucket_ttops.list
    rm bucket_crowns.list
    rm todo.list
    gcloud storage ls gs://$BUCKET/*.laz >>  bucket_laz.list
    gcloud storage ls gs://$BUCKET/*_ttops.gpkg >>  bucket_ttops.list
    gcloud storage ls gs://$BUCKET/*_crowns.gpkg >>  bucket_crowns.list

    python3 compare_file_lists.py \
        --laz_list=bucket_laz.list \
        --ttops_list=bucket_ttops.list \
        --crowns_list=bucket_crowns.list

    LINES="$(cat todo.list | wc -l)"
    LAZ_LINES="$(cat bucket_laz | wc -l)"


    echo "" >> crowns_log.txt
    echo "------------------------------------" >> crowns_log.txt
    echo $(date) >> crowns_log.txt

    if [ "$LINES" -eq "0" ]; then
        echo "Finished: files are complete." >> crowns_log.txt
    else
        if [ $x -lt 3 ]; then
            if [ $1 -eq $LAZ_LINES ]; then # if more files have not shown up
                echo "Finished: but there seem to be new laz files." >> crowns_log.txt
                echo " Trying to run this again..." >> crowns_log.txt
                do_stuff $LAZ_LINES
                # increment x
                ((x++))
            else # if more files have shown up
                echo "Finished: More laz added, trying again."
                do_stuff $LAZ_LINES
                # reset tries
                x=0
            fi
        else
            echo "Finished: there seem to be new laz files, but have"  >> crowns_log.txt
            echo "tried calling self recursively 3 times." >> crowns_log.txt
        fi
    fi
}

do_stuff 0

