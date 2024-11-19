#!/bin/python

# usage:
# GPKG=/home/michael/Work/dual-sensor-tree-mortality/input_files/tile_grid.gpkg
# python3 src/pdal_pipeline_gc.py \
# --grid_gpkg=$GPKG \
# --wrk_dir=output \

# optionally:
# --chm=True
# --bucket_name=chd-trintiy-bucket # yes, the bucket is mispelled!

# TODO:
# * add --ignore_tiles to args: list of tiles to ignore,
#   so as to not unecessarilly redo existing tiles 
# * make smooth_chm() more flexible and make --kernel arg
#   to pass to it.
#  * add --remove_local flag to remove local copies of files

#%%
import pdal
import geopandas as gpd
from tqdm import tqdm
import numpy as np
from skimage.morphology import disk
from skimage.filters import rank
import rasterio
from google.cloud import storage

import argparse
from pathlib import Path
from datetime import datetime
from ast import literal_eval
from time import sleep

from joblib import Parallel, delayed
import textwrap


# %%
def test_args():
    '''Makes args for testing. To use change for local env'''
    args = argparse.Namespace()
    args.grid_gpkg = '/home/michael/Work/dual-sensor-tree-mortality/input_files/output/metadata/tile_grid.gpkg'
    args.n_jobs = 10
    args.resolution = 0.5
    args.bucket_name = 'chd-trintiy-bucket'
    args.wrk_dir=Path('../output')
    args.tif_dir = str(args.wrk_dir / 'tif')
    args.laz_dir = str(args.wrk_dir / 'laz')
    args.log_dir = str(args.wrk_dir / 'log')
    return args


def parse_args():
    '''parses the arguments, returns args'''
    # init parser
    parser = argparse.ArgumentParser(
        prog='pdal_pipeline_gc.py',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''\
        Using contents of tile_grid.gpkg creates and excutes PDAL pipelines
        that:
            + download USGS tiles from urls, filtering noise and witheld
              points
            + merges and clips to new tile boundary
            + writes and smoothes CHM
            + writes normalized pc to laz
            + uploads tiff and laz to gc bucket (optional).

        Skips files that already exist in laz_dir.  Assumes that if files
        with the same stem name as tiles exist in a directory called 'crowns'
        next to laz_dir, the tiles have already been processed and will be skipped.    
         ''')
    )

    parser.add_argument(
        '--grid_gpkg',
        type=str,
        required=True,
        help='path to(custom)tile index file [required]'
    )

    parser.add_argument(
        '--wrk_dir',
        type=str,
        required=True,
        help='working directory [required]'
    )

    parser.add_argument(
        '--bucket_name',
        type=str,
        required=False,
        default=None,
        help='gc bucket name if upload desired [optional]'
    )

    parser.add_argument(
        '--resolution',
        type=float,
        required=False,
        default=0.5,
        help='resolution of chm [optional, default 0.5]'
    )

    parser.add_argument(
        '--n_jobs',
        type=int,
        required=False,
        default=10,
        help='jobs to run in parallel [optional, default 10]'
    )

    parser.add_argument(
        '--chm',
        type=bool,
        required=False,
        default=False,
        help='''
        Bool - If True: normalize pc and make chm,
            if False just write unnormalised pc.
            [optional, default False]
        '''
    )

    parser.add_argument(
        '--slow_your_roll',
        type=int,
        required=False,
        default=1000,
        help='''int - if more than this many laz files are present,
        wait for 30 minutes before processing next tile.
        [optional, default 1000]'''
    )

    args = parser.parse_args()

    # add aditional args
    args.wrk_dir = Path(args.wrk_dir)
    if args.chm:
        args.tif_dir = str(args.wrk_dir / 'tif')
    args.laz_dir = str(args.wrk_dir / 'laz')
    args.log_dir = str(args.wrk_dir / 'log')

    return args


def now():
    '''Returns formatted datetime of now.'''
    return datetime.now().strftime('%Y-%m-%d %H:%M:%S')


def write_to_log():
    '''Writes contents of, and clears, LOG.'''
    with open(Path(ARGS.log_dir) / 'pipeline_log.txt', 'a') as dst:
        _ = [dst.write(msg + '\n') for msg in LOG]
    LOG.clear()


def make_readers(list_of_urls):
    '''Retruns list of basic readers pointing at url'''
    readers = {}
    for url in list_of_urls:
        readers[Path(url).stem.split('_')[-1]] = pdal.Reader.las(url)
    return readers


def execute_single_reader(usgs_tile, reader):
    '''Returns a tuple with point array and message dict'''
    # make filter to remove:
    # noise (7), high-noise (18), misc crap ( > 19), and witheld
    filter = pdal.Filter.range(
        limits='Classification[0:6],Classification[8:17],Withheld[0:0]'
        )
    
    # make pipeline
    pipe = reader.pipeline() | filter

    # execute pipeline with exception handling
    try:
        n = pipe.execute()
        if n > 0:
            message = f'{now()} DOWNLOD SUCCESS: {usgs_tile} with {n} points  '
        else:
            message = f'{now()} EMPTY TILE: {usgs_tile} downloaded with {n} points  '
        arr = pipe.arrays[0]
        srs = pipe.metadata['metadata']['readers.las']['srs']['wkt']

    except Exception as e:
        n = -1
        message = f'{now()} PDAL FAILED: {e}  '
        arr = np.array([np.nan])

    return (arr, {
        'n': n,
        'message': message,
        'srs': srs
        })


def execute_readers(readers):
    '''Execute readers in parallel'''
    return Parallel(n_jobs=ARGS.n_jobs)(
        delayed(execute_single_reader)(usgs_tile, reader)
        for usgs_tile, reader in readers.items()
    )

def single_hag(arr, n_points):
    '''
    Runs hag_nn filter on arr --> hagged,
    returns tuple with hagged and message dict'''
    pipe = pdal.Filter.hag_nn().pipeline(arr)
    try:
        n = pipe.execute()
        if n > 0:
            message = f'{now()} HAG SUCCESS: {n} of {n_points} points made it through.  '
        else:
            message = f'{now()} HAG FAILED: {n_points} points were lost in HAG filter  '
        new_arr = pipe.arrays[0]

    except Exception as e:
        n = -1
        message = f'{now()} HAG FAILED - PDAL pipeline error: {e}  '
        new_arr = np.array([np.nan])

    return (new_arr, {
        'n': n,
        'message': message
        })


def execute_hag_filters(arrs):
    '''Exectues HAG filters in parallel.'''
    return Parallel(n_jobs=ARGS.n_jobs)(
        delayed(single_hag)(arr, len(arr))
        for arr in arrs
    )


def write_tif(arr, tile):
    '''
    Writes hag to tiff using pdal,
    returns message dict.'''
    tiff = Path(ARGS.tif_dir) / f'{tile}_.tiff'

    pipe = pdal.Writer.gdal(
        filename=tiff,
        resolution=ARGS.resolution,
        default_srs=ARGS.srs,
        dimension='HeightAboveGround',
        output_type=['mean']
    ).pipeline(arr)

    try:
        n = pipe.execute()    
        message = f'{now()} TIFF SUCCESS: {tiff}'
        new_arr = np.array([np.nan])

    except Exception as e:
        n = -1
        message = f'{now()} TIFF FAILED - PDAL pipeline error: {e}  '
        new_arr = np.array([np.nan])

    return [(new_arr, {
        'n': n,
        'message': message
        })]


def write_laz(arr, tile):
    '''
    Writes arr to laz using pdal,
    returns message dict.'''
    laz = Path(ARGS.laz_dir) / f'{tile}.laz'

    pipe = pdal.Writer.las(
        filename=laz,
        a_srs=ARGS.srs
    ).pipeline(arr)

    try:
        n = pipe.execute()    
        message = f'{now()} LAZ SUCCESS: {laz}'
        new_arr = np.array([np.nan])

    except Exception as e:
        n = -1
        message = f'{now()} LAZ FAILED - PDAL pipeline error: {e}  '
        new_arr = np.array([np.nan])

    return [(new_arr, {
        'n': n,
        'message': message
        })]


def smooth_chm(input_file, output_file):
    '''Smoothes CHM with 3x3 mean filter, writes new tiff.'''
    try:
        # open input raster
        with rasterio.open(input_file) as src:
            data = src.read(1, masked=True)
            profile = src.profile

        nodata = profile.get('nodata', None)
        data_filled = np.where(data.mask, 0, data)
        data_filled[data_filled < 0] = 0

        # scale data to 0-255
        data_min = np.nanmin(data_filled)
        data_max = np.nanmax(data_filled)
        data_scaled = (data_filled - data_min) / (data_max - data_min)
        data_scaled = (data_scaled * 255).astype(np.uint8)
        
        # mean filter
        filtered_scaled = rank.mean(data_scaled, disk(1))

        # scale data back to original range
        filtered_data = filtered_scaled.astype(np.float32) / 255
        filtered_data = filtered_data * (data_max - data_min) + data_min

        # restore nodata values
        #filtered_data = np.where(np.isnan(data_filled), nodata, filtered_data)

        # update profile
        profile.update(dtype=rasterio.float32, nodata=nodata)

        # write
        with rasterio.open(output_file, 'w', **profile) as dst:
            dst.write(filtered_data.astype(np.float32), 1)
        
        # messages
        n = -1
        message = f'{now()} SMOOTHING SUCCESS: {output_file}  '
        new_arr = np.array([np.nan])

    except Exception as e:
        n = -1
        message = f'{now()} SMOOTHING FAILED: {e}  '
        new_arr = np.array([np.nan])

    return [(new_arr, {
        'n': n,
        'message': message
        })]


def extract_results(results, set_srs=False):
    '''
    Wraps functions, extracts point arrays, adds log messages to LOG,
    returns lists of arrays
    '''
    arrs = []
    for thing in results:
        if thing[1]['n'] > 0:
            arrs.append(thing[0])
        LOG.append(thing[1]['message'])

    if set_srs:
        ARGS.srs = thing[1]['srs'] 

    return arrs


def upload_to_gcs(bucket_name, src, dst):
    '''Uploads a file to Google Cloud Storage.'''
    try:
        # initialize client
        storage_client = storage.Client()

        # get the bucket
        bucket = storage_client.bucket(bucket_name)

        # create a blob
        blob = bucket.blob(dst)

        # Upload the file
        blob.upload_from_filename(src)

        LOG.append(f'{now()} UPLOAD SUCCESS: {src} uploaded to {bucket_name}/{dst}')
    except Exception as e:
        LOG.append(f'{now()} UPLOAD FAILED: {e}')

#%%

# make args global, they will be used and modified by functions
ARGS = parse_args()

# ensure output dirs exist
Path(ARGS.laz_dir).mkdir(parents=True, exist_ok=True)
Path(ARGS.log_dir).mkdir(parents=True, exist_ok=True)
if ARGS.chm:
    Path(ARGS.tif_dir).mkdir(parents=True, exist_ok=True)
    
# see if files are already there (for restarts)
if not ARGS.chm:
    existing_files = [laz.stem for laz in Path(ARGS.laz_dir).glob('*.laz')]
else:
    existing_files = [
        tiff.stem
        for tiff
        in Path(ARGS.tif_dir).glob('*.tiff')
        if not tiff.stem.endswith('_')
        ]

# also check the crown dir for file that were processed in R and deleted
crown_dir = Path(ARGS.laz_dir).parent / 'crowns'
existing_crowns = [cr.stem for cr in crown_dir.glob('*.gpkg')]

# combine existig crowns with exiting laz
existing_files = list(set(existing_files + existing_crowns))

grid = gpd.read_file(ARGS.grid_gpkg)
grid['tiles'] = grid.tiles.apply(literal_eval)

for _, row in tqdm(grid.iterrows(), total=len(grid)):
    try:
        # get tile id
        tile = row['name']
        # start a log list and and tile id
        LOG = [
            '------------------------------------------  ',
            f'# {tile}  '
            ]
    
        # check existing to avoid rework
        if tile in existing_files:
            LOG.append(f'{now()} SKIPPING: {tile} has already been processed')
            write_to_log()
            continue

        # read points
        readers = make_readers(row.tiles)
        points = extract_results(execute_readers(readers), set_srs=True)
        write_to_log()

        # if points is an empty list, write a log message and skip t onext iterr
        if len(points) == 0:
            LOG.append(
                f'{now()} TILE FAILED: No points in any USGS tiles within {tile}.'
                )
            write_to_log()
            continue

        # if chm: hag nn filter points
        if ARGS.chm:
            points = extract_results(execute_hag_filters(points))
            write_to_log()

        # concat all the arrays
        lumped = np.concatenate(points)
        write_to_log()

        # if chm: write chm
        if ARGS.chm:
            # write tiff
            _ = extract_results(write_tif(lumped, tile))
            input_file = Path(ARGS.tif_dir) / f'{tile}_.tiff'

            # smooth tiff
            output_file = Path(ARGS.tif_dir) / f'{tile}.tiff'
            _ = extract_results(smooth_chm(input_file, output_file))

            # delete unsmoothed chm if smoothed succesfully
            if 'SMOOTHING SUCCESS' in LOG[-1]:
                input_file.unlink()
            write_to_log()

            # ferry HeightAboveGround to Z
            lumped['Z'] = lumped['HeightAboveGround']

        # write  pc
        _ = extract_results(write_laz(lumped, tile))
        write_to_log()

        # upload to bucket
        if ARGS.bucket_name is not None:
            src = Path(ARGS.tif_dir) / f'{tile}.tiff'
            dst = f'{tile}.tiff' 
            upload_to_gcs(ARGS.bucket_name, src, dst)
            write_to_log()

            src = Path(ARGS.laz_dir) / f'{tile}.laz'
            dst = f'{tile}.laz' 
            upload_to_gcs(ARGS.bucket_name, src, dst)
            write_to_log()
        else:
            LOG.append(f'{now()} SKIPPING UPLOAD: {tile} - no bucket name provided.')
            write_to_log()

        # put spaces below block of message for tile
        _ = [LOG.append(' ') for i in range(2)]
        write_to_log()
    except Exception as e:
         LOG.append(f'{now()} FAILED: unknown reason not caugth by error handling:  {e}')
         write_to_log()

    laz_count = len([laz for laz in Path(ARGS.laz_dir).glob('*.laz')])

    # if more than ARGS.slow_your_roll laz files have piled up, sleep 30 min
    if  laz_count >= ARGS.slow_your_roll:
        message = f'{now} WAITING: {laz_count} laz files in laz_dir. Waiting 30 minutes.'
        LOG.append(message)
        write_to_log()
        print(message)
        sleep(30 * 60)
#%%

        

