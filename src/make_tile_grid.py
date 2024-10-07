#!bin/python3

# usage:
# python3 make_tile_grid.py \
# --tindex=/home/michael/Work/dual-sensor-tree-mortality/input_files/lpc_index.gpkg \
# --urls=/home/michael/Work/dual-sensor-tree-mortality/input_files/chd_urls.txt
# --grid_size=1000 \
# --buffer=50 \


#%%
import geopandas as gpd
from pathlib import Path
import argparse
from shapely import Polygon
from joblib import Parallel, delayed
from tqdm import tqdm
import numpy as np
import pandas as pd
import pdal

from google.cloud import storage

import requests
import os
from datetime import datetime
import json
from ast import literal_eval

#%%


def parse_arguments():
    '''parses the arguments, returns args'''
    # init parser
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--tindex',
        type=str,
        required=True,
        help='path to tile index file'
    )

    parser.add_argument(
        '--urls',
        type=str,
        required=True,
        help='path to tile url list'
    )

    parser.add_argument(
        '--grid_size',
        type=int,
        required=True,
        help='desired size of output tiles'
    )

    parser.add_argument(
        '--buffer',
        type=int,
        required=True,
        help='desired size of output tiles'
    )

    parser.add_argument(
        '--bucket',
        type=str,
        required=True,
        help='gc bucket name'
    )

    return parser.parse_args()


def make_test_args():
    args = argparse.Namespace()
    args.tindex = Path('/home/michael/Work/dual-sensor-tree-mortality/input_files/lpc_index.gpkg')
    args.urls=Path('/home/michael/Work/dual-sensor-tree-mortality/input_files/chd_urls.txt')
    args.grid_size = 1000
    args.buffer = 30
    return args


@np.vectorize(excluded=['url_list'])
def get_tile_urls(tiles, url_list):
    urls = []
    for t in tiles:
        try:
            urls.append([url for url in url_list if t in url][0])
            return urls
        except:
            return []


def download_file(url, save_dir, tries=5):
    '''
    Download a file from a URL to a specified directory and
    return the file path and log message.'''
    local_filename = os.path.join(save_dir, url.split('/')[-1])
    # no neet to download if it already there
    if os.path.exists(local_filename):
        message = f'{datetime.now()} - ALREADY EXISTS: {local_filename} has already been downloaded\n'
        return None, message
    attempt = 1
    while attempt <= tries:
        try:
            with requests.get(url, stream=True) as r:
                r.raise_for_status()
                with open(local_filename, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
            # return path and message
            message = f'{datetime.now()} - SUCCESS: Downloaded {local_filename}\n'
            return local_filename, message  
        except requests.exceptions.RequestException as e:
            if attempt == tries:
                message = f'{datetime.now()} - FAILURE: Could not download {url}. Error: {e}\n'
                # return None and error message
                return None, message
            else:
                attempt = attempt + 1  


def make_grid(args):
    # read and preprocess tindex
    df = gpd.read_file(args.tindex).to_crs(6339)
    df = df[['location', 'geometry']]

    # read url list
    url_list = pd.read_csv(
        args.urls,
        header=None,
        names=['url']
    ).url.to_list()

    # put urls into the tile lists
    df['tile'] = df['location'].apply(lambda x: [
        url
        for url
        in url_list
        if Path(x).stem
        in url
        ] + ['None']).apply(lambda y: y[0])
    
    # make goejson showing tiles missing urls
    missing_geojson = META_DIR / 'missing_tiles.geojson'
    if missing_geojson.exists():
        missing_geojson.unlink()
    df[df.tile == 'None'].to_file(
        missing_geojson,
        driver='GeoJson')
    
    # drop tilels missing urls
    df = df[df.tile != 'None']
    
    # add centroid column
    df['centroid'] = df.geometry.centroid
    df = df[['tile', 'centroid', 'geometry']]

    # get extent and subdivide into overlapping tiles
    xmin, ymin, xmax, ymax = df.total_bounds
    cols = list(np.arange(xmin, xmax + args.grid_size, args.grid_size))
    rows = list(np.arange(ymin, ymax + args.grid_size, args.grid_size))

    polygons = []
    for x in cols[:-1]:
        for y in rows[:-1]:
            polygons.append(
                Polygon([
                    (x, y),
                    (x + args.grid_size, y),
                    (x + args.grid_size, y + args.grid_size), (x, y + args.grid_size)
                ])
            )

    # make into gdf and set crs to match original tile file
    grid = gpd.GeoDataFrame({'geometry':polygons})
    grid['geometry'] = grid.geometry.buffer(args.buffer, cap_style='square')
    grid = grid.set_crs(df.crs)

    # add centroid as column and make a tile name based on the new centroid
    roids = grid.centroid.get_coordinates()
    roids['x_'] = roids.x.apply(lambda x: str(round(x)))
    roids['y_'] = roids.y.apply(lambda y: str(round(y)))
    grid['name'] = [f'{t.x_}_{t.y_}' for _, t in roids.iterrows()]

    # make a column containing the centroid as string for output file
    roids['x'] = roids.x.apply(lambda x: str(x))
    roids['y'] = roids.y.apply(lambda y: str(y))
    grid['centroid'] = [f'{t.x} {t.y}' for _, t in roids.iterrows()]

    # find the tiles to be merged
    def find_tiles(row, df):
        return df[df.intersects(row.geometry)].tile.to_list()
        
    results = Parallel(n_jobs=-1)(delayed(find_tiles)(
        row,
        df)
    for _, row in tqdm(grid.iterrows(), total=len(df)))

    # put tiles into df
    grid['tiles'] = results
    # drop rows with empty lists
    grid = grid[grid.tiles.str.len() > 0]
    

    # write tile bounds as string for PDAL pipline
    grid['str_bounds'] = [
        f'([{row.minx},{row.maxx}],[{row.miny},{row.maxy}])'
        for _, row
        in grid.bounds.iterrows()
        ]

    # write gpkg
    grid.to_file(
        grid_gpkg,
        driver='GPKG')

    return grid


def download_blob(bucket_name, src, dst):
    '''Downloads a file (src) from the bucket to dst.'''
    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)

    # Construct a client side representation of a blob.
    blob = bucket.blob(src)
    blob.download_to_filename(dst)

    if log:
        return f'{datetime.now()} - SUCCESS: Downloaded {bucket_name}:{src} to {dst}.'


def upload_file(bucket_name, src, dst, log=False):
    '''Uploads a file (src) to the bucket.'''
    
    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(dst)

    # generation_match precondition=0. for dst object that does not yet exist
    
    generation_match_precondition = 0

    blob.upload_from_filename(
        src,
        if_generation_match=generation_match_precondition
        )

    if log:
        return f'{datetime.now()} - SUCCESS: {src} uploaded to {bucket_name}:{dst}.'


#%%

if __name__ == '__main__':

    # get args
    #args = make_test_args()
    args = parse_arguments()

    # make some paths
    SAVE_DIR = Path(args.tindex).parent / 'output'
    MERGED_DIR = SAVE_DIR / 'merged'
    TMP_DIR = Path(args.tindex).parent / 'tmp'
    META_DIR = SAVE_DIR / 'metadata'
    LOG_DIR = SAVE_DIR / 'logs'
    JSON_DIR = LOG_DIR / 'pdal_jsons'
    LOG_FILE = LOG_DIR / 'download_log.txt'

    # make sure dirs exists
    for p in [SAVE_DIR, MERGED_DIR, TMP_DIR, LOG_DIR, JSON_DIR, META_DIR]:
        if not os.path.exists(p):
            os.makedirs(p)

    # init the log messages
    log_messages = ['-----------------------------------\n'
    ]

    # read or make grid
    grid_gpkg = META_DIR / 'tile_grid.gpkg'
    if (grid_gpkg).exists():
        with LOG_FILE.open('a') as log:
            log_messages.append(f'------------{datetime.now()}--------\n',)
            log_messages.append('----------RESTART------------------\n')
            log.writelines(log_messages)
        print('Grid exists.')
        grid = gpd.read_file(grid_gpkg)
        # fix double quoted list
        grid['tiles'] = grid.tiles.apply(literal_eval)
    else:
        with LOG_FILE.open('a') as log:
            log_messages.append(f'------------{datetime.now()}--------\n',)
            log_messages.append('----------START FROM SCRATCH-------\n')
            log.writelines(log_messages)
        print('Making grid. This will take a while.')
        grid = make_grid(args)

    # see if files have already been merged
    already_merged = list(MERGED_DIR.glob('*.laz'))

    for _, row in tqdm(grid.iterrows(), total=len(grid)):
        # download files in parallel, collect paths and messages
        print('Downloading (or skipping if its already there).')
        results = Parallel(n_jobs=-1)(
            delayed(download_file)(url, TMP_DIR) 
            for url
            in row.tiles)

        # separate results into file paths and messages
        file_paths = []
        tile = row['name'] # row has a method colled name :(
        
        if f'{tile}.laz' in already_merged:
            print(f'{tile}.laz is already there. Skipping')
            continue

        log_messages = [f'tile: {tile}']
        for result in results:
            file_path, message = result
            if file_path is not None:
                file_paths.append(file_path)
            log_messages.append(message)


        # write log messages to logfile
        with LOG_FILE.open('a') as log:
            log.writelines(log_messages)

        # create PDAL pipeline
        tile_path = MERGED_DIR / f'{tile}.laz'
        pipe = []
        pipe.extend(file_paths)
        pipe.append({
            'type': 'filters.crop',
            'bounds': row.str_bounds
        })
        pipe.append({
            'type': 'writers.las',
            'filename': str(tile_path)
        })

        # write pipeline to JSON file
        json_path = JSON_DIR / f'{tile}.json'
        with open(json_path, 'w') as json_file:
            json.dump(pipe, json_file, indent=2)

        # execute pipeline
        try:
            pipeline = pdal.Pipeline(json.dumps(pipe))
            count = pipeline.execute()
            metadata = pipeline.metadata
            pdal_log = pipeline.log

                # write metadata and pdal log
            meta_path = META_DIR / f'{tile}_metadata.json'
            with open(meta_path, 'w') as json_file:
                json.dump(metadata, json_file, indent=2)

            pdal_log_path = LOG_DIR / f'{tile}_pdal_pipeline_log.txt'
            with pdal_log_path.open('a') as log:
                log.writelines(pdal_log)

            # and write entry to main log
            log_messages = [
                '\n',
                f'{datetime.now()} - MERGED: {tile} merged with {count} points.\n',
                f'metadata saved to {str(meta_path)}\n',
                f'pdal log saved to {pdal_log_path}\n'
            ]
            
            with LOG_FILE.open('a') as log:
                log.writelines(log_messages)

        except Exception as e:
            log_messages = [
            '\n',
            f'{datetime.now()} - PDAL PIPELINE FAILED: {tile} \n',
            f'{e}\n',
            '\n'
        ]
        
        with LOG_FILE.open('a') as log:
            log.writelines(log_messages)


       

    

   