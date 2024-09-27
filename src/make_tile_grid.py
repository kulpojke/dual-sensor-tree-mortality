#!bin/python3

# usage:
# python3 make_tile_grid.py \
# --tindex=/home/michael/Work/dual-sensor-tree-mortality/input_files/lpc_index.gpkg \
# --grid_size=2000 \
# --buffer=50 \
# --geojson = /home/michael/Work/dual-sensor-tree-mortality/input_files/retile.geojson \
# --txt=/home/michael/Work/dual-sensor-tree-mortality/input_files/retile.txt


import geopandas as gpd
from pathlib import Path
import argparse
from shapely import Polygon
from joblib import Parallel, delayed
from tqdm import tqdm
import numpy as np


def parse_arguments():
    '''parses the arguments, returns args'''

    parser.add_argument(
        '--tindex',
        type=str,
        required=True,
        help='path to tile index file'
    )

    parser.add_argument(
        '--txt,
        type=str,
        required=True,
        help='path to output retile list'
    )

    parser.add_argument(
        '--geojson',
        type=str,
        required=False,
        default=None,
        help='optional - path to output geojson of retile'
    )

    parser.add_argument(
        '--grid_size,
        type=int,
        required=True,
        help='desired size of output tiles'
    )

    parser.add_argument(
        '--buffer,
        type=int,
        required=True,
        help='desired size of output tiles'
    )

    parser = argparse.ArgumentParser()


def make_test_args():
    args = argparse.Namespace()
    args.tindex = Path('/home/michael/Work/dual-sensor-tree-mortality/input_files/lpc_index.gpkg')
    args.grid_size = 2000
    args.buffer = 50
    args.geojson = Path('/home/michael/Work/dual-sensor-tree-mortality/input_files/retile.geojson')
    args.txt = Path('/home/michael/Work/dual-sensor-tree-mortality/input_files/retile.txt')
    return args

if __name__ == '__main__':

    # get args
    #args = make_test_args()
    args = parse_arguments()

    # read and preprocess tindex
    df = gpd.read_file(
        Path('/home/michael/Work/dual-sensor-tree-mortality/input_files/CarrHirzDelta1_index.gpkg')
    ).to_crs(26910)
    df = df[['location', 'geometry']]
    df['tile'] = df['location'].apply(lambda x: Path(x).stem)
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

    grid = gpd.GeoDataFrame({'geometry':polygons})
    grid['geometry'] = grid.geometry.buffer(args.buffer)
    grid = grid.set_crs(df.crs)

    name = grid.centroid.get_coordinates()
    name['x'] = name.x.apply(lambda x: str(round(x)))
    name['y'] = name.y.apply(lambda y: str(round(y)))
    name = [f'{t.x}_{t.y}' for _, t in name.iterrows()]
    grid['name'] = name
    #
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
    # make lists into string
    grid['tiles'] = grid.tiles.apply(lambda x: ', '.join(x))

    if args.geojson:
        grid.to_file(args.geojson, driver='GeoJSON')

    grid[['name', 'tiles']].to_csv(args.txt, header=None, index=False)
