#!bin/python
# This compares the lists of existing files.
# Writes list of files that still need processing.
# Called by start_crown_gc.sh

import argparse


def parse_arguments():
    '''parses the arguments, returns args'''
    # init parser
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--laz_list',
        type=str,
        required=True,
        help='list of laz files in dir'
    )

    parser.add_argument(
        '--ttops_list',
        type=str,
        required=True,
        help='list of ttops files in dir'
    )


    parser.add_argument(
        '--crowns_list',
        type=str,
        required=True,
        help='list of crowns files in dir'
    )

    return parser.parse_args()


if __name__ == '__main__':
    # get args
    args = parse_arguments()

    # read files to lists
    with open(args.laz_list) as src:
        laz = src.readlines()

    with open(args.ttops_list) as src:
        ttops = src.readlines()

    with open(args.crowns_list) as src:
        crowns = src.readlines()

    # ttops and crowns will most likely same len, if not use shorter list
    if len(crowns) > len(ttops):
        done_list = [f.split('_')[0] for f in ttops]
    else:
        done_list = [f.split('_')[0] for f in crowns]


    if len(done_list) > 0:
        todo = [l for l in laz for f in done_list if f not in l]
    else:
        todo = done_list

    with open('todo.list', 'w') as dst:
        dst.writelines(line + '\n' for line in todo)



