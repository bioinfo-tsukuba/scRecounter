#!/usr/bin/env python
# import
## batteries
from __future__ import print_function
import os
import re
import sys
import argparse
import logging
import json


# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

# argparse
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'Load JSON file and print parameter'
epi = """DESCRIPTION:
Load JSON file and print parameter
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
parser.add_argument('json_file', type=str, help='JSON file')
parser.add_argument('--params', type=str, nargs='+', default=[],
                    help='Parameters to export as environment variables',)


# functions
def main(args):
    # load json
    with open(args.json_file, 'r') as f:
        data = json.load(f)

    # export as environment variables
    for key, val in data.items():
        if key in args.params:
            logging.info(f'Exporting: {key}={val}')
            print(f'export {key}={val}')

## script main
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)