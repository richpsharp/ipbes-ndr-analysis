"""Clip warp test."""
import argparse
import logging
import sys

import pygeoprocessing

logging.basicConfig(
    level=logging.DEBUG,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(funcName)s:%(lineno)d] %(message)s'),
    stream=sys.stdout)
LOGGER = logging.getLogger(__name__)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='clip/warp to given vector')
    parser.add_argument('base_raster_path')
    parser.add_argument('base_vector_path')
    parser.add_argument('target_raster_path')
    args = parser.parse_args()

    base_vector_info = pygeoprocessing.get_vector_info(args.base_vector_path)

    pygeoprocessing.warp_raster(
        args.base_raster_path, (300, -300), args.target_raster_path,
        'near', target_bb=base_vector_info['bounding_box'],
        base_projection_wkt=None,
        target_projection_wkt=base_vector_info['projection_wkt'])
