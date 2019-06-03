"""Subtract two command line arguments."""
import itertools
import time
import sys
import logging
import os
import multiprocessing

import numpy
from osgeo import gdal
from osgeo import osr
import pygeoprocessing
import taskgraph

# set a 1GB limit for the cache
gdal.SetCacheMax(2**30)

WORKSPACE_DIR = 'mosaic_workspace'
N_WORKERS = min(8, multiprocessing.cpu_count())
TASKGRAPH_UPDATE_INTERVAL = 5.0

logging.basicConfig(
    level=logging.DEBUG,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(funcName)s:%(lineno)d] %(message)s'),
    stream=sys.stdout)
LOGGER = logging.getLogger(__name__)


def main():
    """Entry point."""
    raster_a_path = sys.argv[1]
    raster_b_path = sys.argv[2]
    target_raster_path = sys.argv[3]

    raster_a_nodata = pygeoprocessing.get_raster_info(raster_a_path)['nodata'][0]
    raster_b_nodata = pygeoprocessing.get_raster_info(raster_b_path)['nodata'][0]

    def subtract(array_a, array_b):
        result = numpy.empty(array_a.shape)
        result[:] = raster_a_nodata
        valid_mask = (
            ~numpy.isclose(array_a, raster_a_nodata) &
            ~numpy.isclose(array_b, raster_b_nodata))
        result[valid_mask] = array_a[valid_mask]-array_b[valid_mask]
        return result

    pygeoprocessing.raster_calculator(
        [(raster_a_path, 1), (raster_b_path, 1)], subtract,
        target_raster_path, gdal.GDT_Float32, raster_a_nodata)


if __name__ == '__main__':
    main()
