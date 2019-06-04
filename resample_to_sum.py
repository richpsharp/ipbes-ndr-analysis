"""Resample inputs."""
import argparse
import os
import sys
import glob
import logging

import numpy
from osgeo import gdal
import pygeoprocessing

gdal.SetCacheMax(2**30)


DEFAULT_GTIFF_CREATION_OPTIONS = (
    'TILED=YES', 'BIGTIFF=YES', 'COMPRESS=DEFLATE',
    'BLOCKXSIZE=256', 'BLOCKYSIZE=256')

logging.basicConfig(
    level=logging.DEBUG,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(funcName)s:%(lineno)d] %(message)s'),
    stream=sys.stdout)
LOGGER = logging.getLogger(__name__)


def main():
    """Entry point."""

    parser = argparse.ArgumentParser(description='Coarsen rasters.')
    parser.add_argument('integers', metavar='N', type=int, nargs='+',
                        help='an integer for the accumulator')
    parser.add_argument('--sum', dest='accumulate', action='store_const',
                        const=sum, default=max,
                        help='sum the integers (default: find the max)')

    factor_list = [10]
    for pattern in sys.argv[1::]:
        for path in glob.glob(pattern):
            for factor in factor_list:
                target_path = '%s_%d.tif' % (
                    os.path.splitext(os.path.basename(path))[0], factor)
                raster_info = pygeoprocessing.get_raster_info(path)
                print(target_path)

                # make a new raster that's / factor as big
                raster_x_size, raster_y_size = raster_info['raster_size']
                nodata = raster_info['nodata'][0]
                n_cols = int(raster_x_size / factor)
                n_rows = int(raster_y_size / factor)
                new_gt = list(raster_info['geotransform'])
                new_gt[1] *= factor
                new_gt[5] *= factor
                driver = gdal.GetDriverByName('GTiff')

                target_raster = driver.Create(
                    target_path, n_cols, n_rows, 1, raster_info['datatype'],
                    options=DEFAULT_GTIFF_CREATION_OPTIONS)
                target_raster.SetProjection(raster_info['projection'])
                target_raster.SetGeoTransform(new_gt)

                target_band = target_raster.GetRasterBand(1)
                target_band.SetNoDataValue(nodata)

                base_raster = gdal.OpenEx(path, gdal.OF_RASTER)
                base_band = base_raster.GetRasterBand(1)

                for xi in range(n_cols):
                    LOGGER.debug('%d of %d', xi, n_cols)
                    for yi in range(n_rows):
                        win_xsize = factor
                        win_ysize = factor
                        if xi*factor + win_xsize > raster_x_size:
                            win_xsize = raster_x_size - xi*factor
                        if yi*factor + win_ysize > raster_y_size:
                            win_ysize = raster_y_size - yi*factor

                        base_array = base_band.ReadAsArray(
                            xoff=xi*factor,
                            yoff=yi*factor,
                            win_xsize=win_xsize,
                            win_ysize=win_ysize)

                        valid_mask = ~numpy.isclose(base_array, nodata)
                        if valid_mask.any():
                            target_band.WriteArray(
                                numpy.array([[numpy.sum(base_array[valid_mask])]]),
                                xoff=xi, yoff=yi)
                        else:
                            target_band.WriteArray(
                                numpy.array([[nodata]]), xoff=xi, yoff=yi)


if __name__ == '__main__':
    main()
