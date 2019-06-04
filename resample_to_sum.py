"""Resample inputs."""
import os
import sys
import glob

from osgeo import gdal
import pygeoprocessing


DEFAULT_GTIFF_CREATION_OPTIONS = (
    'TILED=YES', 'BIGTIFF=YES', 'COMPRESS=DEFLATE',
    'BLOCKXSIZE=256', 'BLOCKYSIZE=256')


def main():
    """Entry point."""
    factor_list = [100]
    for pattern in sys.argv[1::]:
        for path in glob.glob(pattern):
            for factor in factor_list:
                target_path = '%s_%d.tif' % (
                    os.path.splitext(os.path.basename(path))[0], factor)
                raster_info = pygeoprocessing.get_raster_info(path)
                print(target_path)

                # make a new raster that's / factor as big
                n_cols = raster_info['raster_size'][0] / factor
                n_rows = raster_info['raster_size'][1] / factor

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

                for xi in range(n_cols):
                    for yi in range(n_rows):



if __name__ == '__main__':
    main()
