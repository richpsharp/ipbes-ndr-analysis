"""Script to mosaic NDR results into single rasters."""
import sys
import logging
import os

from osgeo import gdal
from osgeo import osr
import pygeoprocessing
import taskgraph

WORKSPACE_DIR = 'mosaic_workspace'
N_WORKERS = -1
TASKGRAPH_UPDATE_INTERVAL = 5.0

logging.basicConfig(
    level=logging.DEBUG,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(funcName)s:%(lineno)d] %(message)s'),
    stream=sys.stdout)
LOGGER = logging.getLogger(__name__)

NDR_DIRECTORY = os.path.join(
    'ipbes_ndr_workspace', 'watershed_processing')

MOSAIC_CELL_SIZE = 1.0
RASTER_SUFFIXES_TO_AGGREGATE = (
    'ssp1_n_export.tif',
    'ssp3_n_export.tif',
    'ssp5_n_export.tif',
    '1850_n_export.tif',
    '1900_n_export.tif',
    '1910_n_export.tif',
    '1945_n_export.tif',
    '1980_n_export.tif',
    'isimip_2015_n_export.tif',
    'worldclim_2015_n_export.tif',
    'esa_2015_n_export.tif',
    'ssp1_modified_load.tif',
    'ssp3_modified_load.tif',
    'ssp5_modified_load.tif',
    '1850_modified_load.tif',
    '1900_modified_load.tif',
    '1910_modified_load.tif',
    '1945_modified_load.tif',
    '1980_modified_load.tif',
    'isimip_2015_modified_load.tif',
    'worldclim_2015_modified_load.tif',
    'esa_2015_modified_load.tif',
    '2015_rural_total_pop_aligned.tif',
    'ssp1_rural_total_pop_aligned.tif',
    'ssp3_rural_total_pop_aligned.tif',
    'ssp5_rural_total_pop_aligned.tif',
    '1850_ag_load_aligned.tif',
    '1900_ag_load_aligned.tif',
    '1920_ag_load_aligned.tif',
    '1945_ag_load_aligned.tif',
    '1980_ag_load_aligned.tif',
    '2015_ag_load_aligned.tif',
    'ssp1_2050_ag_load_aligned.tif',
    'ssp3_2050_ag_load_aligned.tif',
    'ssp5_2050_ag_load_aligned.tif')


def main():
    """Entry point."""
    task_graph = taskgraph.TaskGraph(
        WORKSPACE_DIR, N_WORKERS, TASKGRAPH_UPDATE_INTERVAL)

    try:
        os.makedirs(WORKSPACE_DIR)
    except OSError:
        pass

    global_raster_task_path_map = {}
    for (dirpath, dirnames, filenames) in os.walk(NDR_DIRECTORY):
        if dirnames:
            continue
        for raster_suffix in RASTER_SUFFIXES_TO_AGGREGATE:
            try:
                matching_path = next(iter(
                    (os.path.join(dirpath, file_path) for file_path in filenames
                     if file_path.endswith(raster_suffix))))
            except StopIteration:
                raise ValueError(
                    "Expected to find %s in %s but not found" % (
                        raster_suffix, dirpath))

            base_raster_info = pygeoprocessing.get_raster_info(matching_path)
            target_raster_path = os.path.join(WORKSPACE_DIR, raster_suffix)
            target_token_complete_path = f'''{
                os.path.splitext(target_raster_path)[0]}_{
                    MOSAIC_CELL_SIZE}.TOKEN'''
            LOGGER.debug(target_raster_path)
            make_empty_raster_task = task_graph.add_task(
                func=make_empty_wgs84_raster,
                args=(
                    MOSAIC_CELL_SIZE, base_raster_info['nodata'][0],
                    base_raster_info['datatype'], target_raster_path,
                    target_token_complete_path),
                target_path_list=[target_token_complete_path],
                task_name=f'create empty global {raster_suffix}')
            global_raster_task_path_map[raster_suffix] = (
                make_empty_raster_task, target_raster_path)
        LOGGER.info("found all the raster suffixes in %s", dirpath)
        break

    task_graph.close()
    task_graph.join()


def make_empty_wgs84_raster(
        cell_size, nodata_value, target_datatype, target_raster_path,
        target_token_complete_path):
    """Make a big empty raster in WGS84 projection.

    Parameters:
        cell_size (float): this is the desired cell size in WSG84 degree
            units.
        nodata_value (float): desired nodata avlue of target raster
        target_datatype (gdal enumerated type): desired target datatype.
        target_raster_path (str): this is the target raster that will cover
            [-180, 180), [90, -90) with cell size units with y direction being
            negative.

    Returns:
        None.

    """
    gtiff_driver = gdal.GetDriverByName('GTiff')
    try:
        os.makedirs(os.path.dirname(target_raster_path))
    except OSError:
        pass

    n_cols = int(360.0 / cell_size)
    n_rows = int(180.0 / cell_size)

    geotransform = (-180.0, cell_size, 0.0, 90.0, 0, -cell_size)

    wgs84_srs = osr.SpatialReference()
    wgs84_srs.ImportFromEPSG(4326)

    target_raster = gtiff_driver.Create(
        target_raster_path, n_cols, n_rows, 1, target_datatype,
        options=(
            'TILED=YES', 'BIGTIFF=YES', 'COMPRESS=DEFLATE',
            'BLOCKXSIZE=256', 'BLOCKYSIZE=256'))
    target_raster.SetProjection(wgs84_srs.ExportToWkt())
    target_raster.SetGeoTransform(geotransform)
    target_band = target_raster.GetRasterBand(1)
    target_band.SetNoDataValue(nodata_value)
    target_band.Fill(nodata_value)
    target_band.FlushCache()
    target_band = None
    target_raster = None

    target_raster = gdal.OpenEx(target_raster_path, gdal.OF_RASTER)
    if target_raster:
        with open(target_token_complete_path, 'w') as target_token_file:
            target_token_file.write('complete!')


if __name__ == '__main__':
    main()
