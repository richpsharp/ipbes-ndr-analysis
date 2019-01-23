"""Script to mosaic NDR results into single rasters."""
import sys
import logging
import os

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
    'ssp1_ag_load_aligned.tif',
    'ssp3_ag_load_aligned.tif',
    'ssp5_ag_load_aligned.tif',
    '1850_ag_load_aligned.tif',
    '1900_ag_load_aligned.tif',
    '1910_ag_load_aligned.tif',
    '1945_ag_load_aligned.tif',
    '1980_ag_load_aligned.tif',
    '2015_ag_load_aligned.tif',
    'ssp1_2050_ag_load_aligned.tif',
    'ssp3_2050_ag_load_aligned.tif',
    'ssp5_2050_ag_load_aligned.tif')


def main():
    """Entry point."""
    #task_graph = taskgraph.TaskGraph(
    #    WORKSPACE_DIR, N_WORKERS, TASKGRAPH_UPDATE_INTERVAL)

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
        LOGGER.info("found all the raster suffixes in %d", dirpath)
        break


if __name__ == '__main__':
    main()
