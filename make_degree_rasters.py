"""Script to manage NDR runs for IPBES project."""
import logging
import os
import glob
import sqlite3
import multiprocessing

import numpy
from osgeo import ogr
from osgeo import gdal
from osgeo import osr
import pygeoprocessing


N_CPUS = -1
DRY_RUN = False
TASKGRAPH_REPORTING_FREQUENCY = 60.0

logging.basicConfig(
    format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
    level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')

LOGGER = logging.getLogger('ndr_to_degree')

POSSIBLE_DROPBOX_LOCATIONS = [
    r'D:\Dropbox',
    r'C:\Users\Rich\Dropbox',
    r'C:\Users\rpsharp\Dropbox',
    r'E:\Dropbox']

LOGGER.info("checking dropbox locations")
for dropbox_path in POSSIBLE_DROPBOX_LOCATIONS:
    print dropbox_path
    if os.path.exists(dropbox_path):
        BASE_DROPBOX_DIR = dropbox_path
        break
LOGGER.info("found %s", BASE_DROPBOX_DIR)

TARGET_WORKSPACE = 'ndr_degree_raster_workspace'


def main():
    """Entry point."""
    try:
        os.makedirs(TARGET_WORKSPACE)
    except OSError:
        pass

    database_path = os.path.join(
        'ndr_to_degree_workspace', 'ndr_to_degree.db')

    # create a database connection
    conn = sqlite3.connect(database_path)
    cursor = conn.cursor()

    for summary_field in [
            'cur_n_export', 'ssp1_n_export', 'ssp3_n_export', 'ssp5_n_export',
            'cur_modified_load', 'ssp1_modified_load', 'ssp3_modified_load',
            'ssp5_modified_load']:
        print 'summary %s' % summary_field
        n_export_degree_path = os.path.join(
            TARGET_WORKSPACE, '%s_degree.tif' % summary_field)
        wgs84_sr = osr.SpatialReference()
        wgs84_sr.ImportFromEPSG(4326)
        driver = gdal.GetDriverByName('GTiff')
        print 'create raster'
        summary_raster = driver.Create(
            n_export_degree_path, 361, 181, 1, gdal.GDT_Float64)
        summary_raster.SetProjection(wgs84_sr.ExportToWkt())
        wgs84_gt = [-180.0, 1.0, 0, 90., 0, -1]
        summary_raster.SetGeoTransform(wgs84_gt)
        summary_band = summary_raster.GetRasterBand(1)
        nodata = -1
        summary_band.SetNoDataValue(nodata)
        summary_band.Fill(nodata)
        base_array = numpy.empty((181, 361), dtype=numpy.float32)
        base_array[:] = nodata
        inv_gt = gdal.InvGeoTransform(wgs84_gt)

        cursor.execute(
            "SELECT SUM(%s), SUM(fraction_covered), GRIDCODE "
            "from nutrient_export group by GRIDCODE" % (
                '%s' % summary_field))

        for result in cursor:
            gridcode = int(result[2])
            fraction_covered = float(result[1])
            export_sum = float(result[0])
            ix = (gridcode-1) % 360
            iy = (gridcode-1) // 360
            if fraction_covered != 0.0:
                base_array[iy, ix] = export_sum / fraction_covered
        summary_band.WriteArray(base_array)
    conn.close()


if __name__ == '__main__':
    main()
