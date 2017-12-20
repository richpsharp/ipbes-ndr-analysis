"""Script to manage NDR runs for IPBES project."""
import logging
import os

import rtree.index
import glob
from osgeo import ogr
import pygeoprocessing

logging.basicConfig(
    format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
    level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')

LOGGER = logging.getLogger('ipbes_pollination_analysis')

POSSIBLE_DROPBOX_LOCATIONS = [
    r'D:\Dropbox',
    r'C:\Users\Rich\Dropbox',
    r'C:\Users\rpsharp\Dropbox',
    r'E:\Dropbox']

LOGGER.info("checking dropbox locations")
for path in POSSIBLE_DROPBOX_LOCATIONS:
    print path
    if os.path.exists(path):
        BASE_DROPBOX_DIR = path
        break
LOGGER.info("found %s", BASE_DROPBOX_DIR)

WATERSHED_PATH = os.path.join(
    BASE_DROPBOX_DIR, 'ipbes-data', 'WaterDepletionWG3_watersheds',
    'WaterDepletionWG3.shp')


def main():
    dem_index = rtree.index.Index()

    dem_path_list = glob.glob(os.path.join(
        'D:/', 'dataplatformdata', 'dem_globe_ASTER_1arcsecond', '*.tif'))
    dem_path_index_map = {}
    count = 100
    for dem_id, dem_path in enumerate(dem_path_list):
        raster_info = pygeoprocessing.get_raster_info(dem_path)
        dem_path_index_map[dem_id] = dem_path
        print os.path.basename(dem_path), raster_info['bounding_box']
        dem_index.insert(dem_id, raster_info['bounding_box'])

    watershed_vector = ogr.Open(WATERSHED_PATH)
    watershed_layer = watershed_vector.GetLayer()
    for watershed_id, watershed_feature in enumerate(watershed_layer):
        watershed_geometry = watershed_feature.GetGeometryRef()
        watershed_bb = [
            watershed_geometry.GetEnvelope()[i] for i in [0, 2, 1, 3]]
        overlapping_dem_list = list(dem_index.intersection(watershed_bb))
        print watershed_id, [
            dem_path_index_map[i] for i in overlapping_dem_list]
    # convert form [minx,maxx,miny,maxy] to [minx,miny,maxx,maxy]
    #vector_properties['bounding_box'] = [layer_bb[i] for i in [0, 2, 1, 3]]

if __name__ == '__main__':
    main()
