"""Script to manage NDR runs for IPBES project."""
import logging
import os

import taskgraph
import rtree.index
import glob
from osgeo import ogr
import pygeoprocessing

N_CPUS = -1

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
    'WaterGap3_AllBasins.shp')

TARGET_WORKSPACE = 'ndr_workspace'
TASKGRAPH_DIR = os.path.join(TARGET_WORKSPACE, 'taskgraph_cache')

RTREE_PATH = 'dem_rtree.dat'

DEM_PATH_LIST = glob.glob(os.path.join(
    'D:/', 'dataplatformdata', 'dem_globe_ASTER_1arcsecond', '*.tif'))

def main():

    logging.basicConfig(
        format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
        level=logging.WARN, datefmt='%m/%d/%Y %H:%M:%S ')

    if not os.path.exists(TARGET_WORKSPACE):
        os.makedirs(TARGET_WORKSPACE)

    task_graph = taskgraph.TaskGraph(
        TASKGRAPH_DIR, N_CPUS)

    dem_rtree_path = os.path.join(
        TARGET_WORKSPACE, RTREE_PATH)

    dem_path_index_map = {}
    build_dem_rtree(DEM_PATH_LIST, dem_path_index_map, dem_rtree_path)
    dem_rtree = rtree.index.Index(dem_rtree_path)

    watershed_vector = ogr.Open(WATERSHED_PATH)
    watershed_layer = watershed_vector.GetLayer()
    for watershed_id, watershed_feature in enumerate(watershed_layer):
        watershed_geometry = watershed_feature.GetGeometryRef()
        watershed_bb = [
            watershed_geometry.GetEnvelope()[i] for i in [0, 2, 1, 3]]
        overlapping_dem_list = list(dem_rtree.intersection(watershed_bb))
        print watershed_id, [
            dem_path_index_map[i] for i in overlapping_dem_list]

        watershed_workspace_dir = os.path.join(
            TARGET_WORKSPACE, 'watershed_%d' % watershed_id)
        if not os.path.exists(watershed_workspace_dir):
            os.makedirs(watershed_workspace_dir)
        watershed_dem_path = os.path.join(
            watershed_workspace_dir, 'watershed_dem.tif')
        pygeoprocessing.merge_rasters(
            overlapping_dem_list, watershed_dem_path)
        break
    # convert form [minx,maxx,miny,maxy] to [minx,miny,maxx,maxy]
    #vector_properties['bounding_box'] = [layer_bb[i] for i in [0, 2, 1, 3]]


def build_dem_rtree(dem_path_list, dem_path_index_map, dem_rtree_path):
    """Build RTree indexed by FID for points in `wwwiii_vector_path`."""
    base_dem_rtree_path = os.path.splitext(dem_rtree_path)[0]
    if os.path.exists(dem_rtree_path):
        return
        for ext in ['.dat', '.idx']:
            os.remove(base_dem_rtree_path+ext)
    dem_rtree = rtree.index.Index(base_dem_rtree_path)

    for dem_id, dem_path in enumerate(DEM_PATH_LIST):
        raster_info = pygeoprocessing.get_raster_info(dem_path)
        dem_path_index_map[dem_id] = dem_path
        print os.path.basename(dem_path), raster_info['bounding_box']
        dem_rtree.insert(dem_id, raster_info['bounding_box'])


if __name__ == '__main__':
    main()
