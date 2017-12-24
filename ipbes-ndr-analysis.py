"""Script to manage NDR runs for IPBES project."""
import logging
import os
import pickle

import taskgraph
import rtree.index
import glob
from osgeo import ogr
import pygeoprocessing

N_CPUS = 2

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

WATERSHED_PATH_LIST = glob.glob(
    os.path.join(BASE_DROPBOX_DIR, 'ipbes-data',
        'watersheds_globe_HydroSHEDS_15arcseconds', '*.shp'))

TARGET_WORKSPACE = 'ndr_workspace'
TASKGRAPH_DIR = os.path.join(TARGET_WORKSPACE, 'taskgraph_cache')

RTREE_PATH = 'dem_rtree.dat'


def main():
    """Entry point."""
    logging.basicConfig(
        format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
        level=logging.WARN, datefmt='%m/%d/%Y %H:%M:%S ')

    if not os.path.exists(TARGET_WORKSPACE):
        os.makedirs(TARGET_WORKSPACE)

    task_graph = taskgraph.TaskGraph(
        TASKGRAPH_DIR, N_CPUS)

    dem_rtree_path = os.path.join(TARGET_WORKSPACE, RTREE_PATH)

    dem_path_list = glob.glob(os.path.join(
        BASE_DROPBOX_DIR, 'dataplatform',
        'dem_globe_CGIAR_STRMv41_3arcseconds', '*.tif'))

    dem_path_index_map = {}
    build_dem_rtree_task = task_graph.add_task(
        func=build_dem_rtree,
        args=(dem_path_list, dem_path_index_map, dem_rtree_path),
        target_path_list=[dem_rtree_path])

    for watershed_path in WATERSHED_PATH_LIST:
        watershed_vector = ogr.Open(watershed_path)
        watershed_layer = watershed_vector.GetLayer()
        watershed_basename = os.path.splitext(
            os.path.basename(watershed_path))[0]
        for watershed_id in xrange(watershed_layer.GetFeatureCount()):

            target_dem_path = os.path.join(
                "ws_%s_%d_working_dir" % (watershed_basename, watershed_id),
                'ws_%d_dem.tif' % watershed_id)
            task_graph.add_task(
                func=merge_watershed_dems,
                args=(
                    watershed_path, watershed_id, dem_rtree_path,
                    dem_path_index_map, target_dem_path),
                dependent_task_list=[build_dem_rtree_task],
                task_name='merge_watershed_%d_dems' % watershed_id)

    task_graph.close()
    task_graph.join()


def merge_watershed_dems(
        watershed_path, watershed_id, dem_rtree_path, dem_path_index_map,
        target_dem_path):
    """Find DEMs that overlap the given watershed polyon by id.

    Parameters:
        watershed_path (string):
        watershed_id (int): feature number to index and overlap.
        dem_rtree_path (string): path to a pickled rtree that maps
            bounding boxes to dem ids.
        dem_path_index_map (dict): maps DEM ids from the rtree to filepaths.
        dem_path_id_map (dict): map dem id to dem paths.
        target_dem_path (string): path to file that's created by
            mosaicing all the overlapping dems together, suitable for
            routing in the given watershed.

    Returns:
        None.
    """
    watershed_vector = ogr.Open(watershed_path)
    watershed_layer = watershed_vector.GetLayer()
    watershed_feature = watershed_layer.GetFeature(watershed_id)
    watershed_geometry = watershed_feature.GetGeometryRef()
    watershed_bb = [
        watershed_geometry.GetEnvelope()[i] for i in [0, 2, 1, 3]]
    watershed_geometry = None
    watershed_feature = None
    watershed_vector = None

    with open(dem_rtree_path, 'rb') as dem_rtree_file:
        dem_rtree = pickle.load(dem_rtree_file)

    overlapping_dem_list = list(dem_rtree.intersection(watershed_bb))

    if len(overlapping_dem_list) > 0:
        LOGGER.debug(
            watershed_id, [
                dem_path_index_map[i] for i in overlapping_dem_list])
        workspace_dir = os.path.dirname(target_dem_path)
        if not os.path.exists(workspace_dir):
            os.makedirs(workspace_dir)
        pygeoprocessing.merge_rasters(overlapping_dem_list, target_dem_path)
    else:
        LOGGER.debug(
            "no overlapping dems found for %s wsid %d", watershed_path,
            watershed_id)



def build_dem_rtree(dem_path_list, dem_path_index_map, dem_rtree_path):
    """Build RTree indexed by FID for points in `wwwiii_vector_path`."""
    dem_rtree = rtree.index.Index()
    for dem_id, dem_path in enumerate(dem_path_list):
        raster_info = pygeoprocessing.get_raster_info(dem_path)
        dem_path_index_map[dem_id] = dem_path
        print os.path.basename(dem_path), raster_info['bounding_box']
        dem_rtree.insert(dem_id, raster_info['bounding_box'])
    with open(dem_rtree_path, 'wb') as dem_rtree_file:
        pickle.dump(dem_rtree, dem_rtree_file)


if __name__ == '__main__':
    main()
