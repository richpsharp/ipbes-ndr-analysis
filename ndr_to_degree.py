"""Script to manage NDR runs for IPBES project."""
import heapq
import shutil
import time
import datetime
import threading
import logging
import os
import glob
import math
import sqlite3
import multiprocessing

import taskgraph
import numpy
import pandas
import dill
import rtree.index
from osgeo import ogr
from osgeo import gdal
from osgeo import osr
import pygeoprocessing
import pygeoprocessing.routing


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

WATERSHED_PATH_LIST = glob.glob(
    os.path.join(
        BASE_DROPBOX_DIR, 'ipbes-data',
        'watersheds_globe_HydroSHEDS_15arcseconds', '*.shp'))

DEGREE_GRID_PATH = os.path.join(
    BASE_DROPBOX_DIR, 'ipbes stuff',
    'summary table shapefile', 'degree_basedata', 'grid_1_degree.shp')


TARGET_WORKSPACE = 'ndr_to_degree_workspace'
TASKGRAPH_DIR = os.path.join(TARGET_WORKSPACE, 'taskgraph_cache')

RTREE_PATH = 'watershed_rtree'


def build_watershed_rtree(
        watershed_path_list, watershed_path_index_map_path):
    """Build RTree indexed by FID for points in `wwwiii_vector_path`."""
    LOGGER.info('building rTree %s', watershed_path_index_map_path +'.dat')
    if os.path.exists(watershed_path_index_map_path+'.dat'):
        LOGGER.warn('%s exists so skipping creation.', watershed_path_index_map_path)
        return
    watershed_rtree = rtree.index.Index(watershed_path_index_map_path)
    for global_watershed_path in watershed_path_list:
        print global_watershed_path
        watershed_basename = os.path.splitext(
            os.path.basename(global_watershed_path))[0]
        watershed_vector = gdal.OpenEx(global_watershed_path, gdal.OF_VECTOR)
        watershed_layer = watershed_vector.GetLayer()
        for watershed_id in xrange(watershed_layer.GetFeatureCount()):
            watershed_feature = watershed_layer.GetFeature(watershed_id)
            feature_geom = watershed_feature.GetGeometryRef()
            ws_prefix = 'ws_%s_%d' % (
                watershed_basename, watershed_feature.GetField('BASIN_ID'))

            watershed_area = feature_geom.GetArea()
            if watershed_area < 0.03:
                #  0.04 square degrees is a healthy underapproximation of
                # 100 sq km which is about the minimum watershed size we'd
                # want.
                continue
            x_min, x_max, y_min, y_max = feature_geom.GetEnvelope()

            watershed_rtree.insert(
                0, (x_min, y_min, x_max, y_max), obj=ws_prefix)
            feature_geom = None
            feature_geom = None
            watershed_feature = None
        watershed_layer = None
    watershed_vector = None


def main():
    """Entry point."""
    try:
        os.makedirs(TARGET_WORKSPACE)
    except OSError:
        pass
    watershed_path_index_map_path = os.path.join(
        TARGET_WORKSPACE, 'watershed_rtree')
    print WATERSHED_PATH_LIST
    build_watershed_rtree(
        WATERSHED_PATH_LIST, watershed_path_index_map_path)

    watershed_rtree = rtree.index.Index(watershed_path_index_map_path)

    grid_vector = gdal.OpenEx(DEGREE_GRID_PATH, gdal.OF_VECTOR)
    grid_layer = grid_vector.GetLayer()
    while True:
        grid_feature = grid_layer.GetNextFeature()
        if not grid_feature:
            break
        grid_geometry = grid_feature.GetGeometryRef()
        grid_bounds = grid_geometry.GetEnvelope()
        results = list(watershed_rtree.intersection(
            (grid_bounds[0], grid_bounds[2], grid_bounds[1], grid_bounds[3]),
            objects=True))
        if results:
            for watershed_id in [str(x.object) for x in results]:
                watershed_path = os.path.join(
                    'ndr_workspace', '/'.join(reversed(watershed_id[-4:])),
                    '%s_working_dir' % watershed_id, '%s.shp' % watershed_id)
                print watershed_path
                if os.path.exists(watershed_path):
                    print watershed_path


if __name__ == '__main__':
    main()
