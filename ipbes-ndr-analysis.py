"""Script to manage NDR runs for IPBES project."""
import logging
import os
import glob
import math

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
NODATA = -1

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
    os.path.join(
        BASE_DROPBOX_DIR, 'ipbes-data',
        'watersheds_globe_HydroSHEDS_15arcseconds', '*.shp'))

TARGET_WORKSPACE = 'ndr_workspace'
TASKGRAPH_DIR = os.path.join(TARGET_WORKSPACE, 'taskgraph_cache')

RTREE_PATH = 'dem_rtree'

def length_of_degree(lat, lng):
    """Calcualte the length of a degree in meters."""
    m1 = 111132.92
    m2 = -559.82
    m3 = 1.175
    m4 = -0.0023
    p1 = 111412.84
    p2 = -93.5
    p3 = 0.118

    lat = -19.7860856226 * math.pi / 180

    latlen = (
        m1 + (m2 * math.cos(2 * lat)) + (m3 * math.cos(4 * lat)) + (m4 * math.cos(6 * lat)))
    longlen = abs(
        (p1 * math.cos(lat)) + (p2 * math.cos(3 * lat)) + (p3 * math.cos(5 * lat)))

    return max(latlen, longlen)

def mult_arrays(*array_list):
    """Multiply arrays in array list but block out stacks with NODATA."""
    stack = numpy.stack(array_list)
    valid_mask = (
        numpy.bitwise_and.reduce(stack != NODATA, axis=0))
    n_valid = numpy.count_nonzero(valid_mask)
    broadcast_valid_mask = numpy.broadcast_to(valid_mask, stack.shape)
    valid_stack = stack[broadcast_valid_mask].reshape(
        len(array_list), n_valid)
    result = numpy.empty(array_list[0].shape, dtype=numpy.float32)
    result[:] = NODATA
    result[valid_mask] = numpy.prod(valid_stack, axis=0)
    return result


def main():
    """Entry point."""
    logging.basicConfig(
        format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
        level=logging.WARN, datefmt='%m/%d/%Y %H:%M:%S ')

    if not os.path.exists(TARGET_WORKSPACE):
        os.makedirs(TARGET_WORKSPACE)

    # load biophysical table first, it's used a lot below
    represenative_ndr_biophysical_table_path = os.path.join(
        BASE_DROPBOX_DIR, 'ipbes-data', 'NDR_representative_table.csv')
    biophysical_table = pandas.read_csv(
        represenative_ndr_biophysical_table_path)
    # clean up biophysical table
    biophysical_table = biophysical_table.fillna(0)
    biophysical_table[biophysical_table['load_n'] == 'use raster'] = 0.0
    biophysical_table['load_n'] = biophysical_table['load_n'].apply(
        pandas.to_numeric)

    task_graph = taskgraph.TaskGraph(
        TASKGRAPH_DIR, N_CPUS)

    dem_rtree_path = os.path.join(TARGET_WORKSPACE, RTREE_PATH)

    dem_path_list = glob.glob(os.path.join(
        BASE_DROPBOX_DIR, 'dataplatform',
        'dem_globe_CGIAR_STRMv41_3arcseconds', '*.tif'))

    dem_pixel_size = pygeoprocessing.get_raster_info(
        dem_path_list[0])['pixel_size']

    dem_path_index_map_path = os.path.join(
        TARGET_WORKSPACE, 'dem_path_index_map.dat')
    build_dem_rtree_task = task_graph.add_task(
        func=build_dem_rtree,
        args=(dem_path_list, dem_path_index_map_path, dem_rtree_path),
        target_path_list=[
            dem_rtree_path+'.dat',  # rtree adds a ".dat" file
            dem_path_index_map_path],
        task_name='build_dem_rtree')

    watershed_path = r"C:\Users\Rich\Dropbox\ipbes-data\watersheds_globe_HydroSHEDS_15arcseconds\af_bas_15s_beta.shp"
    watershed_vector = ogr.Open(watershed_path)
    watershed_layer = watershed_vector.GetLayer()

    watershed_id = 85668
    watershed_basename = os.path.splitext(os.path.basename(watershed_path))[0]
    ws_prefix = 'ws_%s_%d' % (watershed_basename, watershed_id)
    ws_working_dir = os.path.join(
        TARGET_WORKSPACE, "%s_working_dir" % ws_prefix)
    watershed_dem_path = os.path.join(
        ws_working_dir, 'ws_%s_dem.tif' % ws_prefix)
    merge_watershed_dems_task = task_graph.add_task(
        func=merge_watershed_dems,
        args=(
            watershed_path, watershed_id, dem_rtree_path,
            dem_path_index_map_path, watershed_dem_path),
        target_path_list=[watershed_dem_path],
        dependent_task_list=[build_dem_rtree_task],
        task_name='merge_watershed_dems_%s' % ws_prefix)

    watershed_feature = watershed_layer.GetFeature(watershed_id)
    feature_geom = watershed_feature.GetGeometryRef()
    feature_centroid = feature_geom.Centroid()

    utm_code = (math.floor((feature_centroid.GetX() + 180)/6) % 60) + 1
    lat_code = 6 if feature_centroid.GetY() > 0 else 7
    epsg_code = int('32%d%02d' % (lat_code, utm_code))
    epsg_srs = osr.SpatialReference()
    epsg_srs.ImportFromEPSG(epsg_code)
    utm_pixel_size = abs(dem_pixel_size[0]) * length_of_degree(
        feature_centroid.GetY(), feature_centroid.GetX())

    feature_centroid = None
    feature_geom = None
    watershed_feature = None

    # clip precip raster to watershed bb
    precip_raster_path = os.path.join(
        BASE_DROPBOX_DIR, 'ipbes-data',
        'precip_globe_WorldClim_30arcseconds.tif')

    landcover_raster_path = os.path.join(
        BASE_DROPBOX_DIR, 'ipbes-data',
        'GLOBIO4_landuse_10sec_tifs_20171207_Idiv', 'Current2015',
        'Globio4_landuse_10sec_2015_cropint.tif')

    base_raster_path_list = [
        watershed_dem_path, precip_raster_path, landcover_raster_path]

    aligned_path_list = [
        os.path.join(
            ws_working_dir, '%s_%s_aligned.tif' % (
                ws_prefix, os.path.splitext(os.path.basename(x))[0]))
        for x in base_raster_path_list]

    # clip dem, precip, & landcover to size of DEM? use 'mode'
    align_resize_task = task_graph.add_task(
        func=pygeoprocessing.align_and_resize_raster_stack,
        args=(
            base_raster_path_list, aligned_path_list,
            ['nearest', 'nearest', 'mode'], dem_pixel_size, 'intersection'),
        target_path_list=aligned_path_list,
        dependent_task_list=[merge_watershed_dems_task],
        task_name='align_resize_task_%s' % ws_prefix)

    utm_dem_path = os.path.join(
        ws_working_dir, '%s_%s_dem.tif' % (ws_prefix, epsg_code))
    utm_precip_path = os.path.join(
        ws_working_dir, '%s_%s_precip.tif' % (ws_prefix, epsg_code))
    utm_landcover_path = os.path.join(
        ws_working_dir, '%s_%s_landcover.tif' % (ws_prefix, epsg_code))

    path_task_id_map = {}
    for raster_id, base_path, target_path in [
            ('dem', aligned_path_list[0], utm_dem_path),
            ('precip', aligned_path_list[1], utm_precip_path),
            ('landcover', aligned_path_list[2], utm_landcover_path)]:
        # determine target pixel size by determining length of degree
        task = task_graph.add_task(
            func=pygeoprocessing.warp_raster,
            args=(
                base_path, (utm_pixel_size, -utm_pixel_size), target_path,
                'nearest'),
            kwargs={'target_sr_wkt': epsg_srs.ExportToWkt()},
            target_path_list=[target_path],
            dependent_task_list=[align_resize_task],
            task_name='warp_raster_%s' % raster_id)
        path_task_id_map[raster_id] = (target_path, task)

    # fill and route dem
    filled_watershed_dem_path = '%s_filled.tif' % ws_prefix
    flow_dir_watershed_dem_path = '%s_flow_dir.tif' % ws_prefix
    fill_pits_task = task_graph.add_task(
        func=pygeoprocessing.routing.fill_pits,
        args=(
            (path_task_id_map['dem'][0], 1), filled_watershed_dem_path,
            flow_dir_watershed_dem_path),
        kwargs={'temp_dir_path': ws_working_dir},
        target_path_list=[
            filled_watershed_dem_path, flow_dir_watershed_dem_path],
        dependent_task_list=[path_task_id_map['dem'][1]],
        task_name='fill_pits_task_%s' % ws_prefix)

    # flow accum dem
    flow_accum_watershed_dem_path = os.path.join(
        ws_working_dir, '%s_flow_accum.tif' % ws_prefix)
    flow_accmulation_task = task_graph.add_task(
        func=pygeoprocessing.routing.flow_accmulation,
        args=(
            (flow_dir_watershed_dem_path, 1), flow_accum_watershed_dem_path),
        kwargs={'temp_dir_path': ws_working_dir},
        target_path_list=[flow_accum_watershed_dem_path],
        dependent_task_list=[fill_pits_task],
        task_name='flow_accmulation_%s' % ws_prefix)

    # reclassify eff_n
    eff_n_lucode_map = dict(
        zip(biophysical_table['ID'], biophysical_table['eff_n']))
    eff_n_raster_path = os.path.join(
        ws_working_dir, '%s_eff_n.tif' % ws_prefix)
    reclassify_eff_n_task = task_graph.add_task(
        func=pygeoprocessing.reclassify_raster,
        args=(
            (path_task_id_map['landcover'][0], 1), eff_n_lucode_map,
            eff_n_raster_path, gdal.GDT_Float32, NODATA),
        target_path_list=[eff_n_raster_path],
        dependent_task_list=[path_task_id_map['landcover'][1]],
        task_name='reclasify_eff_n_%s' % ws_prefix)

    # reclassify load_n
    load_n_raster_path = os.path.join(
        ws_working_dir, '%s_load_n.tif' % ws_prefix)
    load_n_lucode_map = dict(
        zip(biophysical_table['ID'], biophysical_table['load_n']))
    reclassify_load_n_task = task_graph.add_task(
        func=pygeoprocessing.reclassify_raster,
        args=(
            (path_task_id_map['landcover'][0], 1), load_n_lucode_map,
            load_n_raster_path, gdal.GDT_Float32, NODATA),
        target_path_list=[load_n_raster_path],
        dependent_task_list=[path_task_id_map['landcover'][1]],
        task_name='reclasify_load_n_%s' % ws_prefix)

    # calculate modified load (load * precip)
    modified_load_raster_path = os.path.join(
        ws_working_dir, '%s_modified_load.tif' % ws_prefix)
    modified_load_task = task_graph.add_task(
        func=pygeoprocessing.raster_calculator,
        args=(
            [(load_n_raster_path, 1), (path_task_id_map['precip'][0], 1)],
            mult_arrays, modified_load_raster_path, gdal.GDT_Float32, NODATA),
        target_path_list=[modified_load_raster_path],
        dependent_task_list=[
            reclassify_load_n_task, path_task_id_map['precip'][1]],
        task_name='modified_load_%s' % ws_prefix)

    # calculate slope
    target_slope_path = os.path.join(
        ws_working_dir, '%s_slope.tif' % ws_prefix)
    calculate_slope_task = task_graph.add_task(
        func=pygeoprocessing.routing.calculate_slope,
        args=((path_task_id_map['dem'][0], 1), target_slope_path),
        target_path_list=[target_slope_path],
        dependent_task_list=[path_task_id_map['dem'][1]],
        task_name='calculate_slope_%s' % ws_prefix)

    # calculate D_up
    slope_accum_watershed_dem_path = os.path.join(
        ws_working_dir, '%s_s_accum.tif' % ws_prefix)
    slope_accmulation_task = task_graph.add_task(
        func=pygeoprocessing.routing.flow_accmulation,
        args=(
            (flow_dir_watershed_dem_path, 1), slope_accum_watershed_dem_path),
        kwargs={
            'temp_dir_path': ws_working_dir,
            'weight_raster_band_path': (target_slope_path, 1)},
        target_path_list=[slope_accum_watershed_dem_path],
        dependent_task_list=[fill_pits_task, calculate_slope_task],
        task_name='slope_accmulation_%s' % ws_prefix)

    # calculate flow path length down to stream

    # calculate D_dn

    # calculate NDR specific values

    # TODO: make results part of precip name?

    task_graph.close()
    task_graph.join()


def merge_watershed_dems(
        watershed_path, watershed_id, dem_rtree_path, dem_path_index_map_path,
        target_dem_path):
    """Find DEMs that overlap the given watershed polyon by id.

    Parameters:
        watershed_path (string):
        watershed_id (int): feature number to index and overlap.
        dem_rtree_path (string): path to a pickled rtree that maps
            bounding boxes to dem ids.
        dem_path_index_map_path (string): path to a pickled map to maps DEM
            ids from the rtree to filepaths.
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

    LOGGER.debug(watershed_bb)
    dem_rtree = rtree.index.Index(dem_rtree_path)
    LOGGER.debug(dem_rtree.bounds)

    with open(dem_path_index_map_path, 'rb') as dill_file:
        dem_path_index_map = dill.load(dill_file)

    overlapping_dem_list = list(dem_rtree.intersection(watershed_bb))

    if overlapping_dem_list:
        overlapping_dem_path_list = [
            dem_path_index_map[i] for i in overlapping_dem_list]
        LOGGER.debug("%s %s", watershed_id, overlapping_dem_path_list)
        workspace_dir = os.path.dirname(target_dem_path)
        if not os.path.exists(workspace_dir):
            os.makedirs(workspace_dir)
        pygeoprocessing.merge_rasters(
            overlapping_dem_path_list, target_dem_path,
            bounding_box=watershed_bb)
    else:
        LOGGER.debug(
            "no overlapping dems found for %s wsid %d", watershed_path,
            watershed_id)


def build_dem_rtree(dem_path_list, dem_path_index_map_path, dem_rtree_path):
    """Build RTree indexed by FID for points in `wwwiii_vector_path`."""
    LOGGER.info('building rTree %s', dem_rtree_path+'.dat')
    if os.path.exists(dem_rtree_path+'.dat'):
        LOGGER.warn('%s exists so skipping creation.', dem_rtree_path)
        return
    dem_rtree = rtree.index.Index(dem_rtree_path)
    dem_path_index_map = {}
    for dem_id, dem_path in enumerate(dem_path_list):
        raster_info = pygeoprocessing.get_raster_info(dem_path)
        dem_path_index_map[dem_id] = dem_path
        dem_rtree.insert(dem_id, raster_info['bounding_box'])
    dem_rtree.close()
    with open(dem_path_index_map_path, 'wb') as f:
        dill.dump(dem_path_index_map, f)

if __name__ == '__main__':
    main()
