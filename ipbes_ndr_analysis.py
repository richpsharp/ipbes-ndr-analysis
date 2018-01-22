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

# Compiles cython on runtime.
import pyximport
pyximport.install()
import ipbes_ndr_analysis_cython

N_CPUS = 4
DRY_RUN = False
NODATA = -1
IC_NODATA = -9999
USE_AG_LOAD_ID = 999
FLOW_THRESHOLD = 33
RET_LEN = 150.0
K_VAL = 1.0
SCENARIO_LIST = ['cur', 'ssp1', 'ssp3', 'ssp5']
logging.basicConfig(
    format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
    level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')

hdlr = logging.FileHandler('log.txt')
formatter = logging.Formatter(
    '%(asctime)s %(name)-10s %(levelname)-8s %(message)s')
hdlr.setFormatter(formatter)
LOGGER = logging.getLogger('ipbes_ndr_analysis')

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

RESULTS_DIR = os.path.join(
    BASE_DROPBOX_DIR, 'rps_bck_shared_stuff', 'ipbes stuff',
    'ipbes_ndr_results')

WATERSHED_PATH_LIST = glob.glob(
    os.path.join(
        BASE_DROPBOX_DIR, 'ipbes-data',
        'watersheds_globe_HydroSHEDS_15arcseconds', '*.shp'))

TARGET_WORKSPACE = 'ndr_workspace'
TASKGRAPH_DIR = os.path.join(TARGET_WORKSPACE, 'taskgraph_cache')

RTREE_PATH = 'dem_rtree'

DONE_REPORTING = False


def db_to_shapefile(database_path, sleep_time, db_to_shapefile_lock):
    """Converts db to shapefile every `sleep_time` seconds."""
    while True:
        with db_to_shapefile_lock:
            LOGGER.info("reporting results")
            try:
                target_shapefile_path = os.path.join(
                    RESULTS_DIR, 'jan_2018_results.shp')
                if os.path.exists(target_shapefile_path):
                    os.remove(target_shapefile_path)
                try:
                    os.makedirs(RESULTS_DIR)
                except OSError:
                    pass

                # create a shapefile with these fields:
                # ws_id (text)
                # cur_n_exp
                # ssp1_n_exp
                # ssp3_n_exp
                # ssp5_n_exp

                target_sr = osr.SpatialReference()
                target_sr.ImportFromEPSG(4326)

                driver = gdal.GetDriverByName('ESRI Shapefile')
                result_vector = driver.Create(
                    target_shapefile_path, 0, 0, 0, gdal.GDT_Unknown)
                result_layer = result_vector.CreateLayer(
                    os.path.splitext(os.path.basename(target_shapefile_path))[0],
                    target_sr, ogr.wkbPolygon)
                ws_field = ogr.FieldDefn("ws_id", ogr.OFTString)
                ws_field.SetWidth(24)
                result_layer.CreateField(ws_field)

                for scenario in SCENARIO_LIST:
                    scenario_field = ogr.FieldDefn(
                        '%s_n_exp' % scenario, ogr.OFTReal)
                    scenario_field.SetWidth(32)
                    scenario_field.SetPrecision(11)
                    result_layer.CreateField(scenario_field)

                conn = sqlite3.connect(database_path)
                if conn is not None:
                    cursor = conn.cursor()
                    try:
                        cursor.execute(
                            """SELECT ws_prefix_key, geometry_wkt
                                FROM nutrient_export GROUP BY ws_prefix_key;""")
                    except sqlite3.OperationalError:
                        LOGGER.exception('SQL Error in `db_to_shapefile')
                        continue

                    result = cursor.fetchall()
                    for ws_id, ws_geom in result:
                        feature = ogr.Feature(result_layer.GetLayerDefn())
                        feature.SetField('ws_id', ws_id)
                        for scenario in SCENARIO_LIST:
                            cursor.execute(
                                """SELECT total_export FROM nutrient_export
                                WHERE (ws_prefix_key = ? and scenario_key = ?)""",
                                (ws_id, scenario))
                            result = cursor.fetchone()
                            if result is not None:
                                feature.SetField('%s_n_exp' % scenario, result[0])
                        feature.SetGeometry(ogr.CreateGeometryFromWkt(ws_geom))
                        result_layer.CreateFeature(feature)

                result_vector.FlushCache()
                result_layer = None
                result_vector = None

                for old_timestamp_file_path in glob.glob(
                        os.path.join(RESULTS_DIR, '*.txt')):
                    os.remove(old_timestamp_file_path)
                timestring = datetime.datetime.now().strftime("%Y-%m-%d %H_%M_%S")
                timestamp_path = os.path.join(RESULTS_DIR, 'last_update_%s.txt' % (
                    timestring))
                for file_path in glob.glob('results_v_china.*'):
                    dst_file = os.path.join(RESULTS_DIR, file_path)
                    if os.path.exists(dst_file):
                        os.remove(dst_file)
                    shutil.copy2(file_path, dst_file)
                with open(timestamp_path, 'w') as timestamp_file:
                    timestamp_file.write(
                        "Hi, I'm an automatically generated file.\n"
                        "I last updated NDR results on %s.\n" % timestring +
                        "There will be an 'all done.txt' file here when everything "
                        "is done.\n")
            except Exception:
                LOGGER.exception(
                    "There was an exception during results reporting.")
            finally:
                try:
                    result_vector.FlushCache()
                    result_layer = None
                    result_vector = None
                except:
                    pass
                time.sleep(sleep_time)
                if DONE_REPORTING:
                    break


def length_of_degree(lat):
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
        m1 + m2 * math.cos(2 * lat) + m3 * math.cos(4 * lat) +
        m4 * math.cos(6 * lat))
    longlen = abs(
        p1 * math.cos(lat) + p2 * math.cos(3 * lat) + p3 * math.cos(5 * lat))
    return max(latlen, longlen)


class ClampOp(taskgraph.EncapsulatedTaskOp):
    """Clamp non-nodata values to be >= threshold_val."""
    def __init__(self, raster_path_band, threshold_val, target_path):
        super(ClampOp, self).__init__()
        self.raster_path_band = raster_path_band
        self.threshold_val = threshold_val
        self.target_path = target_path

    def __call__(self):
        nodata = pygeoprocessing.get_raster_info(
            self.raster_path_band[0])['nodata'][self.raster_path_band[1]-1]

        def clamp_op(array):
            """Clamp non-nodata in array to >= threshold_val."""
            result = numpy.empty_like(array)
            result[:] = array
            threshold_mask = (array != nodata) & (array <= self.threshold_val)
            result[threshold_mask] = self.threshold_val
            return result

        pygeoprocessing.raster_calculator(
            [self.raster_path_band], clamp_op, self.target_path,
            gdal.GDT_Float32, nodata)


def calculate_ag_load(
        load_n_raster_path, ag_load_raster_path, target_ag_load_path):
    """Add the agricultural load onto the base load.

    Parameters:
        load_n_raster_path (string): path to a base load raster with
            `USE_AG_LOAD_ID` where the pixel should be replaced with the
            managed ag load.
        ag_load_raster_path (string): path to a raster that indicates
            what the ag load is at `USE_AG_LOAD_ID` pixels
        target_ag_load_path (string): generated raster that has the base
            values from `load_n_raster_path` but with the USE_AG_LOAD_IDs
            replaced by `ag_load_raster_path`.

    Returns:
        None.
    """
    def ag_load_op(base_load_n_array, ag_load_array):
        """raster calculator replace USE_AG_LOAD_ID with ag loads."""
        result = numpy.copy(base_load_n_array)
        ag_mask = result == USE_AG_LOAD_ID
        result[ag_mask] = ag_load_array[ag_mask]
        return result

    nodata = pygeoprocessing.get_raster_info(load_n_raster_path)['nodata'][0]

    pygeoprocessing.raster_calculator(
        [(load_n_raster_path, 1), (ag_load_raster_path, 1)],
        ag_load_op, target_ag_load_path,
        gdal.GDT_Float32, nodata)


def result_in_database(database_path, ws_prefix):
    """True if ws_prefix in database."""
    conn = sqlite3.connect(database_path)
    if conn is not None:
        try:
            cursor = conn.cursor()
            for scenario in SCENARIO_LIST:
                cursor.execute(
                    """SELECT total_export FROM nutrient_export
                    WHERE (ws_prefix_key = ? and scenario_key = ?)""", (
                        ws_prefix, scenario))
                result = cursor.fetchone()
                if result is None:
                    return False
            return True
        except sqlite3.OperationalError:
            LOGGER.exception("operational error on %s"% ws_prefix)
            return False
    return False


def aggregate_to_database(
        n_export_raster_path, global_watershed_path, global_watershed_id,
        local_watershed_path, ws_prefix, scenario_key, aggregate_field_name,
        database_lock, target_database_path):
    """Aggregate nutrient load and save to database.

    This function creates a new database if it does not exist and aggregates
        values and inserts a new row if the row does not already exist.

    Parameters:
        n_export_raster_path (string): path to nutrient export raster in
            units of kg/Ha.
        global_watershed_path (string): path to global shapefile where
            `feature_id` is an element. This lets us save the geometry to
            the database.
        global_watershed_id (int): feature id in global watershed that's
            represented in the local watershed.
        local_watershed_path (string): path to a locally projected shapefile
            with the local watershed. This is used for aggregation.
        ws_prefix (string): used to insert in databases
        scenario_key (string): used to insert in database
        aggregate_field_name (string): unique key in local watershed feature
        target_database_path (string): path to SQLite Database. If not exists
            create a table called 'nutrient_export' with fields:
                * ws_prefix_key (string identifying unique global watershed)
                * scenario_key (string identifying which scenario)
                * total_export (float indicating total export from watershed)
                * geometry (geometry saved as Wkt)

    Returns:
        None.
    """
    with database_lock:
        LOGGER.info(
            '********* aggregating results for %s %s', ws_prefix,
            scenario_key)
        # create a database connection
        conn = sqlite3.connect(target_database_path)
        if conn is not None:
            cursor = conn.cursor()
            cursor.execute(
                """SELECT ws_prefix_key, scenario_key FROM nutrient_export
                WHERE (ws_prefix_key = ? and scenario_key = ?)""", (
                    ws_prefix, scenario_key))
            db_result = cursor.fetchone()
            if db_result is not None:
                # already in table, skipping
                LOGGER.debug(
                    "already in table %s %s %s" % (
                        ws_prefix, scenario_key, db_result))
                LOGGER.debug(
                    "result of in table %s" % result_in_database(
                        target_database_path, ws_prefix))
                LOGGER.debug(target_database_path)
                return

            result = pygeoprocessing.zonal_statistics(
                (n_export_raster_path, 1), local_watershed_path,
                aggregate_field_name, polygons_might_overlap=False)

            pixel_area_ha = pygeoprocessing.get_raster_info(
                n_export_raster_path)['mean_pixel_size']**2 * 0.0001
            total_export = result.itervalues().next()['sum'] * pixel_area_ha
            # ran into case where total export was none maybe because of bad data?
            if total_export is None or math.isnan(total_export):
                total_export = -1

            global_watershed_vector = gdal.OpenEx(
                global_watershed_path, gdal.OF_VECTOR)
            global_watershed_layer = global_watershed_vector.GetLayer()
            global_watershed_feature = global_watershed_layer.GetFeature(
                global_watershed_id)
            global_watershed_geom = global_watershed_feature.GetGeometryRef()
            geometry_wkt = global_watershed_geom.ExportToWkt()
            global_watershed_geom = None
            global_watershed_feature = None
            global_watershed_layer = None
            global_watershed_vector = None
            try:
                cursor.execute(
                    """INSERT INTO nutrient_export VALUES (?, ?, ?, ?)""",
                    (ws_prefix, scenario_key, total_export, geometry_wkt))
            except:
                LOGGER.exception('"%s"', total_export)
            conn.commit()
            conn.close()
        else:
            raise IOError(
                "Error! cannot create the database connection.")


def calculate_ndr(downstream_ret_eff_path, ic_path, k_val, target_ndr_path):
    """Calculate NDR raster.

    Parameters:
        downstream_ret_eff_path (string): path to downstream retention
            raster.
        ic_path (string): path to IC raster
        k_val (float): value of k in Eq. 4.
        target_ndr_path (string): path to NDR raster calculated by this func.

    Returns:
        None.
    """
    # calculate ic_0
    ic_raster = gdal.OpenEx(ic_path, gdal.OF_RASTER)
    ic_min, ic_max, _, _ = ic_raster.GetRasterBand(1).GetStatistics(0, 1)
    ic_0 = (ic_max + ic_min) / 2.0

    def ndr_op(downstream_ret_eff_array, ic_array):
        """Calculate NDR from Eq. (4)."""
        result = numpy.empty_like(downstream_ret_eff_array)
        result[:] = NODATA
        valid_mask = (
            downstream_ret_eff_array != NODATA) & (ic_array != IC_NODATA)
        result[valid_mask] = (1 - downstream_ret_eff_array[valid_mask]) / (
            1 + numpy.exp((ic_array[valid_mask] - ic_0) / k_val))
        return result

    pygeoprocessing.raster_calculator(
        [(downstream_ret_eff_path, 1), (ic_path, 1)], ndr_op, target_ndr_path,
        gdal.GDT_Float32, NODATA)


def calc_ic(d_up_array, d_dn_array):
    """Calculate log_10(d_up/d_dn) unless nodata or 0."""
    result = numpy.empty_like(d_up_array)
    result[:] = NODATA
    zero_mask = d_dn_array == 0
    valid_mask = (
        (d_up_array != NODATA) & (d_dn_array != NODATA) & (~zero_mask))
    result[valid_mask] = numpy.log10(
        d_up_array[valid_mask] / d_dn_array[valid_mask])
    result[zero_mask] = 0.0
    return result


def mult_arrays(
        target_raster_path, gdal_type, target_nodata, raster_path_list):
    """Multiply arrays and be careful of nodata values."""
    nodata_array = numpy.array([
        pygeoprocessing.get_raster_info(path)['nodata'][0]
        for path in raster_path_list])

    def _mult_arrays(*array_list):
        """Multiply arrays in array list but block out stacks with NODATA."""
        stack = numpy.stack(array_list)
        valid_mask = (numpy.bitwise_and.reduce(
            [~numpy.isclose(nodata, array)
             for nodata, array in zip(nodata_array, stack)], axis=0))
        n_valid = numpy.count_nonzero(valid_mask)
        broadcast_valid_mask = numpy.broadcast_to(valid_mask, stack.shape)
        valid_stack = stack[broadcast_valid_mask].reshape(
            len(array_list), n_valid)
        result = numpy.empty(array_list[0].shape, dtype=numpy.float64)
        result[:] = NODATA
        result[valid_mask] = numpy.prod(valid_stack, axis=0)
        return result

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in raster_path_list], _mult_arrays,
        target_raster_path, gdal_type, target_nodata)


def div_arrays(num_array, denom_array):
    """Calculate num / denom except when denom = 0 or nodata."""
    result = numpy.empty_like(num_array)
    result[:] = NODATA
    valid_mask = (
        (num_array != NODATA) & (denom_array != NODATA) & (denom_array != 0))
    result[valid_mask] = num_array[valid_mask] / denom_array[valid_mask]
    return result


class MultByScalar(taskgraph.EncapsulatedTaskOp):
    """Multiply raster by a scalar, ignore nodata."""
    def __init__(self, raster_path_band, scalar, target_nodata, target_path):
        super(MultByScalar, self).__init__()
        self.raster_path_band = raster_path_band
        self.scalar = scalar
        self.target_nodata = target_nodata
        self.target_path = target_path

    def __call__(self):
        nodata = pygeoprocessing.get_raster_info(
            self.raster_path_band[0])['nodata'][self.raster_path_band[1]-1]

        def mult_by_scalar(array):
            """Multiply non-nodta values by self.scalar"""
            result = numpy.empty_like(array)
            result[:] = self.target_nodata
            valid_mask = array != nodata
            result[valid_mask] = array[valid_mask] * self.scalar
            return result

        pygeoprocessing.raster_calculator(
            [self.raster_path_band], mult_by_scalar, self.target_path,
            gdal.GDT_Float32, self.target_nodata)

class DUpOp(taskgraph.EncapsulatedTaskOp):
    """Calculate D_up from Equation 7 of NDR user's guide.

    Given a flow accumulation raster, slope accumulation raster, and pixel
    size, we can calculate avg(S)*sqrt(A) for each pixel
        avg(S) = slope_accum / flow_accum
        A = flow_accum * sqrt(flow_accum * pixel_area**2)
    """
    def __init__(
            self, pixel_area, slope_accum_raster_path,
            flow_accum_raster_path, target_d_up_raster_path):
        """Parameters:
            pixel_area (float): area of input raster pixel in m^2.
            slope_accum_raster_path (string): path to slope accumulation
                raster.
            flow_accum_raster_path (string): path to flow accumulation raster.
            target_d_up_raster_path (string): path to target d_up raster path
                created by a call to __call__.
        """
        super(DUpOp, self).__init__()
        self.pixel_area = pixel_area
        self.slope_accum_raster_path = slope_accum_raster_path
        self.flow_accum_raster_path = flow_accum_raster_path
        self.target_d_up_raster_path = target_d_up_raster_path

    def __call__(self):
        flow_accum_nodata = pygeoprocessing.get_raster_info(
            self.flow_accum_raster_path)['nodata'][0]

        def d_up_op(slope_accum_array, flow_accmulation_array):
            """Mult average upslope by sqrt of upslope area."""
            result = numpy.empty_like(slope_accum_array)
            result[:] = NODATA
            valid_mask = flow_accmulation_array != flow_accum_nodata
            result[valid_mask] = (
                slope_accum_array[valid_mask] /
                flow_accmulation_array[valid_mask]) * numpy.sqrt(
                    flow_accmulation_array[valid_mask] * self.pixel_area)
            return result

        pygeoprocessing.raster_calculator(
            [(self.slope_accum_raster_path, 1),
             (self.flow_accum_raster_path, 1)], d_up_op,
            self.target_d_up_raster_path, gdal.GDT_Float32, NODATA)


def main():
    """Entry point."""
    logging.basicConfig(
        format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
        level=logging.WARN, datefmt='%m/%d/%Y %H:%M:%S ')
    task_priority = 0  # we'll use this to schedule deep rather than wide
    try:
        os.makedirs(TARGET_WORKSPACE)
    except OSError:
        pass

    multiprocessing_manager = multiprocessing.Manager()
    database_lock = multiprocessing_manager.Lock()

    # load biophysical table first, it's used a lot below
    represenative_ndr_biophysical_table_path = os.path.join(
        BASE_DROPBOX_DIR, 'ipbes-data', 'NDR_representative_table.csv')
    biophysical_table = pandas.read_csv(
        represenative_ndr_biophysical_table_path)
    # clean up biophysical table
    biophysical_table = biophysical_table.fillna(0)
    biophysical_table.ix[biophysical_table['load_n'] == 'use raster', 'load_n'] = (
        USE_AG_LOAD_ID)
    biophysical_table['load_n'] = biophysical_table['load_n'].apply(
        pandas.to_numeric)
    task_graph = taskgraph.TaskGraph(
        TASKGRAPH_DIR, N_CPUS, 5.0, DRY_RUN)

    database_path = os.path.join(TARGET_WORKSPACE, 'ipbes_ndr_results.db')

    sql_create_projects_table = (
        """ CREATE TABLE IF NOT EXISTS nutrient_export (
            ws_prefix_key TEXT NOT NULL,
            scenario_key TEXT NOT NULL,
            total_export REAL NOT NULL,
            geometry_wkt TEXT NOT NULL,
            PRIMARY KEY (ws_prefix_key, scenario_key)
        ); """)

    # create a database connection
    conn = sqlite3.connect(database_path)
    cursor = conn.cursor()
    cursor.execute(sql_create_projects_table)
    cursor.execute("""CREATE INDEX IF NOT EXISTS idx_nutrient_export
        ON nutrient_export (ws_prefix_key, scenario_key);""")

    db_to_shapefile_lock = threading.Lock()
    db_to_shapefile_thread = threading.Thread(
        target=db_to_shapefile,
        args=(database_path, 60.0, db_to_shapefile_lock))
    db_to_shapefile_thread.daemon = True
    #db_to_shapefile_thread.start()

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
        task_name='build_dem_rtree',
        priority=task_priority)
    task_priority -= 1

    global_watershed_path_list = glob.glob(
        os.path.join(
            BASE_DROPBOX_DIR, 'ipbes-data',
            'watersheds_globe_HydroSHEDS_15arcseconds', '*.shp'))
    watershed_priority_queue = []
    for global_watershed_path in global_watershed_path_list:
        watershed_vector = gdal.OpenEx(global_watershed_path, gdal.OF_VECTOR)
        watershed_layer = watershed_vector.GetLayer()
        for watershed_id in xrange(watershed_layer.GetFeatureCount()):

            watershed_feature = watershed_layer.GetFeature(watershed_id)
            feature_geom = watershed_feature.GetGeometryRef()
            watershed_area = feature_geom.GetArea()
            if watershed_area < 0.03:
                #  0.04 square degrees is a healthy underapproximation of
                # 100 sq km which is about the minimum watershed size we'd
                # want.
                continue
            feature_centroid = feature_geom.Centroid()
            feature_geom = None
            watershed_basename = os.path.splitext(
                os.path.basename(global_watershed_path))[0]
            ws_prefix = 'ws_%s_%d' % (watershed_basename, watershed_id)
            if result_in_database(database_path, ws_prefix):
                LOGGER.info("%s already reported, skipping", ws_prefix)
                continue
            heapq.heappush(watershed_priority_queue, (
                -watershed_area, feature_centroid, global_watershed_path,
                watershed_id, ws_prefix))
            feature_geom = None
            watershed_feature = None
        watershed_layer = None
    watershed_vector = None

    while watershed_priority_queue:
        # iterate from smallest to largest watershed
        (watershed_area, feature_centroid,
         global_watershed_path, watershed_id, ws_prefix) = heapq.heappop(
            watershed_priority_queue)
        LOGGER.debug('****** %f %d %s', watershed_area, watershed_id, global_watershed_path)

        # make a few subdirectories so we don't explode on directories
        last_digits = '%.4d' % watershed_id
        ws_working_dir = os.path.join(
            TARGET_WORKSPACE, last_digits[-1], last_digits[-2],
            last_digits[-3], last_digits[-4],
            "%s_working_dir" % ws_prefix)
        watershed_dem_path = os.path.join(
            ws_working_dir, 'ws_%s_dem.tif' % ws_prefix)
        merge_watershed_dems_task = task_graph.add_task(
            func=merge_watershed_dems,
            args=(
                global_watershed_path, watershed_id, dem_rtree_path,
                dem_path_index_map_path, watershed_dem_path),
            target_path_list=[watershed_dem_path],
            dependent_task_list=[build_dem_rtree_task],
            task_name='merge_watershed_dems_%s' % ws_prefix,
            priority=task_priority)
        task_priority -= 1

        utm_code = (math.floor((feature_centroid.GetX() + 180)/6) % 60) + 1
        lat_code = 6 if feature_centroid.GetY() > 0 else 7
        epsg_code = int('32%d%02d' % (lat_code, utm_code))
        epsg_srs = osr.SpatialReference()
        epsg_srs.ImportFromEPSG(epsg_code)
        utm_pixel_size = 500 #abs(dem_pixel_size[0]) * length_of_degree(feature_centroid.GetY())

        local_watershed_path = os.path.join(ws_working_dir, '%s.shp' % ws_prefix)
        if os.path.exists(local_watershed_path):
            os.remove(local_watershed_path)
        reproject_watershed_task = task_graph.add_task(
            func=reproject_vector_feature,
            args=(
                global_watershed_path, epsg_srs.ExportToWkt(), watershed_id,
                local_watershed_path),
            target_path_list=[local_watershed_path],
            task_name='project_watershed_%s' % ws_prefix,
            priority=task_priority)
        task_priority -= 1

        landcover_raster_path = os.path.join(
            BASE_DROPBOX_DIR, 'ipbes-data',
            'GLOBIO4_landuse_10sec_tifs_20171207_Idiv', 'Current2015',
            'Globio4_landuse_10sec_2015_cropint.tif')

        precip_raster_path_cur = os.path.join(
            BASE_DROPBOX_DIR, 'ipbes-data',
            'precip_globe_WorldClim_30arcseconds.tif')

        precip_raster_path_ssp1 = os.path.join(
            BASE_DROPBOX_DIR, 'ipbes-data',
            'ssp1_3_5_annual_precip_scenarios_for_ndr', 'he26pr50.tif')

        precip_raster_path_ssp3 = os.path.join(
            BASE_DROPBOX_DIR, 'ipbes-data',
            'ssp1_3_5_annual_precip_scenarios_for_ndr', 'he60pr50.tif')

        precip_raster_path_ssp5 = os.path.join(
            BASE_DROPBOX_DIR, 'ipbes-data',
            'ssp1_3_5_annual_precip_scenarios_for_ndr', 'he85pr50.tif')

        ag_load_raster_path_cur = os.path.join(
            BASE_DROPBOX_DIR, 'ipbes-data',
            'ag_load_scenarios', 'ssp1_2015_ag_load.tif')
        ag_load_raster_path_ssp1 = os.path.join(
            BASE_DROPBOX_DIR, 'ipbes-data',
            'ag_load_scenarios', 'ssp1_2050_ag_load.tif')
        ag_load_raster_path_ssp3 = os.path.join(
            BASE_DROPBOX_DIR, 'ipbes-data',
            'ag_load_scenarios', 'ssp3_2050_ag_load.tif')
        ag_load_raster_path_ssp5 = os.path.join(
            BASE_DROPBOX_DIR, 'ipbes-data',
            'ag_load_scenarios', 'ssp5_2050_ag_load.tif')

        base_raster_path_list = [
            watershed_dem_path,
            landcover_raster_path,
            precip_raster_path_cur,
            precip_raster_path_ssp1,
            precip_raster_path_ssp3,
            precip_raster_path_ssp5,
            ag_load_raster_path_cur,
            ag_load_raster_path_ssp1,
            ag_load_raster_path_ssp3,
            ag_load_raster_path_ssp5]

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
                ['nearest', 'mode', 'nearest', 'nearest', 'nearest', 'nearest',
                 'nearest', 'nearest', 'nearest', 'nearest'],
                dem_pixel_size, 'intersection'),
            target_path_list=aligned_path_list,
            dependent_task_list=[merge_watershed_dems_task],
            task_name='align_resize_task_%s' % ws_prefix,
            priority=task_priority)
        task_priority -= 1

        utm_dem_path = os.path.join(
            ws_working_dir, '%s_%s_dem.tif' % (ws_prefix, epsg_code))
        utm_landcover_path = os.path.join(
            ws_working_dir, '%s_%s_landcover.tif' % (ws_prefix, epsg_code))
        utm_precip_scenario_path = os.path.join(
            ws_working_dir, '%s_%s_%s.tif' % (ws_prefix, epsg_code, '%s_precip'))
        utm_ag_load_scenario_path = os.path.join(
            ws_working_dir, '%s_%s_%s.tif' % (ws_prefix, epsg_code, '%s_ag_load'))

        path_task_id_map = {}
        for raster_id, base_path, target_path in [
                ('dem', aligned_path_list[0], utm_dem_path),
                ('landcover', aligned_path_list[1], utm_landcover_path),
                ('precip_cur', aligned_path_list[2],
                 utm_precip_scenario_path % 'cur'),
                ('precip_ssp1', aligned_path_list[3],
                 utm_precip_scenario_path % 'ssp1'),
                ('precip_ssp3', aligned_path_list[4],
                 utm_precip_scenario_path % 'ssp3'),
                ('precip_ssp5', aligned_path_list[5],
                 utm_precip_scenario_path % 'ssp5'),
                ('ag_load_cur', aligned_path_list[6],
                 utm_ag_load_scenario_path % 'cur'),
                ('ag_load_ssp1', aligned_path_list[7],
                 utm_ag_load_scenario_path % 'ssp1'),
                ('ag_load_ssp3', aligned_path_list[8],
                 utm_ag_load_scenario_path % 'ssp3'),
                ('ag_load_ssp5', aligned_path_list[9],
                 utm_ag_load_scenario_path % 'ssp5')]:
            # determine target pixel size by determining length of degree
            task = task_graph.add_task(
                func=pygeoprocessing.warp_raster,
                args=(
                    base_path, (utm_pixel_size, -utm_pixel_size), target_path,
                    'nearest'),
                kwargs={'target_sr_wkt': epsg_srs.ExportToWkt()},
                target_path_list=[target_path],
                dependent_task_list=[align_resize_task],
                task_name='warp_raster_%s' % raster_id,
                priority=task_priority)
            task_priority -= 1
            path_task_id_map[raster_id] = (target_path, task)

        # fill and route dem
        filled_watershed_dem_path = os.path.join(
            ws_working_dir, '%s_filled.tif' % ws_prefix)
        flow_dir_path = os.path.join(
            ws_working_dir, '%s_flow_dir.tif' % ws_prefix)
        fill_pits_task = task_graph.add_task(
            func=pygeoprocessing.routing.fill_pits,
            args=(
                (path_task_id_map['dem'][0], 1),
                filled_watershed_dem_path, flow_dir_path),
            kwargs={'temp_dir_path': ws_working_dir},
            target_path_list=[
                filled_watershed_dem_path, flow_dir_path],
            dependent_task_list=[path_task_id_map['dem'][1]],
            task_name='fill_pits_task_%s' % ws_prefix,
            priority=task_priority)
        task_priority -= 1

        # flow accum dem
        flow_accum_path = os.path.join(
            ws_working_dir, '%s_flow_accum.tif' % ws_prefix)
        flow_accum_task = task_graph.add_task(
            func=pygeoprocessing.routing.flow_accmulation,
            args=(
                (flow_dir_path, 1), flow_accum_path),
            kwargs={'temp_dir_path': ws_working_dir},
            target_path_list=[flow_accum_path],
            dependent_task_list=[fill_pits_task],
            task_name='flow_accmulation_%s' % ws_prefix,
            priority=task_priority)
        task_priority -= 1

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
            task_name='reclasify_eff_n_%s' % ws_prefix,
            priority=task_priority)
        task_priority -= 1

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
            task_name='reclasify_load_n_%s' % ws_prefix,
            priority=task_priority)
        task_priority -= 1

        # calculate slope
        slope_raster_path = os.path.join(
            ws_working_dir, '%s_slope.tif' % ws_prefix)
        calculate_slope_task = task_graph.add_task(
            func=pygeoprocessing.routing.calculate_slope,
            args=((path_task_id_map['dem'][0], 1), slope_raster_path),
            target_path_list=[slope_raster_path],
            dependent_task_list=[path_task_id_map['dem'][1]],
            task_name='calculate_slope_%s' % ws_prefix,
            priority=task_priority)
        task_priority -= 1

        clamp_slope_raster_path = os.path.join(
            ws_working_dir, '%s_clamp_slope.tif' % ws_prefix)
        clamp_slope_task = task_graph.add_task(
            func=ClampOp(
                (slope_raster_path, 1), 0.005, clamp_slope_raster_path),
            target_path_list=[clamp_slope_raster_path],
            dependent_task_list=[calculate_slope_task],
            task_name='clamp_slope_%s' % ws_prefix,
            priority=task_priority)
        task_priority -= 1

        # calculate D_up
        slope_accum_watershed_dem_path = os.path.join(
            ws_working_dir, '%s_s_accum.tif' % ws_prefix)
        slope_accmulation_task = task_graph.add_task(
            func=pygeoprocessing.routing.flow_accmulation,
            args=(
                (flow_dir_path, 1), slope_accum_watershed_dem_path),
            kwargs={
                'temp_dir_path': ws_working_dir,
                'weight_raster_path_band': (clamp_slope_raster_path, 1)},
            target_path_list=[slope_accum_watershed_dem_path],
            dependent_task_list=[fill_pits_task, clamp_slope_task],
            task_name='slope_accmulation_%s' % ws_prefix,
            priority=task_priority)
        task_priority -= 1

        d_up_raster_path = os.path.join(ws_working_dir, '%s_d_up.tif' % ws_prefix)
        d_up_task = task_graph.add_task(
            func=DUpOp(
                utm_pixel_size**2, slope_accum_watershed_dem_path,
                flow_accum_path, d_up_raster_path),
            target_path_list=[d_up_raster_path],
            dependent_task_list=[slope_accmulation_task, flow_accum_task],
            task_name='d_up_%s' % ws_prefix,
            priority=task_priority)
        task_priority -= 1

        # calculate flow path in pixels length down to stream
        pixel_flow_length_raster_path = os.path.join(
            ws_working_dir, '%s_pixel_flow_length.tif' % ws_prefix)
        downstream_flow_length_task = task_graph.add_task(
            func=pygeoprocessing.routing.downstream_flow_length,
            args=(
                (flow_dir_path, 1),
                (flow_accum_path, 1), FLOW_THRESHOLD,
                pixel_flow_length_raster_path),
            kwargs={'temp_dir_path': ws_working_dir},
            target_path_list=[pixel_flow_length_raster_path],
            dependent_task_list=[fill_pits_task, flow_accum_task],
            task_name='downstream_pixel_flow_length_%s' % ws_prefix,
            priority=task_priority)
        task_priority -= 1

        # calculate real flow_path (flow length * pixel size)
        downstream_flow_distance_path = os.path.join(
            ws_working_dir, '%s_m_flow_length.tif' % ws_prefix)
        downstream_flow_distance_task = task_graph.add_task(
            func=MultByScalar(
                (pixel_flow_length_raster_path, 1), utm_pixel_size, NODATA,
                downstream_flow_distance_path),
            target_path_list=[downstream_flow_distance_path],
            dependent_task_list=[downstream_flow_length_task],
            task_name='downstream_m_flow_dist_%s' % ws_prefix,
            priority=task_priority)
        task_priority -= 1

        # calculate downstream distance / downstream slope
        d_dn_per_pixel_path = os.path.join(
            ws_working_dir, '%s_d_dn_per_pixel.tif' % ws_prefix)
        d_dn_per_pixel_task = task_graph.add_task(
            func=pygeoprocessing.raster_calculator,
            args=(
                [(downstream_flow_distance_path, 1),
                 (clamp_slope_raster_path, 1)],
                div_arrays, d_dn_per_pixel_path, gdal.GDT_Float32, NODATA),
            target_path_list=[d_dn_per_pixel_path],
            dependent_task_list=[
                downstream_flow_distance_task, clamp_slope_task],
            task_name='d_dn_per_pixel_%s' % ws_prefix,
            priority=task_priority)
        task_priority -= 1

        # calculate D_dn: downstream sum of distance / downstream slope
        d_dn_raster_path = os.path.join(
            ws_working_dir, '%s_d_dn.tif' % ws_prefix)
        d_dn_task = task_graph.add_task(
            func=pygeoprocessing.routing.downstream_flow_length,
            args=(
                (flow_dir_path, 1),
                (flow_accum_path, 1), FLOW_THRESHOLD,
                d_dn_raster_path),
            kwargs={
                'temp_dir_path': ws_working_dir,
                'weight_raster_path_band': (d_dn_per_pixel_path, 1)
                },
            target_path_list=[d_dn_raster_path],
            dependent_task_list=[
                fill_pits_task, flow_accum_task, d_dn_per_pixel_task],
            task_name='d_dn_%s' % ws_prefix,
            priority=task_priority)
        task_priority -= 1

        # calculate IC
        ic_path = os.path.join(ws_working_dir, '%s_ic.tif' % ws_prefix)
        ic_task = task_graph.add_task(
            func=pygeoprocessing.raster_calculator,
            args=(
                [(d_up_raster_path, 1), (d_dn_raster_path, 1)],
                calc_ic, ic_path, gdal.GDT_Float32, IC_NODATA),
            target_path_list=[ic_path],
            dependent_task_list=[d_up_task, d_dn_task],
            task_name='ic_%s' % ws_prefix,
            priority=task_priority)
        task_priority -= 1

        # calculate eff_i
        downstream_ret_eff_path = os.path.join(
            ws_working_dir, '%s_downstream_ret_eff.tif' % ws_prefix)
        downstream_ret_eff_task = task_graph.add_task(
            func=ipbes_ndr_analysis_cython.calculate_downstream_ret_eff,
            args=(
                (flow_dir_path, 1), (flow_accum_path, 1), (eff_n_raster_path, 1),
                FLOW_THRESHOLD, RET_LEN, downstream_ret_eff_path),
            kwargs={'temp_dir_path': ws_working_dir},
            target_path_list=[downstream_ret_eff_path],
            dependent_task_list=[
                fill_pits_task, flow_accum_task, reclassify_eff_n_task],
            task_name='downstream_ret_eff_%s' % ws_prefix,
            priority=task_priority)
        task_priority -= 1

        # calculate NDR specific values
        ndr_path = os.path.join(ws_working_dir, '%s_ndr.tif' % ws_prefix)
        ndr_task = task_graph.add_task(
            func=calculate_ndr,
            args=(downstream_ret_eff_path, ic_path, K_VAL, ndr_path),
            target_path_list=[ndr_path],
            dependent_task_list=[downstream_ret_eff_task, ic_task],
            task_name='ndr_task_%s' % ws_prefix,
            priority=task_priority)
        task_priority -= 1

        for scenario_key in SCENARIO_LIST:
            # calculate modified load (load * precip)
            # calculate scenario AG load
            scenario_ag_load_path = os.path.join(
                ws_working_dir, '%s_%s_ag_load.tif' % (ws_prefix, scenario_key))
            scenario_load_task = task_graph.add_task(
                func=calculate_ag_load,
                args=(
                    load_n_raster_path,
                    path_task_id_map['ag_load_%s' % scenario_key][0],
                    scenario_ag_load_path),
                target_path_list=[scenario_ag_load_path],
                dependent_task_list=[
                    reclassify_load_n_task,
                    path_task_id_map['ag_load_%s' % scenario_key][1]],
                task_name='scenario_load_%s_%s' % (ws_prefix, scenario_key),
                priority=task_priority)
            task_priority -= 1

            # calculate modified load (load * precip)
            modified_load_raster_path = os.path.join(
                ws_working_dir, '%s_%s_modified_load.tif' % (
                    ws_prefix, scenario_key))
            modified_load_task = task_graph.add_task(
                func=mult_arrays,
                args=(
                    modified_load_raster_path, gdal.GDT_Float32,
                    NODATA, [
                        scenario_ag_load_path,
                        path_task_id_map['precip_%s' % scenario_key][0]]),
                target_path_list=[modified_load_raster_path],
                dependent_task_list=[
                    scenario_load_task,
                    path_task_id_map['precip_%s' % scenario_key][1]],
                task_name='modified_load_%s' % ws_prefix,
                priority=task_priority)
            task_priority -= 1

            n_export_raster_path = os.path.join(
                ws_working_dir, '%s_%s_n_export.tif' % (
                    ws_prefix, scenario_key))
            n_export_task = task_graph.add_task(
                func=mult_arrays,
                args=(
                    n_export_raster_path, gdal.GDT_Float32, NODATA,
                    [modified_load_raster_path, ndr_path]),
                target_path_list=[n_export_raster_path],
                dependent_task_list=[modified_load_task, ndr_task],
                priority=task_priority,
                task_name='n_export_%s_%s' % (scenario_key, ws_prefix))
            task_priority -= 1

            aggregate_result_task = task_graph.add_task(
                func=aggregate_to_database,
                args=(
                    n_export_raster_path, global_watershed_path, watershed_id,
                    local_watershed_path, ws_prefix, scenario_key, 'BASIN_ID',
                    database_lock, database_path),
                dependent_task_list=[n_export_task, reproject_watershed_task],
                priority=task_priority,
                task_name='aggregate_result_%s_%s' % (scenario_key, ws_prefix))
            task_priority -= 1
    task_graph.close()
    task_graph.join()
    # report one last time
    DONE_REPORTING = True
    db_to_shapefile(database_path, 0.0, db_to_shapefile_lock)
    with open(
            os.path.join(RESULTS_DIR, 'ALL_DONE.txt'), 'w') as timestamp_file:
        timestamp_file.write(
            "Hi, I finished the NDR analysis at %s. Check out "
            "results.shp\n" % datetime.datetime.now().strftime(
                "%Y-%m-%d %H_%M_%S"))


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
    watershed_vector = gdal.OpenEx(watershed_path, gdal.OF_VECTOR)
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
        try:
            os.makedirs(workspace_dir)
        except OSError:
            pass
        pygeoprocessing.merge_rasters(
            overlapping_dem_path_list, target_dem_path,
            expected_nodata=-32768.0, bounding_box=watershed_bb)
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


def reproject_vector_feature(
        base_vector_path, target_wkt, feature_id, target_path):
    """Reproject a single OGR DataSource feature.

    Transforms the features of the base vector to the desired output
    projection in a new ESRI Shapefile.

    Parameters:
        base_vector_path (string): Path to the base shapefile to transform.
        target_wkt (string): the desired output projection in Well Known Text
            (by layer.GetSpatialRef().ExportToWkt())
        feature_id (int): the feature to reproject and copy.
        target_path (string): the filepath to the transformed shapefile

    Returns:
        None
    """
    base_vector = gdal.OpenEx(base_vector_path, gdal.OF_VECTOR)

    # if this file already exists, then remove it
    if os.path.isfile(target_path):
        LOGGER.warn(
            "reproject_vector: %s already exists, removing and overwriting",
            target_path)
        os.remove(target_path)

    try:
        os.makedirs(os.path.dirname(target_path))
    except OSError:
        pass

    target_sr = osr.SpatialReference(target_wkt)

    # create a new shapefile from the orginal_datasource
    target_driver = gdal.GetDriverByName('ESRI Shapefile')
    target_vector = target_driver.Create(
        target_path, 0, 0, 0, gdal.GDT_Unknown)

    layer = base_vector.GetLayer()
    layer_dfn = layer.GetLayerDefn()

    # Create new layer for target_vector using same name and
    # geometry type from base vector but new projection
    target_layer = target_vector.CreateLayer(
        layer_dfn.GetName(), target_sr, layer_dfn.GetGeomType())

    # Get the number of fields in original_layer
    original_field_count = layer_dfn.GetFieldCount()

    # For every field, create a duplicate field in the new layer
    for fld_index in xrange(original_field_count):
        original_field = layer_dfn.GetFieldDefn(fld_index)
        target_field = ogr.FieldDefn(
            original_field.GetName(), original_field.GetType())
        target_layer.CreateField(target_field)

    # Get the SR of the original_layer to use in transforming
    base_sr = layer.GetSpatialRef()

    # Create a coordinate transformation
    coord_trans = osr.CoordinateTransformation(base_sr, target_sr)

    # Copy all of the features in layer to the new shapefile
    error_count = 0
    base_feature = layer.GetFeature(feature_id)
    geom = base_feature.GetGeometryRef()
    if geom is None:
        # we encountered this error occasionally when transforming clipped
        # global polygons.  Not clear what is happening but perhaps a
        # feature was retained that otherwise wouldn't have been included
        # in the clip
        raise ValueError("Geometry is None.")

    # Transform geometry into format desired for the new projection
    error_code = geom.Transform(coord_trans)
    if error_code != 0:  # error
        # this could be caused by an out of range transformation
        # whatever the case, don't put the transformed poly into the
        # output set
        raise ValueError("Unable to reproject geometry.")

    # Copy original_datasource's feature and set as new shapes feature
    target_feature = ogr.Feature(target_layer.GetLayerDefn())
    target_feature.SetGeometry(geom)

    # For all the fields in the feature set the field values from the
    # source field
    for fld_index in xrange(target_feature.GetFieldCount()):
        target_feature.SetField(
            fld_index, base_feature.GetField(fld_index))

    target_layer.CreateFeature(target_feature)
    target_feature = None
    base_feature = None
    if error_count > 0:
        LOGGER.warn(
            '%d features out of %d were unable to be transformed and are'
            ' not in the output vector at %s', error_count,
            layer.GetFeatureCount(), target_path)
    layer = None
    base_vector = None

if __name__ == '__main__':
    main()
