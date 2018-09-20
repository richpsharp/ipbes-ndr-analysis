"""Script to manage NDR runs for IPBES project."""
import zipfile
import sys
import shutil
import datetime
import logging
import os
import glob
import math
import sqlite3

import reproduce.utils
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
TASKGRAPH_REPORTING_FREQUENCY = 5.0
TASKGRAPH_DELAYED_START = False
NODATA = -1
IC_NODATA = -9999
USE_AG_LOAD_ID = 999
FLOW_THRESHOLD = 1000
RET_LEN = 150.0
K_VAL = 1.0

logging.basicConfig(
    level=logging.INFO,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(pathname)s.%(funcName)s:%(lineno)d] %(message)s'),
    stream=sys.stdout)
LOGGER = logging.getLogger(__name__)

BUCKET_DOWNLOAD_DIR = 'bucket_sync'
CHURN_DIR = 'churn'
RTREE_PATH = 'dem_rtree'
WATERSHED_PROCESSING_DIR = 'watershed_processing'

# The following paths will be relative to the workspace directory
LANDUSE_DIR = 'globio_landuse_scenarios'
LANDCOVER_RASTER_PATHS = {
    '1850': f"{LANDUSE_DIR}/Globio4_landuse_10sec_1850.tif",
    '1900': f"{LANDUSE_DIR}/Globio4_landuse_10sec_1900.tif",
    '1910': f"{LANDUSE_DIR}/Globio4_landuse_10sec_1910.tif",
    '1945': f"{LANDUSE_DIR}/Globio4_landuse_10sec_1945.tif",
    '1980': f"{LANDUSE_DIR}/Globio4_landuse_10sec_1980.tif",
    '2015': f"{LANDUSE_DIR}/Globio4_landuse_10sec_2015.tif",
    'ssp1': f"{LANDUSE_DIR}/Globio4_landuse_10sec_2050_cropint_SSP1.tif",
    'ssp3': f"{LANDUSE_DIR}/Globio4_landuse_10sec_2050_cropint_SSP3.tif",
    'ssp5': f"{LANDUSE_DIR}/Globio4_landuse_10sec_2050_cropint_SSP5.tif",
}

PRECIP_DIR = 'precip_scenarios'
PRECIP_RASTER_PATHS = {
    '1900': f'{PRECIP_DIR}/precip_1900.tif',
    '1910': f'{PRECIP_DIR}/precip_1910.tif',
    '1945': f'{PRECIP_DIR}/precip_1945.tif',
    '1980': f'{PRECIP_DIR}/precip_1980.tif',
    '2015': f'{PRECIP_DIR}/precip_2015.tif',
    'ssp1': f'{PRECIP_DIR}/ssp1_2050.tif',
    'ssp3': f'{PRECIP_DIR}/ssp3_2050.tif',
    'ssp5': f'{PRECIP_DIR}/ssp5_2050.tif',
    'ssp1_he26pr50': f'{PRECIP_DIR}/he26pr50.tif',
    'ssp3_he60pr50': f'{PRECIP_DIR}/he60pr50.tif',
    'ssp5_he85pr50': f'{PRECIP_DIR}/he85pr50.tif',
}

AG_LOAD_DIR = 'ag_load_scenarios'
AG_RASTER_PATHS = {
    '1850': f'{AG_LOAD_DIR}/1850_ag_load.tif',
    '1900': f'{AG_LOAD_DIR}/1900_ag_load.tif',
    '1920': f'{AG_LOAD_DIR}/1920_ag_load.tif',
    '1945': f'{AG_LOAD_DIR}/1945_ag_load.tif',
    '1980': f'{AG_LOAD_DIR}/1980_ag_load.tif',
    '2015': f'{AG_LOAD_DIR}/2015_ag_load.tif',
    'ssp1': f'{AG_LOAD_DIR}/ssp1_2050_ag_load.tif',
    'ssp3': f'{AG_LOAD_DIR}/ssp3_2050_ag_load.tif',
    'ssp5': f'{AG_LOAD_DIR}/ssp5_2050_ag_load.tif',
}

def db_to_shapefile(database_path):
    """Converts db to shapefile every `sleep_time` seconds."""
    LOGGER.info("reporting results")
    try:
        target_shapefile_path = os.path.join(
            TASKGRAPH_DIR, '25_jan_2018_results_ha.gpkg')
        if os.path.exists(target_shapefile_path):
            os.remove(target_shapefile_path)
        try:
            os.makedirs(TASKGRAPH_DIR)
        except OSError:
            pass

        # create a shapefile with these fields:
        # ws_id (text)
        # cur_n_exp
        # ssp1_n_exp
        # ssp3_n_exp
        # ssp5_n_exp

        wgs84_sr = osr.SpatialReference()
        wgs84_sr.ImportFromEPSG(4326)

        driver = gdal.GetDriverByName('ESRI Shapefile')
        result_vector = driver.Create(
            target_shapefile_path, 0, 0, 0, gdal.GDT_Unknown)
        result_layer = result_vector.CreateLayer(
            os.path.splitext(os.path.basename(target_shapefile_path))[0],
            wgs84_sr, ogr.wkbPolygon)
        ws_field = ogr.FieldDefn("ws_id", ogr.OFTString)
        ws_field.SetWidth(24)
        result_layer.CreateField(ws_field)

        area_field = ogr.FieldDefn(
            'area_ha', ogr.OFTReal)
        area_field.SetWidth(32)
        area_field.SetPrecision(11)
        result_layer.CreateField(area_field)

        for scenario in SCENARIO_LIST:
            scenario_field = ogr.FieldDefn(
                '%snex' % scenario, ogr.OFTReal)
            scenario_field.SetWidth(32)
            scenario_field.SetPrecision(11)
            result_layer.CreateField(scenario_field)

            scenario_field = ogr.FieldDefn(
                '%snexha' % scenario, ogr.OFTReal)
            scenario_field.SetWidth(32)
            scenario_field.SetPrecision(11)
            result_layer.CreateField(scenario_field)

        conn = sqlite3.connect(database_path)
        if conn is not None:
            cursor = conn.cursor()
            scenaro_cursor = conn.cursor()
            try:
                cursor.execute(
                    """SELECT ws_prefix_key, geometry_wkt
                        FROM nutrient_export GROUP BY ws_prefix_key;""")
            except sqlite3.OperationalError:
                LOGGER.exception('SQL Error in `db_to_shapefile')
            for ws_id, ws_geom in cursor:
                feature = ogr.Feature(result_layer.GetLayerDefn())
                feature.SetField('ws_id', ws_id)

                feature_geom = ogr.CreateGeometryFromWkt(ws_geom)
                feature.SetGeometry(feature_geom)
                feature_centroid = feature_geom.Centroid()

                utm_code = (
                    math.floor((feature_centroid.GetX() + 180)/6) % 60) + 1
                lat_code = 6 if feature_centroid.GetY() > 0 else 7
                epsg_code = int('32%d%02d' % (lat_code, utm_code))
                epsg_sr = osr.SpatialReference()
                epsg_sr.ImportFromEPSG(epsg_code)

                local_feature_geom = feature_geom.Clone()
                coord_trans = osr.CoordinateTransformation(wgs84_sr, epsg_sr)
                local_feature_geom.Transform(coord_trans)

                # m^2 to Ha
                feature_area_ha = local_feature_geom.GetArea() * 0.0001

                feature.SetField('area_ha', feature_area_ha)

                for scenario in SCENARIO_LIST:
                    scenaro_cursor.execute(
                        """SELECT total_export FROM nutrient_export
                        WHERE (ws_prefix_key = ? and scenario_key = ?)""",
                        (ws_id, scenario))
                    result = scenaro_cursor.fetchone()
                    if result is not None:
                        feature.SetField('%snex' % scenario, result[0])
                        feature.SetField(
                            '%snexha' % scenario, result[0] / feature_area_ha)

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
    zero_mask = (d_dn_array == 0) | (d_up_array == 0)
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


def main(raw_iam_token_path, raw_workspace_dir):
    """Entry point."""
    iam_token_path = os.path.normpath(raw_iam_token_path)
    workspace_dir = os.path.normpath(raw_workspace_dir)
    downloads_dir = os.path.join(workspace_dir, BUCKET_DOWNLOAD_DIR)
    churn_dir = os.path.join(workspace_dir, CHURN_DIR)
    watershed_processing_dir = os.path.join(
        workspace_dir, WATERSHED_PROCESSING_DIR)

    try:
        os.makedirs(workspace_dir)
    except OSError:
        pass

    task_graph = taskgraph.TaskGraph(
        os.path.join(workspace_dir, 'taskgraph_cache'), N_CPUS,
        TASKGRAPH_REPORTING_FREQUENCY, TASKGRAPH_DELAYED_START)

    ag_load_scenarios_archive_path = os.path.join(
        downloads_dir, 'ag_load_scenarios_blake2b_2c8661957382df98041890e20ede8c93.zip')
    fetch_ag_load_scenarios_task = task_graph.add_task(
        func=reproduce.utils.google_bucket_fetch_and_validate,
        args=(
            'ipbes-ndr-ecoshard-data',
            'ag_load_scenarios_blake2b_2c8661957382df98041890e20ede8c93.zip',
            iam_token_path,
            ag_load_scenarios_archive_path),
        target_path_list=[ag_load_scenarios_archive_path],
        task_name='fetch ag load scenarios')
    ag_load_scenarios_touch_file_path = (
        os.path.join(
            churn_dir, os.path.basename(ag_load_scenarios_archive_path) + '_unzipped'))
    unzip_ag_load_scenarios_task = task_graph.add_task(
        func=unzip_file,
        args=(
            ag_load_scenarios_archive_path, churn_dir,
            ag_load_scenarios_touch_file_path),
        target_path_list=[ag_load_scenarios_touch_file_path],
        dependent_task_list=[fetch_ag_load_scenarios_task],
        task_name=f'unzip ag_load_scenarios')
    ag_load_scenarios_dir_path = os.path.join(
        churn_dir, AG_LOAD_DIR)

    precip_scenarios_archive_path = os.path.join(
        downloads_dir, 'precip_scenarios_for_ndr_blake2b_393c496d9c2a14e47136d51522eea975.zip')
    fetch_precip_scenarios_task = task_graph.add_task(
        func=reproduce.utils.google_bucket_fetch_and_validate,
        args=(
            'ipbes-ndr-ecoshard-data',
            'precip_scenarios_for_ndr_blake2b_393c496d9c2a14e47136d51522eea975.zip',
            iam_token_path,
            precip_scenarios_archive_path),
        target_path_list=[precip_scenarios_archive_path],
        task_name='fetch precip scenarios')
    precip_scenarios_touch_file_path = (
        os.path.join(
            churn_dir, os.path.basename(precip_scenarios_archive_path) + '_unzipped'))
    precip_scenarios_dir_path = os.path.join(
        churn_dir, PRECIP_DIR)
    unzip_precip_scenarios_task = task_graph.add_task(
        func=unzip_file,
        args=(
            precip_scenarios_archive_path, precip_scenarios_dir_path,
            precip_scenarios_touch_file_path),
        target_path_list=[precip_scenarios_touch_file_path],
        dependent_task_list=[fetch_precip_scenarios_task],
        task_name=f'unzip precip_scenarios')

    globio_landuse_archive_path = os.path.join(
        downloads_dir,
        'globio_landuse_historic_and_ssp_blake2b_4153935fd8cbb510d8500d59272e4479.zip')
    fetch_globio_landuse_task = task_graph.add_task(
        func=reproduce.utils.google_bucket_fetch_and_validate,
        args=(
            'ipbes-ndr-ecoshard-data',
            'globio_landuse_historic_and_ssp_blake2b_4153935fd8cbb510d8500d59272e4479.zip',
            iam_token_path, globio_landuse_archive_path),
        target_path_list=[globio_landuse_archive_path],
        task_name='fetch globio landuse')
    globio_landuse_touch_file_path = (
        os.path.join(
            churn_dir, os.path.basename(globio_landuse_archive_path) + '_unzipped'))

    globio_landuse_dir_path = os.path.join(
        churn_dir, LANDUSE_DIR)
    unzip_globio_landuse_task = task_graph.add_task(
        func=unzip_file,
        args=(
            globio_landuse_archive_path, globio_landuse_dir_path,
            globio_landuse_touch_file_path),
        target_path_list=[globio_landuse_touch_file_path],
        dependent_task_list=[fetch_globio_landuse_task],
        task_name=f'unzip globio landuse scenarios')

    watersheds_archive_path = os.path.join(
        downloads_dir,
        'watersheds_globe_HydroSHEDS_15arcseconds_blake2b_14ac9c77d2076d51b0258fd94d9378d4.zip')
    fetch_watersheds_task = task_graph.add_task(
        func=reproduce.utils.google_bucket_fetch_and_validate,
        args=(
            'ipbes-ndr-ecoshard-data', 'watersheds_globe_HydroSHEDS_15arcseconds_blake2b_14ac9c77d2076d51b0258fd94d9378d4.zip',
            iam_token_path, watersheds_archive_path),
        target_path_list=[watersheds_archive_path],
        task_name='download watersheds')
    watersheds_touch_file_path = os.path.join(
        churn_dir, os.path.basename(watersheds_archive_path) + '_unzipped')
    unzip_watersheds_task = task_graph.add_task(
        func=unzip_file,
        args=(
            watersheds_archive_path, churn_dir,
            watersheds_touch_file_path),
        target_path_list=[watersheds_touch_file_path],
        dependent_task_list=[fetch_watersheds_task],
        task_name=f'unzip watersheds_globe_HydroSHEDS_15arcseconds')
    # this will be where all the watersheds unzip
    watersheds_dir_path = os.path.join(
        churn_dir, 'watersheds_globe_HydroSHEDS_15arcseconds')

    biophysical_table_path = os.path.join(
        downloads_dir, 'NDR_representative_table_md5_3b5aeb8dba0f615aed858ad98bbfb755.csv')
    download_biophysical_table_task = task_graph.add_task(
        func=reproduce.utils.google_bucket_fetch_and_validate,
        args=(
            'ipbes-ndr-ecoshard-data',
            'NDR_representative_table_md5_3b5aeb8dba0f615aed858ad98bbfb755.csv',
            iam_token_path, biophysical_table_path),
        target_path_list=[biophysical_table_path],
        task_name='download biophysical table')
    clean_biophysical_table_pickle_path = os.path.join(
        churn_dir, 'biophysical.pickle')
    clean_biophysical_task = task_graph.add_task(
        func=clean_and_pickle_biophysical_table,
        args=(biophysical_table_path, clean_biophysical_table_pickle_path),
        target_path_list=[clean_biophysical_table_pickle_path],
        dependent_task_list=[download_biophysical_table_task])

    global_dem_archive_path = os.path.join(
        downloads_dir, 'global_dem_3s_blake2b_0532bf0a1bedbe5a98d1dc449a33ef0c.zip')
    global_dem_download_task = task_graph.add_task(
        func=reproduce.utils.google_bucket_fetch_and_validate,
        args=(
            'ipbes-ndr-ecoshard-data',
            'global_dem_3s_blake2b_0532bf0a1bedbe5a98d1dc449a33ef0c.zip',
            iam_token_path, global_dem_archive_path),
        target_path_list=[global_dem_archive_path],
        task_name='download dem archive')
    dem_touch_file_path = os.path.join(
        churn_dir, os.path.basename(global_dem_archive_path) + '_unzipped')
    unzip_dem_task = task_graph.add_task(
        func=unzip_file,
        args=(
            global_dem_archive_path, churn_dir,
            dem_touch_file_path),
        target_path_list=[dem_touch_file_path],
        dependent_task_list=[global_dem_download_task],
        task_name=f'unzip global_dem')
    dem_path_dir = os.path.join(churn_dir, 'global_dem_3s')

    dem_rtree_path = os.path.join(churn_dir, RTREE_PATH)
    dem_path_index_map_path = os.path.join(
        churn_dir, 'dem_rtree_path_index_map.dat')
    build_dem_rtree_task = task_graph.add_task(
        func=build_raster_rtree,
        args=(dem_path_dir, dem_path_index_map_path, dem_rtree_path),
        target_path_list=[
            dem_rtree_path+'.dat',  # rtree adds a ".dat" file
            dem_path_index_map_path],
        dependent_task_list=[unzip_dem_task],
        task_name='build_raster_rtree')

    # create a results database
    database_path = os.path.join(workspace_dir, 'ipbes_ndr_results.db')
    sql_create_projects_table = (
        """
        CREATE TABLE IF NOT EXISTS nutrient_export (
            ws_prefix_key TEXT NOT NULL,
            scenario_key TEXT NOT NULL,
            total_export REAL NOT NULL,
            PRIMARY KEY (ws_prefix_key, scenario_key)
        );
        CREATE UNIQUE INDEX IF NOT EXISTS ws_scenario_index
        ON nutrient_export (ws_prefix_key, scenario_key);
        CREATE INDEX IF NOT EXISTS ws_index
        ON nutrient_export (ws_prefix_key);

        CREATE TABLE IF NOT EXISTS geometry_table (
            ws_prefix_key TEXT NOT NULL,
            geometry_wkb BLOB NOT NULL,
            PRIMARY KEY (geometry_id)
        );
        CREATE UNIQUE INDEX IF NOT EXISTS geometry_id_index
        ON geometry_table (geometry_id);
        """)

    conn = sqlite3.connect(database_path)
    cursor = conn.cursor()
    cursor.executescript(sql_create_projects_table)

    unzip_watersheds_task.join()
    build_dem_rtree_task.join()

    biophysical_table = pandas.read_csv(biophysical_table_path)

    eff_n_lucode_map = dict(
        zip(biophysical_table['ID'], biophysical_table['eff_n']))
    load_n_lucode_map = dict(
        zip(biophysical_table['ID'], biophysical_table['load_n']))

    LOGGER.info("scheduling watershed processing")
    global_watershed_path_list = glob.glob(
        os.path.join(watersheds_dir_path, '*.shp'))
    task_id = 0
    for global_watershed_path in global_watershed_path_list:
        watershed_basename = os.path.splitext(
            os.path.basename(global_watershed_path))[0]
        watershed_vector = gdal.OpenEx(global_watershed_path, gdal.OF_VECTOR)
        watershed_layer = watershed_vector.GetLayer()
        for watershed_feature in watershed_layer:
            watershed_fid = watershed_feature.GetFID()
            ws_prefix = 'ws_%s_%d' % (watershed_basename, watershed_fid)
            watershed_geom = watershed_feature.GetGeometryRef()
            watershed_area = watershed_geom.GetArea()
            watershed_geom = None
            if watershed_area < 0.03:
                #  0.03 square degrees is a healthy underapproximation of
                # 100 sq km which is about the minimum watershed size we'd
                # want.
                continue
            schedule_watershed_processing(
                task_graph, task_id, ws_prefix, watershed_fid,
                watershed_feature, dem_rtree_path,
                dem_path_index_map_path,
                eff_n_lucode_map,
                load_n_lucode_map,
                churn_dir, database_path,
                watershed_processing_dir)
            watershed_feature = None
            task_id -= 1
            break
        watershed_layer = None
        watershed_vector = None
        break

    task_graph.close()
    task_graph.join()
    LOGGER.info("all done :)")


def schedule_watershed_processing(
        task_graph, task_id, ws_prefix, watershed_fid,
        base_watershed_feature, dem_rtree_path, dem_path_index_map_path,
        eff_n_lucode_map, load_n_lucode_map, root_data_dir,
        target_result_database_path, workspace_dir):
    """Process a watershed for NDR analysis.

    A successful call to this function will insert any new geometry
        into `geometry_table` and any new export values into
        `nutrient_export` in the `target_result_database_path` database.

    Parameters:
        task_graph (TaskGraph): taskgraph scheduler to schedule against.
        task_id (int): priority to set taskgraph at.
        watershed_fid (integer): FID of the watershed to schedule.
        base_watershed_layer (ogr.Layer): base watershed layer.
        dem_rtree_path (str): path to RTree that can be used to determine
            which DEM tiles intersect a bounding box.
        dem_path_index_map_path (str): path to pickled dictionary that maps
            `dem_rtree_path` IDs to file paths of the DEM tile that matches
            that id.
        eff_n_lucode_map (dict): maps lucodes to NDR efficiency values.
        load_n_lucode_map (dict): maps lucodes to NDR load values.
        root_data_dir (str): path to directory containing all landcover,
            precip, and load rasters referenced relatively in
            LANDCOVER_RASTER_PATHS, PRECIP_RASTER_PATHS, and AG_RASTER_PATHS.
        target_result_database_path (str): A database that has two tables
            'nutrient_export (ws_prefix_key, scenario_key, total_export)'
            'geometry_table (ws_prefix_key, geometry_wkb)'
        workspace_dir (str): path to workspace to create working files.

    Returns:
        None.
    """
    # make a few subdirectories so we don't explode on number of files per
    # directory. The largest watershed is 726k
    last_digits = '%.4d' % watershed_fid
    ws_working_dir = os.path.join(
        workspace_dir, last_digits[-1], last_digits[-2],
        last_digits[-3], last_digits[-4],
        "%s_working_dir" % ws_prefix)

    watershed_dem_path = os.path.join(
        ws_working_dir, 'ws_%s_dem.tif' % ws_prefix)

    watershed_geometry = base_watershed_feature.GetGeometryRef()
    watershed_bb = [
        watershed_geometry.GetEnvelope()[i] for i in [0, 2, 1, 3]]

    merge_watershed_dems_task = task_graph.add_task(
        func=merge_watershed_dems,
        args=(
            watershed_bb, watershed_fid, dem_rtree_path,
            dem_path_index_map_path, watershed_dem_path),
        target_path_list=[watershed_dem_path],
        task_name='merge_watershed_dems_%s' % ws_prefix,
        priority=task_id)

    masked_watershed_dem_path = watershed_dem_path.replace(
        '.tif', '_masked.tif')

    centroid_geom = watershed_geometry.Centroid()
    utm_code = (math.floor((centroid_geom.GetX() + 180)/6) % 60) + 1
    lat_code = 6 if centroid_geom.GetY() > 0 else 7
    epsg_code = int('32%d%02d' % (lat_code, utm_code))
    epsg_srs = osr.SpatialReference()
    epsg_srs.ImportFromEPSG(epsg_code)
    utm_pixel_size = 90.0

    local_watershed_path = os.path.join(ws_working_dir, '%s.gpkg' % ws_prefix)
    watershed_geom_wkb = watershed_geometry.ExportToWkb()

    reproject_watershed_task = task_graph.add_task(
        func=reproject_geometry_to_target,
        args=(
            watershed_geom_wkb,
            watershed_geometry.GetSpatialReference().ExportToWkt(),
            epsg_srs.ExportToWkt(), local_watershed_path),
        target_path_list=[local_watershed_path],
        task_name='project_watershed_%s' % ws_prefix,
        priority=task_id)

    mask_watershed_dem_task = task_graph.add_task(
        func=mask_raster_by_vector,
        args=(
            watershed_dem_path, local_watershed_path,
            masked_watershed_dem_path),
        target_path_list=[masked_watershed_dem_path],
        dependent_task_list=[
            reproject_watershed_task, merge_watershed_dems_task],
        task_name='mask dem %s' % ws_prefix)

    base_raster_path_list = (
        [masked_watershed_dem_path] +
        [os.path.join(root_data_dir, path)
         for path in list(LANDCOVER_RASTER_PATHS.values()) +
         list(PRECIP_RASTER_PATHS.values()) + list(AG_RASTER_PATHS.values())])

    def _base_to_aligned_path_op(base_path):
        """Convert global raster path to local."""
        return os.path.join(
            ws_working_dir, '%s_%s_aligned.tif' % (
                ws_prefix, os.path.splitext(os.path.basename(base_path))[0]))

    aligned_path_list = [
        _base_to_aligned_path_op(path) for path in base_raster_path_list]

    wgs84_sr = osr.SpatialReference()
    wgs84_sr.ImportFromEPSG(4326)
    target_bounding_box = pygeoprocessing.transform_bounding_box(
        watershed_bb, wgs84_sr.ExportToWkt(),
        epsg_srs.ExportToWkt())

    # clip dem, precip, & landcover to size of DEM? use 'mode'
    # we know the input rasters are WGS84 unprojected

    align_resize_task = task_graph.add_task(
        func=pygeoprocessing.align_and_resize_raster_stack,
        args=(
            base_raster_path_list, aligned_path_list,
            ['near'] * len(base_raster_path_list),
            (utm_pixel_size, -utm_pixel_size),
            target_bounding_box),
        kwargs={
            'base_sr_wkt_list': [wgs84_sr.ExportToWkt()] * len(
                base_raster_path_list),
            'target_sr_wkt': epsg_srs.ExportToWkt()
            },
        target_path_list=aligned_path_list,
        dependent_task_list=[mask_watershed_dem_task],
        task_name='align resize %s' % ws_prefix,
        priority=task_id)

    # fill and route dem
    filled_watershed_dem_path = os.path.join(
        ws_working_dir, '%s_dem_filled.tif' % ws_prefix)
    flow_dir_path = os.path.join(
        ws_working_dir, '%s_flow_dir.tif' % ws_prefix)

    fill_pits_task = task_graph.add_task(
        func=pygeoprocessing.routing.fill_pits,
        args=(
            (aligned_path_list[0], 1),
            filled_watershed_dem_path),
        kwargs={'working_dir': ws_working_dir},
        target_path_list=[
            filled_watershed_dem_path],
        dependent_task_list=[align_resize_task],
        task_name='fill pits %s' % ws_prefix,
        priority=task_id)

    flow_dir_task = task_graph.add_task(
        func=pygeoprocessing.routing.flow_dir_mfd,
        args=(
            (filled_watershed_dem_path, 1),
            flow_dir_path),
        kwargs={'working_dir': ws_working_dir},
        target_path_list=[
            flow_dir_path],
        dependent_task_list=[fill_pits_task],
        task_name='flow dir %s' % ws_prefix,
        priority=task_id)

    # flow accum dem
    flow_accum_path = os.path.join(
        ws_working_dir, '%s_flow_accum.tif' % ws_prefix)
    flow_accum_task = task_graph.add_task(
        func=pygeoprocessing.routing.flow_accumulation_mfd,
        args=(
            (flow_dir_path, 1), flow_accum_path),
        target_path_list=[flow_accum_path],
        dependent_task_list=[flow_dir_task],
        task_name='flow accmulation %s' % ws_prefix,
        priority=task_id)

    # calculate slope
    slope_raster_path = os.path.join(
        ws_working_dir, '%s_slope.tif' % ws_prefix)
    calculate_slope_task = task_graph.add_task(
        func=pygeoprocessing.calculate_slope,
        args=((filled_watershed_dem_path, 1), slope_raster_path),
        target_path_list=[slope_raster_path],
        dependent_task_list=[fill_pits_task],
        task_name='calculate_slope_%s' % ws_prefix,
        priority=task_id)

    clamp_slope_raster_path = os.path.join(
        ws_working_dir, '%s_clamp_slope.tif' % ws_prefix)
    clamp_slope_task = task_graph.add_task(
        func=ClampOp(
            (slope_raster_path, 1), 0.005, clamp_slope_raster_path),
        target_path_list=[clamp_slope_raster_path],
        dependent_task_list=[calculate_slope_task],
        task_name='clamp_slope_%s' % ws_prefix,
        priority=task_id)

    # calculate D_up
    slope_accum_watershed_dem_path = os.path.join(
        ws_working_dir, '%s_s_accum.tif' % ws_prefix)
    slope_accmulation_task = task_graph.add_task(
        func=pygeoprocessing.routing.flow_accumulation_mfd,
        args=(
            (flow_dir_path, 1), slope_accum_watershed_dem_path),
        kwargs={
            'weight_raster_path_band': (clamp_slope_raster_path, 1)},
        target_path_list=[slope_accum_watershed_dem_path],
        dependent_task_list=[fill_pits_task, clamp_slope_task],
        task_name='slope_accmulation_%s' % ws_prefix,
        priority=task_id)

    d_up_raster_path = os.path.join(ws_working_dir, '%s_d_up.tif' % ws_prefix)
    d_up_task = task_graph.add_task(
        func=DUpOp(
            utm_pixel_size**2, slope_accum_watershed_dem_path,
            flow_accum_path, d_up_raster_path),
        target_path_list=[d_up_raster_path],
        dependent_task_list=[slope_accmulation_task, flow_accum_task],
        task_name='d_up_%s' % ws_prefix,
        priority=task_id)

    # calculate the flow channels
    channel_path = os.path.join(ws_working_dir, '%s_channel.tif' % ws_prefix)
    threshold_flow_task = task_graph.add_task(
        func=threshold_flow_accumulation,
        args=(
            flow_accum_path, FLOW_THRESHOLD, channel_path),
        target_path_list=[channel_path],
        dependent_task_list=[flow_accum_task],
        task_name='threshold flow accum %s' % ws_prefix,
        priority=task_id)

    # calculate flow path in pixels length down to stream
    pixel_flow_length_raster_path = os.path.join(
        ws_working_dir, '%s_pixel_flow_length.tif' % ws_prefix)
    downstream_flow_length_task = task_graph.add_task(
        func=pygeoprocessing.routing.distance_to_channel_mfd,
        args=(
            (flow_dir_path, 1), (channel_path, 1),
            pixel_flow_length_raster_path),
        target_path_list=[pixel_flow_length_raster_path],
        dependent_task_list=[fill_pits_task, threshold_flow_task],
        task_name='downstream_pixel_flow_length_%s' % ws_prefix,
        priority=task_id)

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
        priority=task_id)

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
        priority=task_id)

    # calculate D_dn: downstream sum of distance / downstream slope
    d_dn_raster_path = os.path.join(
        ws_working_dir, '%s_d_dn.tif' % ws_prefix)
    d_dn_task = task_graph.add_task(
        func=pygeoprocessing.routing.distance_to_channel_mfd,
        args=(
            (flow_dir_path, 1), (channel_path, 1), d_dn_raster_path),
        kwargs={
            'weight_raster_path_band': (d_dn_per_pixel_path, 1)
            },
        target_path_list=[d_dn_raster_path],
        dependent_task_list=[
            fill_pits_task, flow_accum_task, d_dn_per_pixel_task],
        task_name='d_dn_%s' % ws_prefix,
        priority=task_id)

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
        priority=task_id)

    for landcover_id, global_landcover_path in LANDCOVER_RASTER_PATHS.items():
        local_landcover_path = os.path.join(
            ws_working_dir, '%s_%s_aligned.tif' % (
                ws_prefix, os.path.splitext(os.path.basename(
                    global_landcover_path))[0]))
        eff_n_raster_path = os.path.join(
            ws_working_dir, '%s_eff_n.tif' % ws_prefix)
        reclassify_eff_n_task = task_graph.add_task(
            func=pygeoprocessing.reclassify_raster,
            args=(
                (local_landcover_path, 1), eff_n_lucode_map,
                eff_n_raster_path, gdal.GDT_Float32, NODATA),
            target_path_list=[eff_n_raster_path],
            dependent_task_list=[align_resize_task],
            task_name='reclassify_eff_n_%s' % ws_prefix,
            priority=task_id)

        load_n_raster_path = os.path.join(
            ws_working_dir, '%s_load_n.tif' % ws_prefix)
        reclassify_load_n_task = task_graph.add_task(
            func=pygeoprocessing.reclassify_raster,
            args=(
                (local_landcover_path, 1), load_n_lucode_map,
                load_n_raster_path, gdal.GDT_Float32, NODATA),
            target_path_list=[load_n_raster_path],
            dependent_task_list=[align_resize_task],
            task_name='reclasify_load_n_%s' % ws_prefix,
            priority=task_id)

    LOGGER.warn("don't forget the rest!")
    return

    # reclassify eff_n
    for x in landcover_raster_path_lists:
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


def merge_watershed_dems(
        watershed_bb, watershed_id, dem_rtree_path, dem_path_index_map_path,
        target_dem_path):
    """Find DEMs that overlap the given watershed polyon by id.

    Parameters:
        watershed_bb (string): watershed bounding box
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
    os.makedirs(os.path.dirname(target_dem_path), exist_ok=True)
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
            "no overlapping dems found for %s wsid %d", target_dem_path,
            watershed_id)


def build_raster_rtree(
        raster_dir_path, raster_path_index_map_path, raster_rtree_path):
    """Build RTree for list of rasters if RTree does not exist.

    If the RTree already exists, this function logs a warning and returns.

    Paramters:
        raster_dir_path (str): path to a directory of GIS rasters.
        raster_path_index_map_path (str): this is a path to a pickle file
            generated by this function that will index RTree indexes to
            the rasters on disk that will intersect the RTree.
        raster_rtree_path (str): this is the target path to the saved rTree
            generated by this call.

    Returns:
        None.
    """
    LOGGER.info('building rTree %s', raster_rtree_path+'.dat')
    if os.path.exists(raster_rtree_path+'.dat'):
        LOGGER.warn('%s exists so skipping rTree creation.', raster_rtree_path)
        return
    raster_rtree = rtree.index.Index(raster_rtree_path)
    raster_path_index_map = {}
    raster_path_list = glob.glob(os.path.join(raster_dir_path, '*.tif'))
    for raster_id, raster_path in enumerate(raster_path_list):
        raster_info = pygeoprocessing.get_raster_info(raster_path)
        raster_path_index_map[raster_id] = raster_path
        raster_rtree.insert(raster_id, raster_info['bounding_box'])
    raster_rtree.close()
    with open(raster_path_index_map_path, 'wb') as f:
        dill.dump(raster_path_index_map, f)


def reproject_geometry_to_target(
        geom_wkb, base_sr_wkt, target_sr_wkt, target_path):
    """Reproject a single OGR DataSource feature.

    Transforms the features of the base vector to the desired output
    projection in a new ESRI Shapefile.

    Parameters:
        geom_wkt (str): base geometry in wkb.
        base_sr_wkt (str): the spatial reference in wkt for `geom_wkb`.
        target_sr_wkt (str): the desired output projection in Well Known Text
            (by layer.GetSpatialRef().ExportToWkt())
        feature_id (int): the feature to reproject and copy.
        target_path (str): the filepath to the transformed shapefile

    Returns:
        None
    """
    os.makedirs(os.path.dirname(target_path), exist_ok=True)
    # if this file already exists, then remove it
    if os.path.isfile(target_path):
        LOGGER.warn(
            "reproject_vector: %s already exists, removing and overwriting",
            target_path)
        os.remove(target_path)

    target_sr = osr.SpatialReference(target_sr_wkt)

    # create a new shapefile from the orginal_datasource
    target_driver = gdal.GetDriverByName('GPKG')
    target_vector = target_driver.Create(
        target_path, 0, 0, 0, gdal.GDT_Unknown)
    layer_name = os.path.splitext(os.path.basename(target_path))[0]
    base_geom = ogr.CreateGeometryFromWkb(geom_wkb)
    target_layer = target_vector.CreateLayer(
        layer_name, target_sr, base_geom.GetGeometryType())

    # Create a coordinate transformation
    base_sr = osr.SpatialReference(base_sr_wkt)
    coord_trans = osr.CoordinateTransformation(base_sr, target_sr)

    # Transform geometry into format desired for the new projection
    error_code = base_geom.Transform(coord_trans)
    if error_code != 0:  # error
        # this could be caused by an out of range transformation
        # whatever the case, don't put the transformed poly into the
        # output set
        raise ValueError(
            "Unable to reproject geometry on %s." % target_path)

    # Copy original_datasource's feature and set as new shapes feature
    target_feature = ogr.Feature(target_layer.GetLayerDefn())
    target_feature.SetGeometry(base_geom)
    target_layer.CreateFeature(target_feature)
    target_feature = None


def unzip_file(zipfile_path, target_dir, touchfile_path):
    """Unzip contents of `zipfile_path`.

    Parameters:
        zipfile_path (string): path to a zipped file.
        target_dir (string): path to extract zip file to.
        touchfile_path (string): path to a file to create if unzipping is
            successful.

    Returns:
        None.

    """
    with zipfile.ZipFile(zipfile_path, 'r') as zip_ref:
        zip_ref.extractall(target_dir)

    with open(touchfile_path, 'w') as touchfile:
        touchfile.write(f'unzipped {zipfile_path}')


def clean_and_pickle_biophysical_table(
        biophysical_table_path, clean_biophysical_table_pickle_path):
    """Clean out nans and set replacement lucodde and pickle table."""

    biophysical_table = pandas.read_csv(biophysical_table_path)
    # clean up biophysical table
    biophysical_table = biophysical_table.fillna(0)
    biophysical_table.ix[
        biophysical_table['load_n'] == 'use raster', 'load_n'] = (
            USE_AG_LOAD_ID)
    biophysical_table['load_n'] = biophysical_table['load_n'].apply(
        pandas.to_numeric)

    with open(clean_biophysical_table_pickle_path, 'wb') as clean_bpt_file:
        dill.dump(biophysical_table, clean_bpt_file)


def mask_raster_by_vector(
        base_raster_path, vector_path, target_masked_raster_path):
    """Mask out values in base raster that don't overlap with vector_path."""
    base_raster_info = pygeoprocessing.get_raster_info(base_raster_path)
    pygeoprocessing.new_raster_from_base(
        base_raster_path, target_masked_raster_path,
        base_raster_info['datatype'], base_raster_info['nodata'])

    pygeoprocessing.rasterize(
        vector_path, target_masked_raster_path, [1], None)

    target_raster = gdal.OpenEx(
        target_masked_raster_path, gdal.OF_RASTER | gdal.GA_Update)
    target_band = target_raster.GetRasterBand(1)
    base_raster = gdal.OpenEx(
        base_raster_path, gdal.OF_RASTER)
    base_band = base_raster.GetRasterBand(1)

    for offset_dict in pygeoprocessing.iterblocks(
            base_raster_path, offset_only=True):
        target_array = target_band.ReadAsArray(**offset_dict)
        mask_array = numpy.isclose(target_array, 1)
        base_array = base_band.ReadAsArray(**offset_dict)
        target_array[mask_array] = base_array[mask_array]
        target_band.WriteArray(
            target_array, xoff=offset_dict['xoff'],
            yoff=offset_dict['yoff'])


def threshold_flow_accumulation(
        flow_accum_path, flow_threshold, target_channel_path):
    """Calculate channel raster by thresholding flow accumulation.

    Parameters:
        flow_accum_path (str): path to a single band flow accumulation raster.
        flow_threshold (float): if the value in `flow_accum_path` is less
            than or equal to this value, the pixel will be classified as a
            channel.
        target_channel_path (str): path to target raster that will contain
            pixels set to 1 if they are a channel, 0 if not, and possibly
            between 0 and 1 if a partial channel. (to be defined).

    Returns:
        None.
    """
    nodata = pygeoprocessing.get_raster_info(flow_accum_path)['nodata'][0]
    channel_nodata = -1.0

    def threshold_op(flow_val):
        valid_mask = ~numpy.isclose(flow_val, nodata)
        result = numpy.empty(flow_val.shape, dtype=numpy.float32)
        result[:] = channel_nodata
        result[valid_mask] = flow_val[valid_mask] >= flow_threshold
        return result

    pygeoprocessing.raster_calculator(
        [(flow_accum_path, 1)], threshold_op, target_channel_path,
        gdal.GDT_Float32, channel_nodata)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        LOGGER.error(
            "usage: python %s iam_token_path workspace_dir", sys.argv[0])
        sys.exit(-1)
    raw_iam_token_path = sys.argv[1]
    raw_workspace_dir = sys.argv[2]
    if not os.path.isfile(raw_iam_token_path):
        LOGGER.error(
            '%s is not a file, should be an iam token', raw_workspace_dir)
        sys.exit(-1)
    if os.path.isfile(raw_workspace_dir):
        LOGGER.error(
            '%s is supposed to be the workspace directory but points to an '
            'existing file' % raw_workspace_dir)
        sys.exit(-1)
    main(raw_iam_token_path, raw_workspace_dir)
