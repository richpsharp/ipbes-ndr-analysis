"""Script to manage NDR runs for IPBES project."""
import logging
import os
import glob
import sqlite3

import rtree.index
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


def result_in_database(database_path, gridcode, ws_id):
    """True if ws_prefix in database."""
    conn = sqlite3.connect(database_path)
    if conn is not None:
        try:
            cursor = conn.cursor()
            print gridcode, ws_id
            cursor.execute(
                """SELECT cur_export FROM nutrient_export
                WHERE (GRIDCODE = ? and WS_ID = ?)""", (
                    gridcode, ws_id))
            result = cursor.fetchone()
            if result is None:
                return False
            return True
        except sqlite3.OperationalError:
            LOGGER.exception("operational error on %s"% ws_prefix)
            return False
    return False


def build_watershed_rtree(
        watershed_path_list, watershed_path_index_map_path):
    """Build RTree indexed by FID for points in `wwwiii_vector_path`."""
    LOGGER.info('building rTree %s', watershed_path_index_map_path + '.dat')
    if os.path.exists(watershed_path_index_map_path+'.dat'):
        LOGGER.warn(
            '%s exists so skipping creation.', watershed_path_index_map_path)
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
            ws_prefix = 'ws_%s_%.4d' % (
                watershed_basename, watershed_id)

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


def analyze_grid(
        grid_fid, watershed_path_index_map_path, database_path):
    conn = sqlite3.connect(database_path)
    cursor = conn.cursor()

    watershed_rtree = rtree.index.Index(watershed_path_index_map_path)

    grid_vector = gdal.OpenEx(DEGREE_GRID_PATH, gdal.OF_VECTOR)
    grid_layer = grid_vector.GetLayer()

    wgs84_srs = osr.SpatialReference()
    wgs84_srs.ImportFromEPSG(4326)
    esri_driver = gdal.GetDriverByName('ESRI Shapefile')

    grid_feature = grid_layer.GetFeature(grid_fid)
    if not grid_feature:
        return
    grid_code = grid_feature.GetField('GRIDCODE')
    print grid_code
    grid_geometry = grid_feature.GetGeometryRef()
    grid_bounds = grid_geometry.GetEnvelope()
    results = list(watershed_rtree.intersection(
        (grid_bounds[0], grid_bounds[2], grid_bounds[1], grid_bounds[3]),
        objects=True))
    if results:
        for watershed_id in [str(x.object) for x in results]:
            if result_in_database(database_path, grid_code, watershed_id):
                LOGGER.debug("%s %s in database", grid_code, watershed_id)
                continue
            # this truncates zeros
            shp_id = '%s%s' % (watershed_id[0:-4], int(watershed_id[-4:]))
            watershed_path = os.path.join(
                'ndr_workspace', '/'.join(reversed(watershed_id[-4:])),
                '%s_working_dir' % watershed_id, '%s.shp' % shp_id)
            if os.path.exists(watershed_path):
                watershed_vector = gdal.OpenEx(
                    watershed_path, gdal.OF_VECTOR)
                watershed_layer = watershed_vector.GetLayer()
                watershed_feature = watershed_layer.GetNextFeature()
                watershed_geometry = watershed_feature.GetGeometryRef()

                watershed_srs = (
                    watershed_geometry.GetSpatialReference())
                utm_to_wgs84 = osr.CoordinateTransformation(
                    watershed_srs, wgs84_srs)
                wgs84_to_utm = osr.CoordinateTransformation(
                    wgs84_srs, watershed_srs)
                watershed_geometry.Transform(utm_to_wgs84)
                watershed_intersect_geom = (
                    watershed_geometry.Intersection(grid_geometry))
                fraction_covered = watershed_intersect_geom.GetArea() / (
                    grid_geometry.GetArea())
                if not watershed_intersect_geom.IsEmpty():
                    local_clip_path = os.path.join(
                        os.path.dirname(watershed_path),
                        'grid_clipped%s.shp' % grid_code)
                    if os.path.exists(local_clip_path):
                        os.remove(local_clip_path)
                    watershed_clip_vector = esri_driver.CreateCopy(
                        local_clip_path, watershed_vector)
                    watershed_clip_layer = watershed_clip_vector.GetLayer()
                    watershed_clip_feature = (
                        watershed_clip_layer.GetNextFeature())
                    watershed_intersect_geom.Transform(wgs84_to_utm)
                    watershed_clip_feature.SetGeometry(
                        watershed_intersect_geom)

                    watershed_clip_layer.SetFeature(
                        watershed_clip_feature)
                    watershed_clip_layer.SyncToDisk()
                    watershed_clip_vector.FlushCache()
                    watershed_intersect_geom = None
                    watershed_clip_feature = None
                    watershed_clip_layer = None
                    watershed_clip_vector = None

                    export_values = {}
                    for scenario_id in ['cur', 'ssp1', 'ssp3', 'ssp5']:
                        export_path = os.path.join(
                            os.path.dirname(watershed_path),
                            '%s_%s_n_export.tif' % (
                                watershed_id, scenario_id))
                        export_stats = pygeoprocessing.zonal_statistics(
                            (export_path, 1), local_clip_path, 'BASIN_ID')
                        if export_stats:
                            export_values[scenario_id] = (
                                export_stats.itervalues().next()['sum'])
                        else:
                            break
                    if not export_values:
                        continue
                    try:
                        cursor.execute(
                            """INSERT INTO nutrient_export VALUES
                             (?, ?, ?, ?, ?, ?, ?)""", (
                                 grid_code, watershed_id,
                                 export_values['cur'],
                                 export_values['ssp1'],
                                 export_values['ssp3'],
                                 export_values['ssp5'],
                                 fraction_covered))
                        conn.commit()
                    except:
                        LOGGER.exception('"%s"', shp_id)
                watershed_intersect_geom = None
                watershed_clip_feature = None
                watershed_clip_layer = None
                watershed_clip_vector = None
                watershed_vector = None
    conn.close()


def main():
    """Entry point."""
    try:
        os.makedirs(TARGET_WORKSPACE)
    except OSError:
        pass

    database_path = os.path.join(
        TARGET_WORKSPACE, 'ndr_to_degree.db')

    sql_create_projects_table = (
        """ CREATE TABLE IF NOT EXISTS nutrient_export (
            GRIDCODE TEXT NOT NULL,
            WS_ID TEXT NOT NULL,
            cur_export REAL NOT NULL,
            ssp1_export REAL NOT NULL,
            ssp3_export REAL NOT NULL,
            ssp5_export REAL NOT NULL,
            fraction_covered REAL NOT NULL
        ); """)

    # create a database connection
    conn = sqlite3.connect(database_path)
    cursor = conn.cursor()
    cursor.execute(sql_create_projects_table)
    conn.close()

    watershed_path_index_map_path = os.path.join(
        TARGET_WORKSPACE, 'watershed_rtree')
    print WATERSHED_PATH_LIST
    build_watershed_rtree(
        WATERSHED_PATH_LIST, watershed_path_index_map_path)

    grid_vector = gdal.OpenEx(DEGREE_GRID_PATH, gdal.OF_VECTOR)
    grid_layer = grid_vector.GetLayer()

    while True:
        grid_feature = grid_layer.GetNextFeature()
        if not grid_feature:
            break

        analyze_grid(
            grid_feature.GetFID(), watershed_path_index_map_path,
            database_path)


if __name__ == '__main__':
    main()
