"""Convert Sqlite3 database into shapefile."""
import os
import datetime

import sqlite3
from osgeo import ogr
from osgeo import osr

def db_to_shapefile(database_path, target_dir):
    database_path = r"D:\rpsharp_documents\bitbucket_repos\ipbes-ndr-analysis\ndr_workspace\ipbes_ndr_results.db"
    target_dir = r'D:\Dropbox\rps_bck_shared_stuff\ipbes stuff\ipbes_ndr_results'
    target_shapefile = os.path.join(target_dir, 'results.shp')

    if os.path.exists(target_shapefile):
        os.remove(target_shapefile)

    # create a shapefile with these fields:
    # ws_id (text)
    # cur_n_exp
    # ssp1_n_exp
    # ssp3_n_exp
    # ssp5_n_exp

    target_sr = osr.SpatialReference()
    target_sr.ImportFromEPSG(4326)

    driver = ogr.GetDriverByName('ESRI Shapefile')
    result_vector = driver.CreateDataSource(target_shapefile)
    result_layer = result_vector.CreateLayer(
        os.path.splitext(os.path.basename(target_shapefile))[0],
        target_sr, ogr.wkbPolygon)
    ws_field = ogr.FieldDefn("ws_id", ogr.OFTString)
    ws_field.SetWidth(24)
    result_layer.CreateField(ws_field)

    scenario_list = ['cur', 'ssp1', 'ssp3', 'ssp5']
    for scenario in scenario_list:
        scenario_field = ogr.FieldDefn('%s_n_exp' % scenario, ogr.OFTReal)
        scenario_field.SetWidth(24)
        scenario_field.SetPrecision(11)
        result_layer.CreateField(scenario_field)

    conn = sqlite3.connect(database_path)
    if conn is not None:
        cursor = conn.cursor()
        cursor.execute(
            """SELECT ws_prefix_key, geometry_wkt FROM nutrient_export
               GROUP BY ws_prefix_key;""")
        result = cursor.fetchall()
        for ws_id, ws_geom in result:
            feature = ogr.Feature(result_layer.GetLayerDefn())
            feature.SetField('ws_id', ws_id)
            print ws_id
            for scenario in scenario_list:
                cursor.execute(
                    """SELECT total_export FROM nutrient_export
                    WHERE (ws_prefix_key = ? and scenario_key = ?)""", (
                        ws_id, scenario))
                feature.SetField('%s_n_exp' % scenario, cursor.fetchone()[0])
            feature.SetGeometry(ogr.CreateGeometryFromWkt(ws_geom))
            result_layer.CreateFeature(feature)

    timestring = datetime.datetime.now().strftime("%Y-%m-%d %H_%M_%S")
    timestamp_path = os.path.join(target_dir, 'last_update_%s.txt' % (
        timestring))
    with open(timestamp_path, 'w') as timestamp_file:
        timestamp_file.write(
            "Hi, I'm an automatically generated file.\n"
            "I last updated NDR results on %s.\n" % timestring +
            "There will be an 'all done.txt' file here when everything "
            "is done.\n")


if __name__ == '__main__':
    main()
