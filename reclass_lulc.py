"""Reclass the ESA lulc."""
import sys
import logging

from osgeo import gdal
import natcap.invest.utils
import pygeoprocessing

# set a 1GB limit for the cache
gdal.SetCacheMax(2**30)

logging.basicConfig(
    level=logging.DEBUG,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(funcName)s:%(lineno)d] %(message)s'),
    stream=sys.stdout)
LOGGER = logging.getLogger(__name__)


def main():
    """Entry point."""
    lulc_raster_path = 'workspace_ipbes_ndr/churn/globio_landuse_scenarios/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7_md5_1254d25f937e6d9bdee5779d377c5aa4.tif'

    lulc_types_raster_path = 'esa_lulc_types.tif'
    table_path = 'esa_biophysical_table.csv'

    biophysical_table = natcap.invest.utils.build_lookup_from_csv(
        table_path, 'esa_id', to_lower=True, warn_if_missing=True)
    LOGGER.debug(biophysical_table)
    value_map = dict([
        (lucode, int(biophysical_table[lucode]['nathab']))]
        for lucode in biophysical_table)
    LOGGER.debug(value_map)
    pygeoprocessing.reclassify_raster(
        (lulc_raster_path, 1), value_map, lulc_types_raster_path, gdal.GDT_Byte,
        128, values_required=True)


if __name__ == '__main__':
    main()
