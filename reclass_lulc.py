"""Reclass the ESA lulc."""
import sys
import logging
import pandas

from osgeo import gdal
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


def build_lookup_from_csv(
        table_path, key_field, to_lower=True, warn_if_missing=True):
    """Read a CSV table into a dictionary indexed by `key_field`.

    Creates a dictionary from a CSV whose keys are unique entries in the CSV
    table under the column named by `key_field` and values are dictionaries
    indexed by the other columns in `table_path` including `key_field` whose
    values are the values on that row of the CSV table.

    Parameters:
        table_path (string): path to a CSV file containing at
            least the header key_field
        key_field: (string): a column in the CSV file at `table_path` that
            can uniquely identify each row in the table.
        to_lower (bool): if True, converts all unicode in the CSV,
            including headers and values to lowercase, otherwise uses raw
            string values.
        warn_if_missing (bool): If True, warnings are logged if there are
            empty headers or value rows.

    Returns:
        lookup_dict (dict): a dictionary of the form {
                key_field_0: {csv_header_0: value0, csv_header_1: value1...},
                key_field_1: {csv_header_0: valuea, csv_header_1: valueb...}
            }

        if `to_lower` all strings including key_fields and values are
        converted to lowercase unicode.
    """
    # Check if the file encoding is UTF-8 BOM first, related to issue
    # https://bitbucket.org/natcap/invest/issues/3832/invest-table-parsing-does-not-support-utf
    encoding = None
    with open(table_path) as file_obj:
        first_line = file_obj.readline()
        if first_line.startswith('\xef\xbb\xbf'):
            encoding = 'utf-8-sig'
    table = pandas.read_csv(
        table_path, sep=None, engine='python', encoding=encoding)
    header_row = list(table)
    key_field = str(key_field)
    if to_lower:
        key_field = key_field.lower()
        header_row = [
            x if not isinstance(x, str) else x.lower()
            for x in header_row]

    if key_field not in header_row:
        raise ValueError(
            '%s expected in %s for the CSV file at %s' % (
                key_field, header_row, table_path))
    if warn_if_missing and '' in header_row:
        LOGGER.warn(
            "There are empty strings in the header row at %s", table_path)

    key_index = header_row.index(key_field)
    lookup_dict = {}
    for index, row in table.iterrows():
        if to_lower:
            row = pandas.Series([
                x if not isinstance(x, str) else x.lower()
                for x in row])
        # check if every single element in the row is null
        if row.isnull().values.all():
            LOGGER.warn(
                "Encountered an entirely blank row on line %d", index+2)
            continue
        if row.isnull().values.any():
            row = row.fillna('')
        lookup_dict[row[key_index]] = dict(zip(header_row, row))
    return lookup_dict


def main():
    """Entry point."""
    lulc_raster_path = 'workspace_ipbes_ndr/churn/globio_landuse_scenarios/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7_md5_1254d25f937e6d9bdee5779d377c5aa4.tif'

    lulc_types_raster_path = 'esa_lulc_types.tif'
    table_path = 'esa_biophysical_table.csv'

    biophysical_table = build_lookup_from_csv(
        table_path, 'esa_id', to_lower=True, warn_if_missing=True)
    LOGGER.debug(biophysical_table)
    value_map = dict([
        (lucode, int(biophysical_table[lucode]['nathab'])
        for lucode in biophysical_table)])
    LOGGER.debug(value_map)
    pygeoprocessing.reclassify_raster(
        (lulc_raster_path, 1), value_map, lulc_types_raster_path, gdal.GDT_Byte,
        128, values_required=True)


if __name__ == '__main__':
    main()
