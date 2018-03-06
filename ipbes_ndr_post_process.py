"""
The tasks in 3/2/2018 of this design doc https://docs.google.com/document/d/1Pdk-1tcC2TZKdOgb4NZTDf4Qpdu_6dGIsJ2NHiCkY3Y/edit
"""
import os
import logging

from osgeo import gdal
import numpy
import pygeoprocessing
import taskgraph

logging.basicConfig(
    format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
    level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')
LOGGER = logging.getLogger('ipbes-ndr-post-process')
LOGGER.setLevel(logging.DEBUG)

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


WORKSPACE_DIR = os.path.join(
    BASE_DROPBOX_DIR, 'ipbes stuff', 'ipbes_ndr_results',
    'ndr_final_set_of_files')

NODATA = -1.0


class Divide(object):
    """num/denom"""
    def __init__(self, num_path, denom_path):
        self.num_nodata = pygeoprocessing.get_raster_info(
            num_path)['nodata'][0]
        self.denom_nodata = pygeoprocessing.get_raster_info(
            denom_path)['nodata'][0]

    def __call__(self, num, denom):
        result = numpy.empty_like(num)
        result[:] = NODATA
        valid_mask = (
            (num != self.num_nodata) & (denom != self.denom_nodata) &
            (denom != 0))
        result[valid_mask] = num[valid_mask] / denom[valid_mask]
        result[denom == 0.0] = 0.0
        return result


class Mult(object):
    """vala*valb"""
    def __init__(self, vala_path, valb_path):
        self.vala_nodata = pygeoprocessing.get_raster_info(
            vala_path)['nodata'][0]
        self.valb_nodata = pygeoprocessing.get_raster_info(
            valb_path)['nodata'][0]

    def __call__(self, vala, valb):
        result = numpy.empty_like(vala)
        result[:] = NODATA
        valid_mask = (
            (vala != self.vala_nodata) & (valb != self.valb_nodata))
        result[valid_mask] = vala[valid_mask] * valb[valid_mask]
        return result


class PropDiff(object):
    """(fut-cur)/cur"""
    def __init__(self, cur_path, fut_path):
        self.cur_nodata = pygeoprocessing.get_raster_info(
            cur_path)['nodata'][0]
        self.fut_nodata = pygeoprocessing.get_raster_info(
            fut_path)['nodata'][0]

    def __call__(self, cur, fut):
        result = numpy.empty_like(cur)
        result[:] = NODATA
        valid_mask = (
            (cur != self.cur_nodata) & (fut != self.fut_nodata) &
            (cur != 0))
        result[valid_mask] = (
            fut[valid_mask] - cur[valid_mask]) / cur[valid_mask]
        result[cur == 0.0] = 0.0
        return result


class Log(object):
    """log(val)"""
    def __init__(self, path):
        self.nodata = pygeoprocessing.get_raster_info(path)['nodata'][0]

    def __call__(self, val):
        result = numpy.empty_like(val)
        result[:] = NODATA
        valid_mask = (val != self.nodata) &  (val > 0)
        result[valid_mask] = numpy.log(val[valid_mask])
        return result


def main():
    """
    SvNEx_[cur|ssp[1|3|5]] = [cur|ssp[1|3|5]]_service / [cur|ssp[1|3|5]]_n_export
    cSvNEx_[ssp[1|3|5] = (SvNEx_ssp[1|3|5]] -SvNEx_cur) /SvNEx_cur
    pSvNEx = SvNEx_[cur|ssp[1|3|5]] * [cur|ssp[1|3|5]]_gpwpop_rural_degree
    cpSvNEx_[cur|ssp[1|3|5]] = pSvNEx_ssp[1|3|5]] - pSvNEx_cur) / pSvNEx_cur
    pNEx_[cur|ssp[1|3|5]]  = [cur|ssp[1|3|5]]_n_export * [cur|ssp[1|3|5]]_gpwpop_rural_degree
    cpNEx_ssp[1|3|5] = pNEx_ssp[1|3|5]] - pNEx_cur) / pNEx_cur
    """
    try:
        os.makedirs(WORKSPACE_DIR)
    except OSError:
        pass

    task_graph = taskgraph.TaskGraph(os.path.join(
        WORKSPACE_DIR, 'task_graph_dir'), -1)

    gpw_rescale_path_task_map = {}
    for scenario_id in ['cur', 'ssp1', 'ssp3', 'ssp5']:
        gpwpop_path = os.path.join(
            BASE_DROPBOX_DIR, 'ipbes stuff', 'population stuff',
            'rural_scenario_population',
            '%s_gpwpop_rural_degree.tif' % scenario_id)
        gpwpop_rescale_path = os.path.join(
            WORKSPACE_DIR, 'gpwpop_%s_rescale.tif' % scenario_id)

        LOGGER.debug(os.path.exists(gpwpop_path))

        task = task_graph.add_task(
            func=pygeoprocessing.warp_raster,
            args=(
                gpwpop_path, [1.0, -1.0], gpwpop_rescale_path, 'nearest'),
            kwargs={'target_bb': [-180.0, -91.0, 181.0, 90.0]},
            target_path_list=[gpwpop_rescale_path])
        gpw_rescale_path_task_map[scenario_id] = (gpwpop_rescale_path, task)

    # SvNEx_[cur|ssp[1|3|5]] = [cur|ssp[1|3|5]]_service / [cur|ssp[1|3|5]]_n_export_degree
    svnex_path_tasks = {}
    for scenario_id in ['cur', 'ssp1', 'ssp3', 'ssp5']:
        service_path = os.path.join(
            BASE_DROPBOX_DIR, 'ipbes stuff', 'ipbes_ndr_results',
            'n_export_and_n_load_tifs_at_degree',
            '%s_service.tif' % scenario_id)
        export_path = os.path.join(
            BASE_DROPBOX_DIR, 'ipbes stuff', 'ipbes_ndr_results',
            'n_export_and_n_load_tifs_at_degree',
            '%s_n_export_degree.tif' % scenario_id)
        svnex_path = os.path.join(WORKSPACE_DIR, 'SvNEx_%s.tif' % scenario_id)

        task = task_graph.add_task(
            func=pygeoprocessing.raster_calculator,
            args=(
                [(service_path, 1), (export_path, 1)], Divide(
                    service_path, export_path), svnex_path, gdal.GDT_Float32,
                NODATA),
            target_path_list=[svnex_path],
            task_name='svnex_%s' % scenario_id)

        svnex_path_tasks[scenario_id] = (svnex_path, task)

    # [cur|ssp[1|3|5]]_service / [cur|ssp[1|3|5]]_n_export
    for scenario_id in ['ssp1', 'ssp3', 'ssp5']:
        svnex_cur_path = svnex_path_tasks['cur'][0]
        svnex_fut_path = svnex_path_tasks[scenario_id][0]
        cSvNEx_path = os.path.join(
            WORKSPACE_DIR, 'cSvNEx_%s.tif' % scenario_id)

        task_graph.add_task(
            func=pygeoprocessing.raster_calculator,
            args=(
                [(svnex_cur_path, 1), (svnex_fut_path, 1)], PropDiff(
                    svnex_cur_path, svnex_fut_path), cSvNEx_path,
                gdal.GDT_Float32, NODATA),
            target_path_list=[cSvNEx_path],
            dependent_task_list=[
                svnex_path_tasks['cur'][1], svnex_path_tasks[scenario_id][1]],
            task_name='cSvNEx_%s' % scenario_id)

    #pSvNEx = SvNEx_[cur|ssp[1|3|5]] * [cur|ssp[1|3|5]]_gpwpop_rural_degree
    pSvNEx_path_tasks = {}
    for scenario_id in ['cur', 'ssp1', 'ssp3', 'ssp5']:
        svnex_path = svnex_path_tasks[scenario_id][0]
        gpwpop_path = gpw_rescale_path_task_map[scenario_id][0]
        pSvNEx_path = os.path.join(
            WORKSPACE_DIR, 'pSvNEx_%s.tif' % scenario_id)

        task = task_graph.add_task(
            func=pygeoprocessing.raster_calculator,
            args=(
                [(svnex_path, 1), (gpwpop_path, 1)], Mult(
                    svnex_path, gpwpop_path), pSvNEx_path,
                gdal.GDT_Float32, NODATA),
            target_path_list=[pSvNEx_path],
            dependent_task_list=[
                svnex_path_tasks[scenario_id][1],
                gpw_rescale_path_task_map[scenario_id][1]],
            task_name='cSvNEx_%s' % scenario_id)

        pSvNEx_path_tasks[scenario_id] = (pSvNEx_path, task)

    #cpSvNEx_[cur|ssp[1|3|5]] = pSvNEx_ssp[1|3|5]] - pSvNEx_cur) / pSvNEx_cur
    for scenario_id in ['ssp1', 'ssp3', 'ssp5']:
        pSvNEx_cur_path = pSvNEx_path_tasks['cur'][0]
        pSvNEx_fut_path = pSvNEx_path_tasks[scenario_id][0]
        cpSvNEx_path = os.path.join(WORKSPACE_DIR, 'cpSvNEx_%s.tif' % scenario_id)

        task_graph.add_task(
            func=pygeoprocessing.raster_calculator,
            args=(
                [(pSvNEx_cur_path, 1), (pSvNEx_fut_path, 1)], PropDiff(
                    pSvNEx_cur_path, pSvNEx_fut_path), cpSvNEx_path,
                gdal.GDT_Float32, NODATA),
            target_path_list=[cpSvNEx_path],
            dependent_task_list=[
                pSvNEx_path_tasks['cur'][1], pSvNEx_path_tasks[scenario_id][1]],
            task_name='cpSvNEx_%s' % scenario_id)

    #pNEx_[cur|ssp[1|3|5]]  = [cur|ssp[1|3|5]]_n_export * [cur|ssp[1|3|5]]_gpwpop_rural_degree
    pNEx_path_tasks = {}
    for scenario_id in ['cur', 'ssp1', 'ssp3', 'ssp5']:
        n_export_path = os.path.join(
            BASE_DROPBOX_DIR, 'ipbes stuff', 'ipbes_ndr_results',
            'n_export_and_n_load_tifs_at_degree',
            '%s_n_export_degree.tif' % scenario_id)
        gpwpop_path = gpw_rescale_path_task_map[scenario_id][0]
        pNEx_path = os.path.join(WORKSPACE_DIR, 'pNEx_%s.tif' % scenario_id)

        task = task_graph.add_task(
            func=pygeoprocessing.raster_calculator,
            args=(
                [(n_export_path, 1), (gpwpop_path, 1)], Mult(
                    n_export_path, gpwpop_path), pNEx_path,
                gdal.GDT_Float32, NODATA),
            target_path_list=[pNEx_path],
            dependent_task_list=[gpw_rescale_path_task_map[scenario_id][1]],
            task_name='pNEx_%s' % scenario_id)

        pNEx_path_tasks[scenario_id] = (pNEx_path, task)

    #cpNEx_[cur|ssp[1|3|5]] = pNEx_ssp[1|3|5]] - pNEx_cur) / pNEx_cur
    for scenario_id in ['ssp1', 'ssp3', 'ssp5']:
        pNEx_cur_path = pNEx_path_tasks['cur'][0]
        pNEx_fut_path = pNEx_path_tasks[scenario_id][0]
        cpNEx_path = os.path.join(WORKSPACE_DIR, 'cpNEx_%s.tif' % scenario_id)

        task_graph.add_task(
            func=pygeoprocessing.raster_calculator,
            args=(
                [(pNEx_cur_path, 1), (pNEx_fut_path, 1)], PropDiff(
                    pNEx_cur_path, pNEx_fut_path), cpNEx_path,
                gdal.GDT_Float32, NODATA),
            target_path_list=[cpNEx_path],
            dependent_task_list=[
                pNEx_path_tasks['cur'][1], pNEx_path_tasks[scenario_id][1]],
            task_name='cpSvNEx_%s' % scenario_id)

    # logrurpop_[cur|ssp[1|3|5]] = log(cur|ssp[1|3|5]_gpwpop_rural_degree)
    logrurpop_path_task_map = {}
    for scenario_id in ['cur', 'ssp1', 'ssp3', 'ssp5']:
        logrurpop_path = os.path.join(
            WORKSPACE_DIR, 'logrurpop_%s' % scenario_id)

        gwppop_rural_path = gpw_rescale_path_task_map[scenario_id][0]
        task = task_graph.add_task(
            func=pygeoprocessing.raster_calculator,
            args=(
                [(gwppop_rural_path, 1)], Log(gwppop_rural_path),
                logrurpop_path, gdal.GDT_Float32, NODATA),
            target_path_list=[logrurpop_path],
            dependent_task_list=[gpw_rescale_path_task_map[scenario_id][1]],
            task_name='logrurpop_%s' % scenario_id)
        logrurpop_path_task_map[scenario_id] = (
            logrurpop_path, task)

    # clogrurpop = (logrurpop_ssp[1|3|5] - logrurpop_cur)/logrurpop_cur
    for scenario_id in ['ssp1', 'ssp3', 'ssp5']:
        pass
    # cNEx_ssp[1|3|5] = NEx_ssp[1|3|5]] - NEx_cur) / NEx_cur
    # pcNEx__[ssp[1|3|5] = logrurpop_[cur|ssp[1|3|5]] * cNEx__[ssp[1|3|5]
    # pcSvNEx_[ssp[1|3|5]= logrurpop_[cur|ssp[1|3|5]] * cSvNEx_[ssp[1|3|5]

    task_graph.close()
    task_graph.join()

if __name__ == '__main__':
    main()
