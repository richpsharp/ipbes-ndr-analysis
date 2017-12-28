# distutils: language=c++
"""Pitfilling module."""
import logging
import shutil
import tempfile
import time

import numpy
import pygeoprocessing
from osgeo import gdal

cimport numpy
cimport cython
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libc.math cimport exp
from libc.stdlib cimport malloc
from libc.stdlib cimport free
from libcpp.list cimport list
from libcpp.pair cimport pair
from libcpp.stack cimport stack

# This module expects rasters with a memory xy block size of 2**BLOCK_BITS
cdef int BLOCK_BITS = 8
NODATA = -1
cdef int _NODATA = NODATA

cdef bint isclose(double a, double b):
    return abs(a - b) <= (1e-5 + 1e-7 * abs(b))

cdef struct FlowPixel:
    int n_i
    int xi
    int yi
    double ret_eff

cdef extern from "LRUCache.h" nogil:
    cdef cppclass LRUCache[KEY_T, VAL_T]:
        LRUCache(int)
        void put(KEY_T&, VAL_T&, list[pair[KEY_T,VAL_T]]&)
        list[pair[KEY_T,VAL_T]].iterator begin()
        list[pair[KEY_T,VAL_T]].iterator end()
        bint exist(KEY_T &)
        VAL_T get(KEY_T &)

ctypedef pair[int, double*] BlockBufferPair

# TODO: make managed raster band aware
cdef class ManagedRaster:
    cdef LRUCache[int, double*]* lru_cache
    cdef int block_xsize
    cdef int block_ysize
    cdef int raster_x_size
    cdef int raster_y_size
    cdef int block_nx
    cdef int block_ny
    cdef int write_mode
    cdef int block_bits
    cdef char* raster_path

    def __cinit__(self, char* raster_path, n_blocks, write_mode):
        raster_info = pygeoprocessing.get_raster_info(raster_path)
        self.raster_x_size, self.raster_y_size = raster_info['raster_size']
        self.block_xsize, self.block_ysize = raster_info['block_size']
        self.block_bits = BLOCK_BITS
        self.block_xsize = 1<<self.block_bits
        self.block_ysize = 1<<self.block_bits
        """if (self.block_xsize != (1 << self.block_bits) or
                self.block_ysize != (1 << self.block_bits)):
            error_string = (
                "Expected block size that was %d bits wide, got %s" % (
                    self.block_bits, raster_info['block_size']))
            raise ValueError(error_string)
        """
        self.block_nx = (
            self.raster_x_size + (self.block_xsize) - 1) / self.block_xsize
        self.block_ny = (
            self.raster_y_size + (self.block_ysize) - 1) / self.block_ysize

        self.lru_cache = new LRUCache[int, double*](n_blocks)
        self.raster_path = raster_path
        self.write_mode = write_mode

    def __dealloc__(self):
        cdef int xi_copy, yi_copy
        cdef numpy.ndarray[double, ndim=2] block_array = numpy.empty(
            (self.block_ysize, self.block_xsize))
        cdef double *double_buffer
        cdef int block_xi
        cdef int block_yi
        # initially the win size is the same as the block size unless
        # we're at the edge of a raster
        cdef int win_xsize
        cdef int win_ysize

        # we need the offsets to subtract from global indexes for cached array
        cdef int xoff
        cdef int yoff

        cdef list[BlockBufferPair].iterator it = self.lru_cache.begin()
        cdef list[BlockBufferPair].iterator end = self.lru_cache.end()
        if not self.write_mode:
            while it != end:
                # write the changed value back if desired
                free(deref(it).second)
                inc(it)
            return

        raster = gdal.Open(self.raster_path, gdal.GA_Update)
        raster_band = raster.GetRasterBand(1)

        while it != end:
            # write the changed value back if desired
            double_buffer = deref(it).second

            if self.write_mode:
                block_index = deref(it).first
                block_xi = block_index % self.block_nx
                block_yi = block_index / self.block_nx

                # we need the offsets to subtract from global indexes for cached array
                xoff = block_xi * self.block_xsize
                yoff = block_yi * self.block_ysize

                win_xsize = self.block_xsize
                win_ysize = self.block_ysize

                # load a new block
                if xoff+win_xsize > self.raster_x_size:
                    win_xsize = win_xsize - (xoff+win_xsize - self.raster_x_size)
                if yoff+win_ysize > self.raster_y_size:
                    win_ysize = win_ysize - (yoff+win_ysize - self.raster_y_size)

                for xi_copy in xrange(win_xsize):
                    for yi_copy in xrange(win_ysize):
                        block_array[yi_copy, xi_copy] = (
                            double_buffer[yi_copy * self.block_xsize + xi_copy])
                raster_band.WriteArray(
                    block_array[0:win_ysize, 0:win_xsize],
                    xoff=xoff, yoff=yoff)
            free(double_buffer)
            inc(it)

        raster_band.FlushCache()
        raster_band = None
        raster = None

    cdef void set(self, int xi, int yi, double value):
        cdef int block_xi = xi >> self.block_bits
        cdef int block_yi = yi >> self.block_bits
        # this is the flat index for the block
        cdef int block_index = block_yi * self.block_nx + block_xi
        if not self.lru_cache.exist(block_index):
            self.load_block(block_index)
        cdef int xoff = block_xi << self.block_bits
        cdef int yoff = block_yi << self.block_bits
        self.lru_cache.get(
            block_index)[(yi-yoff)*self.block_xsize+xi-xoff] = value

    cdef double get(self, int xi, int yi):
        cdef int block_xi = xi >> self.block_bits
        cdef int block_yi = yi >> self.block_bits
        # this is the flat index for the block
        cdef int block_index = block_yi * self.block_nx + block_xi
        if not self.lru_cache.exist(block_index):
            self.load_block(block_index)
        cdef int xoff = block_xi << self.block_bits
        cdef int yoff = block_yi << self.block_bits
        return self.lru_cache.get(
            block_index)[(yi-yoff)*self.block_xsize+xi-xoff]

    cdef void get_block(self, int xc, int yc, double* block, int border) except *:
        """Load a block from the raster into `block`.

        Parameters:
            xc (int): x coordinate at center of block
            yy (int): y coordinate at center of block
            block (double*): pre-allocated block of size `(1+2*border)**2`
                after the call this array will contain either random pixels
                if the block lies outside of the raster, or the pixel values
                that overlap the block.
            border (int): number of pixels around `xc`, `yc` to load into
                `block`.

        Returns:
            None
        """
        cdef int block_xi
        cdef int block_yi
        cdef int block_index
        cdef int xi, yi
        cdef int xoff
        cdef int yoff

        for xi in xrange(xc-border, xc+border+1):
            if xi < 0 or xi >= self.raster_x_size:
                continue
            for yi in xrange(yc-border, yc+border+1):
                if yi < 0 or yi >= self.raster_y_size:
                    continue

                block_xi = xi >> self.block_bits
                block_yi = yi >> self.block_bits
                block_index = block_yi * self.block_nx + block_xi
                # this is the flat index for the block
                if not self.lru_cache.exist(block_index):
                    self.load_block(block_index)
                xoff = block_xi << self.block_bits
                yoff = block_yi << self.block_bits
                block[(yi-(yc-border))*(1+border*2)+xi-(xc-border)] = (
                    self.lru_cache.get(
                        block_index)[(yi-yoff)*self.block_xsize+xi-xoff])

    cdef void load_block(self, int block_index) except *:
        cdef int block_xi = block_index % self.block_nx
        cdef int block_yi = block_index / self.block_nx

        # we need the offsets to subtract from global indexes for cached array
        cdef int xoff = block_xi * self.block_xsize
        cdef int yoff = block_yi * self.block_ysize

        cdef int xi_copy, yi_copy
        cdef numpy.ndarray[double, ndim=2] block_array
        cdef double *double_buffer
        cdef list[BlockBufferPair] removed_value_list

        # determine the block aligned xoffset for read as array

        # initially the win size is the same as the block size unless
        # we're at the edge of a raster
        cdef int win_xsize = self.block_xsize
        cdef int win_ysize = self.block_ysize

        # load a new block
        if xoff+win_xsize > self.raster_x_size:
            win_xsize = win_xsize - (xoff+win_xsize - self.raster_x_size)
        if yoff+win_ysize > self.raster_y_size:
            win_ysize = win_ysize - (yoff+win_ysize - self.raster_y_size)

        raster = gdal.Open(self.raster_path)
        raster_band = raster.GetRasterBand(1)
        block_array = raster_band.ReadAsArray(
            xoff=xoff, yoff=yoff, win_xsize=win_xsize,
            win_ysize=win_ysize).astype(
            numpy.float64)
        raster_band = None
        raster = None
        double_buffer = <double*>malloc(
            sizeof(double) * self.block_xsize * win_ysize)
        for xi_copy in xrange(win_xsize):
            for yi_copy in xrange(win_ysize):
                double_buffer[yi_copy*self.block_xsize+xi_copy] = (
                    block_array[yi_copy, xi_copy])
        self.lru_cache.put(
            <int>block_index, <double*>double_buffer, removed_value_list)

        if self.write_mode:
            raster = gdal.Open(self.raster_path, gdal.GA_Update)
            raster_band = raster.GetRasterBand(1)

        block_array = numpy.empty(
            (self.block_ysize, self.block_xsize), dtype=numpy.double)
        while not removed_value_list.empty():
            # write the changed value back if desired
            double_buffer = removed_value_list.front().second

            if self.write_mode:
                block_index = removed_value_list.front().first
                block_xi = block_index % self.block_nx
                block_yi = block_index / self.block_nx

                # we need the offsets to subtract from global indexes for cached array
                xoff = block_xi * self.block_xsize
                yoff = block_yi * self.block_ysize

                win_xsize = self.block_xsize
                win_ysize = self.block_ysize

                # load a new block
                if xoff+win_xsize > self.raster_x_size:
                    win_xsize = win_xsize - (xoff+win_xsize - self.raster_x_size)
                if yoff+win_ysize > self.raster_y_size:
                    win_ysize = win_ysize - (yoff+win_ysize - self.raster_y_size)

                for xi_copy in xrange(win_xsize):
                    for yi_copy in xrange(win_ysize):
                        block_array[yi_copy, xi_copy] = (
                            double_buffer[yi_copy * self.block_xsize + xi_copy])
                raster_band.WriteArray(
                    block_array[0:win_ysize, 0:win_xsize],
                    xoff=xoff, yoff=yoff)
            free(double_buffer)
            removed_value_list.pop_front()

        if self.write_mode:
            raster_band.FlushCache()
            raster_band = None
            raster = None


def calculate_downstream_ret_eff(
        flow_dir_raster_path_band, flow_accum_raster_path_band,
        ret_eff_raster_path_band, double flow_threshold, double ret_len,
        target_downstream_retention_raster_path, temp_dir_path=None):
    """Calculates downstream retention by Equation 5 in UG.

    Parameters:
        flow_dir_raster_path_band (tuple): a path, band tuple to a flow
            direction raster calculated by `routing.fill_pits`.

            indicating the D8 flow direction raster with direction convention:
             321
             4 0
             567

        flow_accum_raster_path_band (tuple): a path/band tuple to a raster
            of same dimensions as the flow_dir_raster indicating values where
            flow flow length should terminate/start.
        ret_eff_raster_path_band (tuple): a raster path band tuple showing
            per pixel nutrient retention efficiency.
        flow_threshold (float): value in flow_accum_raster to indicate a
            sink for flow length if value is >= `flow_threshold`.
        target_downstream_retention_raster_path (string): path to output
            raster for effective downstream retention efficiency.
        temp_dir_path (string): if not None, a path to a directory where
            temporary files can be constructed. Otherwise uses system tempdir.

    Returns:
        None.
    """
    cdef numpy.ndarray[numpy.uint8_t, ndim=2] buffer_array
    cdef int raster_x_size, raster_y_size
    cdef int xi_n, yi_n, i
    cdef int xi, yi, win_ysize, win_xsize
    cdef int xoff, yoff
    cdef int flow_direction_nodata, n_dir, flow_dir
    cdef double s_i
    cdef double cell_size
    cdef double ret_eff_i, ret_eff_nodata
    cdef stack[FlowPixel] flow_stack
    cdef FlowPixel fp
    cdef ManagedRaster flow_accum_managed_raster
    cdef ManagedRaster flow_dir_managed_raster
    cdef ManagedRaster ret_eff_managed_raster

    logger = logging.getLogger(
        'pygeoprocessing.routing.downstream_flow_length')
    logger.addHandler(logging.NullHandler())  # silence logging by default

    # flow direction scheme is
    # 321
    # 4 0
    # 567
    # each flow direction is encoded as 1 << n, n in [0, 7]

    # use this to have offsets to visit neighbor pixels, pick 2 at a time to
    # add to a (xi, yi) tuple
    cdef int* OFFSET_ARRAY = [
        1, 0,  # 0
        1, -1,  # 1
        0, -1,  # 2
        -1, -1,  # 3
        -1, 0,  # 4
        -1, 1,  # 5
        0, 1,  # 6
        1, 1  # 7
        ]

    # this is used to set flow direction on a neighbor by indicating which
    # neighbor it flows to
    cdef int* REVERSE_FLOW_DIR = [
        4, # 0
        5, # 1
        6, # 2
        7, # 3
        0, # 4
        1, # 5
        2, # 6
        3 # 7
    ]

    # make an interesting temporary directory that has the time/date and
    # 'flow_accumulation' on it so we can figure out what's going on if we
    # ever run across it in a temp dir.
    temp_dir_path = tempfile.mkdtemp(
        dir=temp_dir_path, prefix='downstream_flow_length_',
        suffix=time.strftime('%Y-%m-%d_%H_%M_%S', time.gmtime()))

    pygeoprocessing.new_raster_from_base(
        flow_dir_raster_path_band[0],
        target_downstream_retention_raster_path, gdal.GDT_Float64,
        [_NODATA], fill_value_list=[_NODATA],
        gtiff_creation_options=(
            'TILED=YES', 'BIGTIFF=IF_SAFER', 'COMPRESS=LZW',
            'BLOCKXSIZE=%d' % (1<<BLOCK_BITS),
            'BLOCKYSIZE=%d' % (1<<BLOCK_BITS)))
    cell_size = pygeoprocessing.get_raster_info(
        target_downstream_retention_raster_path)['mean_pixel_size']

    # these are used to determine if a sample is within the raster
    flow_direction_raster_info = pygeoprocessing.get_raster_info(
        flow_dir_raster_path_band[0])
    flow_direction_nodata = flow_direction_raster_info['nodata'][
        flow_dir_raster_path_band[1]-1]
    raster_x_size, raster_y_size = flow_direction_raster_info['raster_size']

    flow_dir_managed_raster = ManagedRaster(
        flow_dir_raster_path_band[0], 2**9, 0)
    flow_accum_managed_raster = ManagedRaster(
        flow_accum_raster_path_band[0], 2**9, 0)
    ret_eff_managed_raster = ManagedRaster(
        ret_eff_raster_path_band[0], 2**9, 0)
    ret_eff_nodata = pygeoprocessing.get_raster_info(
        ret_eff_raster_path_band[0])['nodata'][
            ret_eff_raster_path_band[1]-1]
    downstream_retention_managed_raster = ManagedRaster(
        target_downstream_retention_raster_path, 2**9, 1)

    cdef double s_i_flat = exp(-5 * cell_size / ret_len)
    cdef double s_i_diag = exp(-5 * cell_size * 1.4142135 / ret_len)
    print s_i_flat
    print s_i_diag

    flow_direction_raster = gdal.Open(flow_dir_raster_path_band[0])
    flow_direction_band = flow_direction_raster.GetRasterBand(
        flow_dir_raster_path_band[1])

    logger.info('finding drains')
    start_drain_time = time.time()
    for offset_dict in pygeoprocessing.iterblocks(
            flow_dir_raster_path_band[0], offset_only=True,
            largest_block=0):
        # statically type these for later
        win_xsize = offset_dict['win_xsize']
        win_ysize = offset_dict['win_ysize']
        xoff = offset_dict['xoff']
        yoff = offset_dict['yoff']

        # make a buffer big enough to capture block and boundaries around it
        buffer_array = numpy.empty(
            (offset_dict['win_ysize']+2, offset_dict['win_xsize']+2),
            dtype=numpy.uint8)
        buffer_array[:] = flow_direction_nodata

        # default numpy array boundaries
        buffer_off = {
            'xa': 1,
            'xb': -1,
            'ya': 1,
            'yb': -1
        }
        # check if we can widen the border to include real data from the
        # raster
        for a_buffer_id, b_buffer_id, off_id, win_size_id, raster_size in [
                ('xa', 'xb', 'xoff', 'win_xsize', raster_x_size),
                ('ya', 'yb', 'yoff', 'win_ysize', raster_y_size)]:
            if offset_dict[off_id] > 0:
                # in thise case we have valid data to the left (or up)
                # grow the window and buffer slice in that direction
                buffer_off[a_buffer_id] = None
                offset_dict[off_id] -= 1
                offset_dict[win_size_id] += 1

            if offset_dict[off_id] + offset_dict[win_size_id] < raster_size:
                # here we have valid data to the right (or bottom)
                # grow the right buffer and add 1 to window
                buffer_off[b_buffer_id] = None
                offset_dict[win_size_id] += 1

        # read in the valid memory block
        buffer_array[
            buffer_off['ya']:buffer_off['yb'],
            buffer_off['xa']:buffer_off['xb']] = (
                flow_direction_band.ReadAsArray(
                    **offset_dict).astype(numpy.int8))

        # irrespective of how we sampled the DEM only look at the block in
        # the middle for valid
        for yi in xrange(1, win_ysize+1):
            for xi in xrange(1, win_xsize+1):
                flow_dir = (buffer_array[yi, xi])
                if isclose(flow_dir, flow_direction_nodata):
                    continue
                n_dir = buffer_array[
                    yi+OFFSET_ARRAY[2*flow_dir+1],
                    xi+OFFSET_ARRAY[2*flow_dir]]
                if (isclose(n_dir, flow_direction_nodata) or
                        flow_accum_managed_raster.get(
                            xi-1+xoff, yi-1+yoff) >= flow_threshold):
                    # it flows to nodata (or edge) so it's a seed
                    ret_eff_i = ret_eff_managed_raster.get(
                        xi-1+xoff, yi-1+yoff)
                    if isclose(ret_eff_i, ret_eff_nodata):
                        ret_eff_i = 0.0
                    eff_val = ret_eff_i * (1 - s_i_flat)
                    flow_stack.push(
                        FlowPixel(0, xi-1+xoff, yi-1+yoff, eff_val))
                    downstream_retention_managed_raster.set(
                        xi-1+xoff, yi-1+yoff, eff_val)

    logger.info("drains detected in %fs", time.time()-start_drain_time)
    while not flow_stack.empty():
        fp = flow_stack.top()
        flow_stack.pop()

        if (fp.xi == 0 or fp.xi == (raster_x_size-1) or
                fp.yi == 0 or fp.yi == (raster_y_size-1)):
            check_bounds_top = 1
        else:
            check_bounds_top = 0

        for i in xrange(fp.n_i, 8):
            # neighbor x,y indexes
            xi_n = fp.xi+OFFSET_ARRAY[2*i]
            yi_n = fp.yi+OFFSET_ARRAY[2*i+1]
            if check_bounds_top:
                if (xi_n < 0 or yi_n < 0 or
                        xi_n >= raster_x_size or yi_n >= raster_y_size):
                    continue
            if flow_dir_managed_raster.get(
                    xi_n, yi_n) == REVERSE_FLOW_DIR[i]:
                downstream_retention = (
                    downstream_retention_managed_raster.get(xi_n, yi_n))
                if isclose(downstream_retention, _NODATA):
                    # read upstream flow ret and see if we should update it
                    s_i = s_i_flat if i % 2 == 0 else s_i_diag
                    ret_eff_i = ret_eff_managed_raster.get(xi_n, yi_n)
                    if isclose(ret_eff_i, ret_eff_nodata):
                        ret_eff_i = 0.0
                    if ret_eff_i <= fp.ret_eff:
                        eff_val = fp.ret_eff
                    else:
                        eff_val = fp.ret_eff * s_i + ret_eff_i * (1 - s_i)
                    downstream_retention_managed_raster.set(
                        xi_n, yi_n, eff_val)
                    flow_stack.push(FlowPixel(0, xi_n, yi_n, eff_val))
    shutil.rmtree(temp_dir_path)
