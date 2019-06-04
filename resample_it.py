"""Resample inputs."""
import os
import sys
import subprocess
import glob

import pygeoprocessing


def main():
    factor = 10
    for pattern in sys.argv[1::]:
        for path in glob.glob(pattern):
            target_path = '%s_%d.tif' % (
                os.path.splitext(os.path.basename(path))[0], factor)
            raster_info = pygeoprocessing.get_raster_info(path)
            print(target_path)
            subprocess.check_call([
                'gdalwarp', path, target_path, '-r', 'average', '-tr',
                str(raster_info['pixel_size'][0]*factor),
                str(raster_info['pixel_size'][1]*factor)])


if __name__ == '__main__':
    main()
