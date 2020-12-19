import os
import sys
import logging

import ipbes_ndr_analysis
LOGGER = logging.getLogger(__name__)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        LOGGER.error(
            "usage: python %s workspace_dir", sys.argv[0])
        sys.exit(-1)
    raw_workspace_dir = sys.argv[1]
    if os.path.isfile(raw_workspace_dir):
        LOGGER.error(
            '%s is supposed to be the workspace directory but points to an '
            'existing file' % raw_workspace_dir)
        sys.exit(-1)
    ipbes_ndr_analysis.main(raw_workspace_dir)
