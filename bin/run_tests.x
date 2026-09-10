!/bin/bash
set -e  # Exit on error
set -u  # Exit on unset variables
set -x  # Print commands as they are executed

bin/test_betweenness
bin/test_cellgeometry
bin/test_concom
bin/test_conts
bin/test_handle
bin/test_maxflow
bin/test_scc
bin/test_vtuio

