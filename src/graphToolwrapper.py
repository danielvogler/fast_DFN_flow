#! /usr/bin/python

# Copyright:
# 	Alex Hobe
# 	Daniel Vogler
# 	Martin P. Seybold

# TODO
# - explicitly wrap used functions.
#   - caution, this may make other functions unavailable

import graph_tool.all as gtAll # https://graph-tool.skewed.de/static/doc/index.html
import graph_tool.topology as topo
    

class graphToolwrapper(object):
    """graphToolwrapper: Used to minimizes changes required when graphtool changes.

    """
    def all():
        return gtAll 

    def topology():
        return topo
