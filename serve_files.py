#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Creates a HTTP file server using npm http-server
"""

import argparse
from ExternalTools import executeCommandLine

parser = argparse.ArgumentParser(description='Serve files rooted at the output parameter. Uses npm http-server')
parser.add_argument("output", help="Output folder served.")
parser.add_argument("--port", type=int, help="The port to use", default=19974)
args=parser.parse_args()

try:
    executeCommandLine("http-server --cors -p %d %s"%(args.port, args.output))
except KeyboardInterrupt:
    pass