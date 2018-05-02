#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Creates a HTTP file server accepting CORS and POST
99% Based on https://gist.github.com/jlesquembre/3922805 (all I did was refactoring it a little bit)
"""

from twisted.web.server import Site
from twisted.web.static import File
from twisted.internet import reactor
import argparse

class ResponseFile(File):
    def render(self, request):
        print request
        request.setHeader('Access-Control-Allow-Origin', '*')
        request.setHeader('Access-Control-Allow-Methods', 'POST, GET, OPTIONS')
        request.setHeader('Access-Control-Allow-Credentials', 'true')
        request.setHeader('Access-Control-Allow-Headers', 'origin, x-requested-with, accept, content-type, range')

        return File.render(self, request)

    def render_OPTIONS(self, request):
        return self.render_GET(request)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Serve files rooted at the output parameter. Accepts CORS and POST.')
    parser.add_argument("--output", help="Output folder served.", required=True)
    parser.add_argument("--port", type=int, help="The port to use", default=19974)
    args=parser.parse_args()

    root = ResponseFile(args.output)
    reactor.listenTCP(args.port, Site(root))
    print "Server is starting, now you can see the IGV plots"
    reactor.run()



