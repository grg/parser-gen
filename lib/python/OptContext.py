#!/usr/bin/env python

"""Optimization run result set"""

from abc import ABCMeta, abstractmethod

class OptContext(object):
    """Optimal chain graph node"""
    __metaclass__ = ABCMeta

    runs = 0

    def __init__(self):
        OptContext.runs += 1
        self.run = OptContext.runs
        self.reset()

    def reset(self):
        self.results = {}
        self.count = 0
        self.exploreResults = {}
        self.exploreResultCount = {}
        self.edgesPerCluster = {}
        self.shortestLenByCluster = {}
        self.fringeByCluster = {}
        self.fringeWithMinCoverLen = {}
        self.clustersFromNode = {}
        self.coveredClustersFromNode = {}
        self.dagOrder = {}
        self.dagOrderList = []
        self.dagParents = {}
        self.dag = None
        self.dagNodes = None
        self.optNodes = set()

    def __str__(self):
        return "%s: [id=%d]" % (self.__class__.__name__, self.run)

# Basic test code
if __name__ == '__main__':
    hdr = Header.Header('TestHeader')
    hdr.addField('test', 32)
    chainNode1 = DAGChainNode(HeaderNode(hdr, 1), 1, 3, 3)
    chainNode2 = DAGChainNode(HeaderNode(hdr, 2), 0, 2, 2)

    chain = DAGChain()
    chain.add(chainNode1)
    chain.add(chainNode2)

    optNode = OptNode(chain, 8, chain, None)
    print optNode

