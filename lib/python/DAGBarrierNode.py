#!/usr/bin/env python

"""A header node for use in a DAG"""

import Header
from HeaderLib import getHeaderLengths
from DAGNode import DAGNode

class BarrierNode(DAGNode):
    """Barrier node for use in a DAG"""
    
    def __init__(self, hdr, inst, barrierLoc):
        self.barrierLoc = barrierLoc
        super(self.__class__, self).__init__(hdr, inst, 0)

    def getName(self):
        return '%s-bar%d' % (self.hdr.name, self.barrierLoc)

    def getDecisionBytes(self):
        """Get all decision byte positions"""
        return []

    def getExtractBytes(self):
        """Get all extract byte positions"""
        return []

    def getFields(self):
        """Get all fields within the header"""
        return []

# Basic test code
if __name__ == '__main__':
    hdr = Header.Header('TestHeader')
    hdr.addField('f1', 8)
    hdr.addField('f2', 16)
    hdr.addField('f3', 8)
    node = BarrierNode(hdr, 1, 2)
    print node
    print 'Total length:', node.getTotalLength()
