#!/usr/bin/env python

"""A header node for use in a DAG"""

import Header
from HeaderLib import getHeaderLengths
from DAGNode import DAGNode

class PadNode(DAGNode):
    """Pad node for use in a DAG"""
    
    def getName(self):
        return self.hdr.name + '-pad-l%d' % self.length

    def getTotalLength(self):
        """Get the total length of the header"""
        return self.length + self.hdr.length()[0] / 8

    def getDecisionBytes(self):
        """Get all decision byte positions"""
        return []

    def getExtractBytes(self):
        """Get all extract byte positions"""
        return []

    def getFields(self):
        """Get all fields within the header"""
        return []

    def _recalcStr(self):
        super(PadNode, self)._recalcStr()
        self._cmpName = '%s:%03d:%03d' % (self.getName(), self.inst, self.getLength())

    #def getCmpName(self):
    #    return '%s:%03d:%03d' % (self.getName(), self.inst, self.getLength())


# Basic test code
if __name__ == '__main__':
    hdr = Header.Header('TestHeader')
    hdr.addField('test', 8)
    node = PadNode(hdr, 1, 4)
    print node
    print 'Total length:', node.getTotalLength()
