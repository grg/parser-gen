#!/usr/bin/env python

"""A header node for use in a DAG"""

import Header
from HeaderLib import getHeaderLengths
from DAGNode import DAGNode
from DAGHeaderNode import HeaderNode

class SubHeaderNode(HeaderNode):
    """Sub-header node for use in a DAG"""
    
    def __init__(self, hdr, inst, startPos, endPos):
        self.startPos = startPos
        self.endPos = endPos
        super(self.__class__, self).__init__(hdr, inst)
        self.length = endPos - startPos + 1
        self._recalcStr()
#
#    def __str__(self):
#        return "HeaderNode: [HdrName='%s', Inst=%d, Length=%d]" % (self.hdr.name, self.inst, self.length)

    def getName(self):
        return '%s-sub-%d:%d' % (self.hdr.name, self.startPos, self.endPos)

    def getProcessLength(self):
        """Get the length of the header including any fields that need to be
        read beyond the end of header"""
        return self.length

    def getLengths(self):
        """Get all valid lengths"""
        return [self.length]

    def getDecisionBytes(self):
        """Get all decision byte positions"""
        decBytes = self.hdr.getDecisionBytes()[0]
        decBytes = [byte - self.startPos for byte in decBytes]
        while len(decBytes) > 0 and decBytes[-1] >= self.length:
            decBytes.pop()
        while len(decBytes) > 0 and decBytes[0] < 0:
            decBytes.pop(0)
        return decBytes

    def getExtractBytes(self):
        """Get all extract byte positions"""
        extBytes = self.hdr.getExtractBytes()
        extBytes = [byte - self.startPos for byte in extBytes]
        while len(extBytes) > 0 and extBytes[-1] >= self.length:
            extBytes.pop()
        while len(extBytes) > 0 and extBytes[0] < 0:
            extBytes.pop(0)
        return extBytes

    def getFields(self):
        """Get all fields within the header"""
        raise RuntimeError('This function needs vetting')
        return self.hdr.getLookupLengthFields()


# Basic test code
if __name__ == '__main__':
    hdr = Header.Header('TestHeader')
    hdr.addField('f0', 8)
    hdr.addField('f1', 8)
    hdr.addField('f2', 16)
    hdr.addField('f3', 8)
    hdr.addField('f4', 8)
    hdr.setCalcLength(['f1', '*', 2, '*', 'f3'])
    node = SubHeaderNode(hdr, 1, 1, 4)
    print "Node:", node
    print "Decision bytes:", node.getDecisionBytes()
