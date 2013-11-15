#!/usr/bin/env python

"""A node in a DAG"""

from abc import ABCMeta, abstractmethod
import Header

class DAGNode(object):
    __metaclass__ = ABCMeta
    """A node in a DAG"""
    
    def __init__(self, hdr, inst, length):
        self.hdr = hdr
        self.inst = inst
        self.length = length
        self.topoOrder = -1
        self.topoPos = -1
        self.nxt = set()
        self._recalcStr()

    def _recalcStr(self):
        self._str = "%s: [HdrName='%s', Inst=%d, Length=%d]" % (self.__class__.__name__, self.getName(), self.inst, self.length)
        self._shortStr = "%s-%d (l:%d)" % (self.getName(), self.inst, self.length)
        self._cmpName = '%s:%03d' % (self.getName(), self.inst)

    def __str__(self):
        return self._str

    def getName(self):
        return self.hdr.name

    def shortStr(self):
        #return "%s: [%s-%d (l:%d)]" % (self.__class__.__name__, self.getName(), self.inst, self.length)
        return self._shortStr

    def getLength(self):
        """Get the length of the node"""
        return self.length

    def getTotalLength(self):
        """Get the total length of the header"""
        return self.length

    def getProcessLength(self):
        """Get the length of the header including any fields that need to be
        read beyond the end of header"""
        return self.length

    def setTopo(self, order, pos):
        self.topoOrder = order
        self.topoPos = pos

    @abstractmethod
    def getDecisionBytes(self): pass

    @abstractmethod
    def getExtractBytes(self): pass

    @abstractmethod
    def getFields(self): pass

    def getCmpName(self):
        #return '%s:%03d' % (self.getName(), self.inst)
        return self._cmpName

    def __cmp__(self, other):
        if isinstance(other, DAGNode):
            return cmp(self.getCmpName(), other.getCmpName())
        else:
            return -1


# Basic test code
if __name__ == '__main__':
    class MyDAGNode(DAGNode):
        def getDecisionBytes(self):
            return []

    hdr = Header.Header('TestHeader')
    node = MyDAGNode(hdr, 1, 0)
    print node
