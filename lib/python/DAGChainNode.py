#!/usr/bin/env python

"""A node in a DAGChain (contains DAGNode objects)"""

import Header
from DAGHeaderNode import HeaderNode

class DAGChainNode(object):
    """A node in a DAGChain"""
    
    def __init__(self, dagNode, startPos, consumed, read):
        """
        Create a new node

        Parameters:
          dagNode - dag node object
          startPos - starting position within the dagNode
          consumed - number of bytes (from the startPos) consumed
          read - number of bytes (from the startPos) read
        """
        dagLen = dagNode.getLength()
        maxLen = dagNode.getProcessLength()

        if dagLen != 0 and startPos >= maxLen:
            raise ValueError('Start position %d is beyond node end (%s)' %
                    (startPos, str(dagNode)))

        if startPos + consumed - 1 > dagLen:
            consumed = dagLen - startPos
        if consumed < 0:
            consumed = 0

        if startPos + read - 1 > maxLen:
            read = maxLen - startPos
        if read < 0:
            read = 0

        if read < consumed:
            raise ValueError('Read bytes %d cannot be less than consumed bytes %d (%s)' %
                    (startPos, str(dagNode)))


        self._dagNode = dagNode
        self._startPos = startPos
        self._consumed = consumed
        self._read = read
        self._str = None
        self.__recalcHash()

    @property
    def dagNode(self):
        """Property for storing the original dagNode"""
        return self._dagNode

    @dagNode.setter
    def dagNode(self, value):
        self._dagNode = value
        self.__recalcHash()

    @property
    def startPos(self):
        """Property for storing the original startPos"""
        return self._startPos

    @startPos.setter
    def startPos(self, value):
        self._startPos = value
        self.__recalcHash()

    @property
    def consumed(self):
        """Property for storing the original consumed"""
        return self._consumed

    @consumed.setter
    def consumed(self, value):
        self._consumed = value
        self.__recalcHash()

    @property
    def read(self):
        """Property for storing the original read"""
        return self._read

    @read.setter
    def read(self, value):
        self._read = value
        self.__recalcHash()

    def __str__(self):
        #return "%s: [DagNode='%s', StartPos=%d, Span=%d]" % \
        #        (self.__class__.__name__, self.dagNode.shortStr(), self.startPos, self.span)
        if self._str:
            return self._str

        if self.consumed > 0 and self.read > 0:
            self._str = "[%s %d:%d/%d]" % \
                    (self.dagNode.shortStr(), self.startPos,
                            self.startPos + self.consumed - 1,
                            self.startPos + self.read - 1)
        elif self.read > 0:
            self._str = "[%s %d:--/%d]" % \
                    (self.dagNode.shortStr(), self.startPos,
                            self.startPos + self.read - 1)
        else:
            self._str = "[%s %d:--/--]" % \
                    (self.dagNode.shortStr(), self.startPos)
        return self._str

    def __hash__(self):
        return self._hash

    def __recalcHash(self):
        self._str = None
        self._hash = hash(self.__str__())

    def __eq__(self, other):
        #if isinstance(other, DAGChainNode):
        return self.__hash__() == other.__hash__()

    def __cmp__(self, other):
        #if isinstance(other, DAGChainNode):
        if type(other) == DAGChainNode:
            if self.__hash__() == other.__hash__():
                return 0

            c = cmp(self.dagNode, other.dagNode)
            if c != 0:
                return c

            c = cmp(self.startPos, other.startPos)
            if c != 0:
                return c

            c = cmp(self.consumed, other.consumed)
            if c != 0:
                return c

            c = cmp(self.read, other.read)
            if c != 0:
                return c

            return 0
        else:
            return -1

    def unconsumed(self):
        return self.dagNode.getLength() - self.startPos - self.consumed

    def unread(self):
        procLength = self.dagNode.getProcessLength()

        return procLength - self.startPos - self.read

    def unprocLookupBytes(self):
        """
        Lookup bytes that haven't been read (indexed from start + read)
        """
        unproc = []
        lookupBytes = self.dagNode.getDecisionBytes()
        for lb in lookupBytes:
            if lb >= self.startPos + self.read:
                unproc.append(lb - self.startPos - self.read)

        return unproc

    def unprocExtractBytes(self):
        """
        Extract bytes that haven't been processed (indexed from start + read)
        """
        unproc = []
        extractBytes = self.dagNode.getExtractBytes()
        for eb in extractBytes:
            if eb >= self.startPos + self.read:
                unproc.append(eb - self.startPos - self.read)

        return unproc


# Basic test code
if __name__ == '__main__':
    hdr = Header.Header('TestHeader')
    hdr.addField('test', 32)
    dagNode1 = HeaderNode(hdr, 1)
    dagNode2 = HeaderNode(hdr, 2)
    chainNode1 = DAGChainNode(dagNode1, 1, 4, 4)
    chainNode2 = DAGChainNode(dagNode1, 1, 4, 4)
    chainNode3 = DAGChainNode(dagNode2, 1, 4, 4)
    print "Chain node 1:", chainNode1
    print "Chain node 2:", chainNode2
    print "Chain node 3:", chainNode3
    print "Node 1 == Node 2?", chainNode1 == chainNode2
    print "Node 1 == Node 3?", chainNode1 == chainNode3
