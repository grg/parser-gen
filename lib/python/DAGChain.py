#!/usr/bin/env python

"""A chain of nodes from a DAG"""

import Header
from DAGNode import DAGNode
from DAGHeaderNode import HeaderNode
from DAGSubHeaderNode import SubHeaderNode
from DAGBarrierNode import BarrierNode
from DAGPadNode import PadNode
from DAGChainNode import DAGChainNode
import copy

class DAGChain(object):
    """A chain of DAG nodes"""

    singleEntryForInternal = True
    
    def __init__(self):
        self._chain = []
        self._done = False
        self._patterns = 1
        self._hash = None
        self._successors = None
        self.__recalcHash()

    @property
    def chain(self):
        """Property containing the list of chain nodes"""
        return self._chain

    @chain.setter
    def chain(self, value):
        self._chain = value
        self.__recalcHash()

    @property
    def done(self):
        """Property indicating if the chain is done"""
        return self._done

    @done.setter
    def done(self, value):
        self._done = value
        self.__recalcHash()

    @property
    def patterns(self):
        """Property indicating the number of patterns"""
        return self._patterns

    @patterns.setter
    def patterns(self, value):
        self._patterns = value
        self.__recalcHash()

    def __str__(self):
        chainStr = ""
        first = True
        for chainNode in self.chain:
            if not first:
                chainStr += " -> "
            chainStr += str(chainNode)
            first = False

        return "[%s, done=%s, patterns=%d]" % \
                (chainStr, self.done, self.patterns)
        #return "%s: [[%s], done=%s, patterns=%d]" % \
        #        (self.__class__.__name__, chainStr, self.done, self.patterns)

    def shortStr(self):
        chainStr = ""
        first = True
        for node in self.chain:
            if not first:
                chainStr += '  '
            chainStr += node.dagNode.shortStr()
            chainStr += '-p%02d-c%02d-r%02d' % (node.startPos, node.consumed, node.read)
            first = False
        return "%s: [%s]" % (self.__class__.__name__, chainStr)

    def __key(self):
        chainStr = ""
        first = True
        for node in self.chain:
            if not first:
                chainStr += '-'
            chainStr += node.dagNode.shortStr()
            chainStr += '-p%02d-c%02d-r%02d' % (node.startPos, node.consumed, node.read)
            first = False

        return (self.__class__, chainStr, self.done, self.patterns)

    def key(self):
        return self.__key()

    def __hash__(self):
        return self._hash

    def __recalcHash(self):
        self._successors = None
        self._hash = hash(self.__key())

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

    def __cmp__(self, other):
        if isinstance(other, DAGChain):
            if self.__hash__() == other.__hash__():
                return 0

            for i in xrange(min(len(self.chain), len(other.chain))):
                c = cmp(self.chain[i], other.chain[i])
                if c != 0:
                    return c

            if len(self.chain) < len(other.chain):
                return -1
            elif len(self.chain) > len(other.chain):
                return 1

            c = cmp(self.done, other.done)
            if c != 0:
                return c

            c = cmp(self.patterns, other.patterns)
            if c != 0:
                return c

            return 0
        else:
            return -1

    def chainPop(self):
        """Pop a node from the end of the chain"""
        self.chain.pop()
        self.__recalcHash()

    def add(self, chainNode, startPos=None, consumed=None, read=None):
        """Add a chain node to this chain"""
        if chainNode is None:
            self.done = True
        else:
            if self.done:
                raise RuntimeError('Attempt to add a node to a done DAGChain')

            if isinstance(chainNode, DAGChainNode):
                self.chain.append(chainNode)
            elif isinstance(chainNode, DAGNode):
                if consumed is None:
                    consumed = chainNode.getLength() - startPos
                if read is None:
                    if consumed == chainNode.getLength() - startPos:
                        read = chainNode.getProcessLength() - startPos
                    else:
                        read = consumed
                self.chain.append(DAGChainNode(chainNode, startPos, consumed, read))
            else:
                raise TypeError('Incorrect type %s' % type(chainNodes))
        self.calcPatterns()

    def trim(self):
        """Trim any nodes with consumed/read 0 from the end of the chain"""
        while len(self.chain) > 0 and self.chain[-1].consumed == 0 and self.chain[-1].read == 0:
            self.chain.pop()
        self.__recalcHash()

    def extend(self, chainNodes):
        """Add multiple chain nodes to this chain"""
        if self.done:
            raise RuntimeError('Attempt to extend a chain whose done flag is set')
        elif type(chainNodes) == list:
            self.chain.extend(chainNodes)
        elif isinstance(chainNodes, DAGChain):
            if len(self.chain) > 0 and len(chainNodes.chain) > 0:
                myLastNode = self.chain[-1]
                theirFirstNode = chainNodes.chain[0]
                if myLastNode.dagNode == theirFirstNode.dagNode:
                    if myLastNode.startPos + myLastNode.consumed != theirFirstNode.startPos:
                        raise RuntimeError('Attempt to extend a chain with another chain that doesn\'t align to end')
                    else:
                        newLastNode = copy.copy(myLastNode)
                        newLastNode.consumed += theirFirstNode.consumed
                        newLastNode.read = theirFirstNode.read + theirFirstNode.startPos - myLastNode.startPos
                        self.chain[-1] = newLastNode
                        self.chain.extend(chainNodes.chain[1:])
                else:
                    self.chain.extend(chainNodes.chain)
            else:
                self.chain.extend(chainNodes.chain)
            self.done = chainNodes.done
        else:
            raise TypeError('Incorrect type %s' % type(chainNodes))

    def totalConsumed(self):
        """Return the total consumed of all nodes"""
        consumed = 0
        for node in self.chain:
            consumed += node.consumed
        return consumed

    def totalRead(self):
        """Return the total read of all nodes"""
        read = self.totalConsumed()
        try:
        #if len(self.chain) > 0:
            read += self.chain[-1].read - self.chain[-1].consumed
        except:
            pass
        return read

    def nextDAGNodes(self):
        """What DAG nodes can follow this node?"""
        if len(self.chain) > 0 and not self.done:
            return self.chain[-1].dagNode.nxt
        else:
            return []

    def copy(self):
        """Create a copy of this node"""
        chainCopy = DAGChain()
        for cn in self.chain:
            cnCopy = copy.copy(cn)
            chainCopy.add(cnCopy)
        #chainCopy.chain = copy.copy(self.chain)
        chainCopy.done = self.done
        chainCopy.patterns = self.patterns
        return chainCopy

    def copyToDepth(self, depth):
        """Create a copy of this node up to a specified depth"""
        chainCopy = DAGChain()
        pos = 0
        index = 0
        subset = False
        while pos < depth:
            cn = copy.copy(self.chain[index])
            if pos + cn.consumed > depth:
                cn.consumed = depth - pos
                cn.read = depth - pos
                subset = True
            chainCopy.add(cn)
            pos += cn.consumed
            index += 1
            if index == len(self.chain):
                pos = depth
        if not subset:
            chainCopy.done = self.done
            chainCopy.patterns = self.patterns
        else:
            chainCopy.calcPatterns()
        print "ChainCopy: %s   (Depth: %d)" % (chainCopy, depth)

        return chainCopy

    def isDone(self):
        return self.done

    def lastNode(self, ignoreDone=False):
        if self.done and not ignoreDone:
            return None
        else:
            return self.chain[-1]

    def nodeCount(self):
        return len(self.chain)

    def getLookupBytes(self):
        """Get all lookup bytes for this chain"""
        # Get a list of lookup bytes across the chain, indexed from the first startPos
        offset = 0
        lookupBytes = set()
        for node in self.chain:
            nodeLookupBytes = node.dagNode.getDecisionBytes()
            for lb in nodeLookupBytes:
                if lb >= node.startPos and lb < node.startPos + node.read:
                    lookupBytes.add(lb - node.startPos + offset)
            offset += node.consumed

        lookupBytes = sorted(lookupBytes)
        #print "LookupBytes:", lookupBytes

        return lookupBytes

    def getExtractBytes(self):
        """Get all extract bytes for this chain"""
        # Get a list of extract bytes across the chain, indexed from the first startPos
        offset = 0
        extractBytes = set()
        for node in self.chain:
            nodeExtractBytes = node.dagNode.getExtractBytes()
            for lb in nodeExtractBytes:
                if lb >= node.startPos and lb < node.startPos + node.read:
                    extractBytes.add(lb - node.startPos + offset)
            offset += node.consumed

        extractBytes = sorted(extractBytes)
        #print "ExtractBytes:", extractBytes

        return extractBytes

    def lookupsUsed(self, lookupWidth, firstLookupAtZero):
        """
        How many lookups are used up by this chain?
        
        Return: lookupsUsed, bytesRemainingInFinalLookup
        """
        # Get a list of lookup bytes across the chain, indexed from the first startPos
        lookupBytes = self.getLookupBytes()
        totalCons = self.totalConsumed()
        totalRead = self.totalRead()

        # Calculate how many lookups are performed and how many unused lookup
        # bytes remain in the final lookup
        unusedLookupBytes = 0
        lookups = 0

        if firstLookupAtZero:
            while len(lookupBytes) > 0 and lookupBytes[0] < lookupWidth:
                lookupBytes.pop(0)
            lookups += 1
            unusedLookupBytes = lookupWidth - totalCons

        while len(lookupBytes) > 0:
            firstLookup = lookupBytes.pop(0)
            while len(lookupBytes) > 0 and lookupBytes[0] < firstLookup + lookupWidth:
                lookupBytes.pop(0)
            lookups += 1
            unusedLookupBytes = lookupWidth - (totalRead - firstLookup)

        if unusedLookupBytes < 0:
            unusedLookupBytes = 0

        #print "%s   Lookups: %d   UnusedLookupBytes: %d" % (self, lookups, unusedLookupBytes)
        #print "FIXME: lookupsUsed"

        return lookups, unusedLookupBytes

    def totalExtracted(self):
        """Number of bytes extracted by the chain"""
        # Get a list of bytes extracted in the chain
        offset = 0
        extractBytes = []
        for node in self.chain:
            nodeExtractBytes = node.dagNode.getExtractBytes()
            for eb in nodeExtractBytes:
                if eb >= node.startPos and eb < node.startPos + node.read:
                    extractBytes.append(eb - node.startPos + offset)
            offset += node.consumed

        #print "ExtractBytes:", extractBytes

        return len(extractBytes)

    def hdrsFound(self):
        """Number of headers found so far in this chain"""
        hdrCount = 0
        for node in self.chain[1:]:
            if isinstance(node, HeaderNode):
                hdrCount += 1

        return hdrCount

    #def getDecisionCombos(self):
    #    """How many decision combinations are there in this chain"""
    #    decCombos = 1
    #    nxtHdrName = None
    #    tarLength = Header.ANY
    #    for pos in xrange(len(self.chain) - 1, -1, -1):
    #        cn = self.chain[pos]
    #        dagNode = cn.dagNode
    #        if isinstance(dagNode, HeaderNode):
    #            decLength = cn.startPos + cn.read
    #            if isinstance(tarLength, int):
    #                tarLength += cn.length
    #            #print "Lookup: %s->%s L: %d" % (dagNode.hdr.name, nxtHdrName, length)
    #            decCombos *= dagNode.hdr.getDecisionCombos(decLength, nxtHdrName, tarLength)
    #            nxtHdrName = dagNode.hdr.name
    #            length = 0
    #        elif isinstance(dagNode, PadNode):
    #            if pos > 0:
    #                tarLength = cn.length
    #            else:
    #                minHdrLen = cn.hdr.length()[0]
    #                tarLength = cn.length + minHdrLen
    #                decCombos *= dagNode.hdr.getDecisionCombos(minHdrLen, nxtHdrName, tarLength)
    #        else:
    #            raise TypeError('Unknown type of dagNode: ' + type(dagNode))
    #        
    #    return decCombos

    def calcPatterns(self):
        """
        Calculate the number of patterns that match the chain
        """
        #print "calcPatterns: %s" % self
        # Identify the header lengths
        totalCons = self.totalConsumed()
        totalRead = self.totalRead()

        cons = totalCons
        read = totalRead
        if len(self.chain) > 0:
            cons += self.chain[0].startPos
            read += self.chain[0].startPos
        consLens = []
        readLens = []
        for node in self.chain:
            #if isinstance(node.dagNode, HeaderNode):
            if node.dagNode.__class__ == HeaderNode:
                maxConsLen = node.dagNode.getLength()
                maxReadLen = node.dagNode.getProcessLength()
            else:
                maxConsLen = node.dagNode.length
                maxReadLen = node.dagNode.length
            consLens.append(min(cons, maxConsLen))
            readLens.append(min(read, maxReadLen))
            cons -= maxConsLen
            read -= maxConsLen # Deliberately CONS not READ
        #print "ConsLens: %s: %s" % (self, consLens)
        #print "ReadLens: %s: %s" % (self, readLens)

        # Work out the number of patterns that match the chain
        patterns = 1
        cons = totalCons
        read = totalRead
        if len(self.chain) > 0:
            cons += self.chain[0].startPos
            read += self.chain[0].startPos
            #if isinstance(self.chain[0].dagNode, SubHeaderNode):
            if self.chain[0].dagNode.__class__ == SubHeaderNode:
                cons += self.chain[0].dagNode.startPos
                read += self.chain[0].dagNode.startPos
            #elif isinstance(self.chain[0].dagNode, BarrierNode):
            elif self.chain[0].dagNode.__class__ == BarrierNode:
                #if len(self.chain) > 1 and isinstance(self.chain[1].dagNode, SubHeaderNode):
                if len(self.chain) > 1 and self.chain[1].dagNode.__class__ == SubHeaderNode:
                    cons += self.chain[1].dagNode.startPos
                    read += self.chain[1].dagNode.startPos
        if self.done:
            nxtHdrName = None
        else:
            nxtHdrName = Header.ANY
        #skip = 0
        seenLastHdr = self.done
        #print "Cons: %d   Read: %d   nxtHdrName: %s" % (cons, read, nxtHdrName)
        pos = len(self.chain) - 1
        while pos >= 0:
        #for pos in xrange(len(self.chain) - 1, -1, -1):
            node = self.chain[pos]
            dagNode = node.dagNode
            hdrNode = dagNode.hdr
            consLen = consLens[pos]
            readLen = readLens[pos]

            #if skip:
            #    skip -= 1
            #    continue

            # Skip over barrier nodes
            #if isinstance(dagNode, BarrierNode):
            if dagNode.__class__ == BarrierNode:
                if pos == len(self.chain) - 1:
                    if node.consumed != 0:
                        raise RuntimeError('Barriers at end of chain should not consume any data: %s' % self)
                    elif pos > 0 and (self.chain[pos-1].dagNode.hdr != hdrNode or self.chain[pos-1].dagNode.inst != dagNode.inst):
                        nxtHdrName = hdrNode.name
                    pos -= 1
                    continue
                elif pos != 0:
                    raise RuntimeError('Barriers may only appear at beginning/end of chain: %s' % self)
                #pos -= 1
                #continue

            #if isinstance(dagNode, SubHeaderNode):
            if dagNode.__class__ == SubHeaderNode:
                prevPos = pos - 1
                while prevPos >= 0 and self.chain[prevPos].dagNode.hdr == hdrNode and self.chain[prevPos].dagNode.inst == dagNode.inst:
                    #consLen += consLens[prevPos]
                    #readLen += readLens[prevPos]
                    prevPos -= 1
                    #skip += 1
                    pos -= 1
                consLen = dagNode.startPos + node.consumed
                readLen = dagNode.startPos + node.read
                if not seenLastHdr:
                    hdrLen = Header.ANY
                else:
                    hdrLen = consLen
                decCombos = hdrNode.getDecisionCombos(
                        min(readLen, read), nxtHdrName, hdrLen)
                #print "Dec combos: %d (%s %s %s %s)" % (decCombos, node, min(readLen, read), nxtHdrName, hdrLen)
                nxtHdrName = hdrNode.name
            #elif isinstance(dagNode, HeaderNode):
            elif dagNode.__class__ == HeaderNode:
                if not seenLastHdr:
                    hdrLen = Header.ANY
                else:
                    hdrLen = consLen
                decCombos = hdrNode.getDecisionCombos(
                        min(readLen, read), nxtHdrName, hdrLen)
                #print "Dec combos: %d (%s %s %s %s)" % (decCombos, node, min(readLen, read), nxtHdrName, hdrLen)
                nxtHdrName = hdrNode.name
            elif dagNode.__class__ == PadNode or dagNode.__class__ == BarrierNode:
                # If we arrive via a BarrierNode then we must be at the
                # beginning of the chain. Also, the barrier node must be at
                # the end of a header.
                decLen = hdrNode.length()[0]/8 + consLen
                #print "Pad:", decLen, nxtHdrName, decLen
                decCombos = hdrNode.getDecisionCombos(
                        decLen, nxtHdrName, decLen)
                #print "Dec combos: %d (%s %s %s %s)" % (decCombos, node, decLen, nxtHdrName, decLen)
                prevPos = pos - 1
                while prevPos >= 0 and self.chain[prevPos].dagNode.hdr == hdrNode and self.chain[prevPos].dagNode.inst == dagNode.inst:
                    consLen += consLens[prevPos]
                    readLen += readLens[prevPos]
                    prevPos -= 1
                    #skip += 1
                    pos -= 1

            cons -= consLen
            read -= consLen
            #print "Dec combos: %d (%s)" % (decCombos, node)
            patterns *= decCombos
            seenLastHdr = True
            #print hdrNode
            #print hdrNode.getDecisionBytes()
            #print hdrNode.getLengthBytes_i()
            #print hdrNode.getLengthFields()
            #print hdrNode.getLengthVarValues()
            #print hdrNode.getDecisionCombos(hdrNode.length())
            pos -= 1

        # Subtract out the initial decision count
        #if not DAGChain.singleEntryForInternal and len(self.chain) > 0:
        #    firstCN = self.chain[0]
        #    dagNode = firstCN.dagNode
        #    hdr = dagNode.hdr
        #    hdrLen = Header.ANY
        #    #if isinstance(dagNode, SubHeaderNode):
        #    if dagNode.__class__ == SubHeaderNode:
        #        prevPos = dagNode.startPos
        #    #elif isinstance(dagNode, BarrierNode):
        #    elif dagNode.__class__ == BarrierNode:
        #        prevPos = dagNode.barrierLoc
        #    #elif isinstance(dagNode, PadNode):
        #    elif dagNode.__class__ == PadNode:
        #        prevPos = hdr.length()[0]
        #        hdrLen = dagNode.getTotalLength()
        #    else:
        #        prevPos = firstCN.startPos
        #    prevDecCombos = hdrNode.getDecisionCombos(prevPos, Header.ANY, hdrLen)
        #    #print "calcPatterns (sub): %s %d %d %d" % (self, patterns, prevDecCombos, patterns / prevDecCombos)
        #    patterns /= prevDecCombos
        #    if patterns < 1:
        #        patterns = 1
        if DAGChain.singleEntryForInternal and len(self.chain) > 0:
            firstCN = self.chain[0]
            firstDagNode = firstCN.dagNode
            firstHdr = firstDagNode.hdr

            lastCN = self.chain[-1]
            lastDagNode = lastCN.dagNode
            lastHdr = lastDagNode.hdr

            if firstHdr == lastHdr and firstDagNode.inst == lastDagNode.inst:
                hdrLen = firstHdr.length()[0]

                if firstDagNode.__class__ == SubHeaderNode:
                    startPos = firstDagNode.startPos + firstCN.startPos
                elif firstDagNode.__class__ == BarrierNode:
                    startPos = firstDagNode.barrierLoc
                elif firstDagNode.__class__ == PadNode:
                    startPos = firstHdr.length()[0] + firstCN.startPos
                    hdrLen = firstDagNode.getTotalLength()
                else:
                    startPos = firstCN.startPos

                if lastDagNode.__class__ == SubHeaderNode:
                    endPos = lastDagNode.startPos + lastCN.startPos + lastCN.consumed - 1
                elif lastDagNode.__class__ == BarrierNode:
                    endPos = lastDagNode.barrierLoc - 1
                elif lastDagNode.__class__ == PadNode:
                    endPos = firstHdr.length()[0] + lastCN.startPos + lastCN.consumed - 1
                    hdrLen = lastDagNode.getTotalLength()
                else:
                    endPos = lastCN.startPos + lastCN.consumed - 1

                #decBytes = firstHdr.getDecisionBytes()[0]
                #while len(decBytes) > 0 and decBytes[0] < startPos:
                #    decBytes.pop(0)
                #while len(decBytes) > 0 and decBytes[-1] > endPos:
                #    decBytes.pop()
                decBefore = firstHdr.getDecisionCombos(startPos, Header.ANY, Header.ANY)
                decAfter = firstHdr.getDecisionCombos(endPos + 1, Header.ANY, Header.ANY)
                #if len(decBytes) == 0:
                if decBefore == decAfter:
                    #print "calcPatterns (sefi): %s   Patterns: %d -> %d" % (self, patterns, 1)
                    patterns = 1

        self.patterns = patterns
        #print "calcPatterns (Done): %s" % (self)
        #print ""

    #def findLookupSubchains(self):
    #    """
    #    Identify all 'sensible' lookup subchains.

    #    A lookup subchain is a subchain that ends in a lookup, or contains no lookups.
    #    """
    #    lookupBytes = self.getLookupBytes()

    #    subchains = []
    #    if len(lookupBytes) == 0:
    #        return [self]
    #    else:
    #        #lookupBytes.append(lookupBytes[-1] + 1)
    #        if lookupBytes[0] > 0:
    #            lookupBytes.insert(0, lookupBytes[0] - 1)
    #        for lb in lookupBytes:
    #            print lb
    #            newChain = self.copyToDepth(lb+1)
    #            newChain.done = False
    #            newChain.calcPatterns()
    #            subchains.append(newChain)

    #    return subchains

    def findLookupSubchains(self):
        """
        Identify all 'sensible' lookup subchains.

        A lookup subchain is a subchain that ends in a lookup, or contains no lookups.
        Lookups are generated for each lookup byte in the header (and for zero
        lookup bytes for the first header)
        """
        lookupBytes = self.getLookupBytes()

        #print "findLookupSubchains:", self, len(self.chain), lookupBytes
        #if len(self.chain) > 0:
        #    print type(self.chain[0].dagNode)

        subchains = []
        # If we start with a barrier node, then we need to add subchains
        # for every chain up until the first containing a decision byte
        if len(self.chain) > 0: # and isinstance(self.chain[0].dagNode, BarrierNode):
            for i in xrange(0, len(self.chain)):
                cn = self.chain[i]
                cnLookupBytes = cn.dagNode.getDecisionBytes()
                cnLookupBytes = [lb - self.chain[i].startPos for lb in cnLookupBytes]
                while len(cnLookupBytes) > 0 and cnLookupBytes[0] < 0:
                    cnLookupBytes.pop(0)
                while len(cnLookupBytes) > 0 and cnLookupBytes[-1] >= self.chain[i].read:
                    cnLookupBytes.pop()
                if len(cnLookupBytes) == 0:
                    newChain = self.copy()
                    newChain.chain = newChain.chain[0:i + 1]
                    newChain.done = False
                    newChain.calcPatterns()
                    subchains.append(newChain)
                else:
                    break
            if len(subchains) > 1 and \
                    self.chain[0].dagNode.hdr == self.chain[1].dagNode.hdr and \
                    self.chain[0].dagNode.inst == self.chain[1].dagNode.inst and \
                    isinstance(self.chain[0].dagNode, BarrierNode) and \
                    not isinstance(self.chain[1].dagNode, PadNode):
                subchains.pop(0)

        if len(lookupBytes) == 0:
            ## We should return every subsequence of headers
            #for i in xrange(len(self.chain)):
            #    newChain = self.copy()
            #    newChain.chain = newChain.chain[0:i + 1]
            #    newChain.done = False
            #    newChain.calcPatterns()
            #    subchains.append(newChain)
            #return subchains
            #return [self]
            if len(subchains) == 0:
                subchains.append(self)
            #print "Subchains:"
            #for s in subchains:
            #    print "  %s" % s
            return subchains
        else:
            #print "findLookupSubchains:", self

            seenLookupBytes = False
            firstDone = False
            for i in xrange(len(self.chain)):
                cn = self.chain[i]
                cnLookupBytes = cn.dagNode.getDecisionBytes()
                cnLookupBytes = [lb - self.chain[i].startPos for lb in cnLookupBytes]
                while len(cnLookupBytes) > 0 and cnLookupBytes[0] < 0:
                    cnLookupBytes.pop(0)
                while len(cnLookupBytes) > 0 and cnLookupBytes[-1] >= self.chain[i].read:
                    cnLookupBytes.pop()
                if len(cnLookupBytes) == 0:
                    continue

                if not firstDone and cnLookupBytes[0] > 0:
                    newChain = self.copy()
                    newChain.chain = newChain.chain[0:i + 1]
                    newChain.done = False
                    lastCN = newChain.chain[-1]
                    lastCN.consumed -= lastCN.consumed - cnLookupBytes[0]
                    lastCN.read = lastCN.consumed
                    newChain.calcPatterns()
                    subchains.append(newChain)
                    #print "a", newChain, self

                for lookup in cnLookupBytes:
                    newChain = self.copy()
                    newChain.chain = newChain.chain[0:i + 1]
                    newChain.done = False
                    lastCN = newChain.chain[-1]
                    if lookup + 1 < lastCN.consumed:
                        lastCN.consumed -= lastCN.consumed - lookup - 1
                    if lookup + 1 < lastCN.read:
                        lastCN.read -= lastCN.read - lookup - 1
                    newChain.calcPatterns()
                    subchains.append(newChain)
                    #print "b", lookup + lastCN.startPos, newChain, self

                firstDone = True

            ##lookupBytes.append(lookupBytes[-1] + 1)
            #if lookupBytes[0] > 0:
            #    lookupBytes.insert(0, lookupBytes[0] - 1)
            #for lb in lookupBytes:
            #    print lb
            #    newChain = self.copyToDepth(lb+1)
            #    newChain.done = False
            #    newChain.calcPatterns()
            #    subchains.append(newChain)

        #print "Subchains:"
        #for s in subchains:
        #    print "  %s" % s

        return subchains

    ##### def findLookupSubchains(self):
    #####     """
    #####     Identify all 'sensible' lookup subchains.

    #####     A lookup subchain is a subchain that ends in a lookup, or contains no lookups.
    #####     Lookups span the maximum number of lookup bytes in each header.
    #####     """
    #####     lookupBytes = self.getLookupBytes()

    #####     subchains = []
    #####     if len(lookupBytes) == 0:
    #####         return [self]
    #####     else:
    #####         #print "findLookupSubchains:", self
    #####         for i in xrange(len(self.chain)):
    #####             cn = self.chain[i]
    #####             cnLookupBytes = cn.dagNode.getDecisionBytes()
    #####             cnLookupBytes = [lb - self.chain[i].startPos for lb in cnLookupBytes]
    #####             while len(cnLookupBytes) > 0 and cnLookupBytes[0] < 0:
    #####                 cnLookupBytes.pop(0)
    #####             while len(cnLookupBytes) > 0 and cnLookupBytes[-1] >= self.chain[i].read:
    #####                 cnLookupBytes.pop()
    #####             if len(cnLookupBytes) == 0:
    #####                 continue

    #####             lastLookup = cnLookupBytes[-1]
    #####             newChain = self.copy()
    #####             newChain.chain = newChain.chain[0:i + 1]
    #####             newChain.done = False
    #####             lastCN = newChain.chain[-1]
    #####             if lastLookup + 1 < lastCN.consumed:
    #####                 lastCN.consumed -= lastCN.consumed - lastLookup - 1
    #####             if lastLookup + 1 < lastCN.read:
    #####                 lastCN.read -= lastCN.read - lastLookup - 1
    #####             newChain.calcPatterns()
    #####             subchains.append(newChain)
    #####             #print "b", lookup + lastCN.startPos, newChain, self


    #####         ##lookupBytes.append(lookupBytes[-1] + 1)
    #####         #if lookupBytes[0] > 0:
    #####         #    lookupBytes.insert(0, lookupBytes[0] - 1)
    #####         #for lb in lookupBytes:
    #####         #    print lb
    #####         #    newChain = self.copyToDepth(lb+1)
    #####         #    newChain.done = False
    #####         #    newChain.calcPatterns()
    #####         #    subchains.append(newChain)

    #####     #print "Subchains:"
    #####     #for s in subchains:
    #####     #    print "  %s" % s

    #####     return subchains

    def successors(self):
        """Identify all successor nodes from this node"""
        if self._successors:
            return self._successors

        if self.done:
            return []

        self._successors = set()
        last = self.chain[-1]
        if last.unread() > 0 or isinstance(last.dagNode, BarrierNode):
            self._successors.add(DAGChainNode(last.dagNode, last.startPos +
                last.consumed, 0, 0))
        else:
            for nxt in last.dagNode.nxt:
                if nxt:
                    self._successors.add(DAGChainNode(nxt, 0, 0, 0))
                else:
                    self._successors.add(nxt)


        return self._successors

    def ingress(self):
        """
        Identify an ingress node for this chain (corresponding to another
        chain's successor)
        """
        if len(self.chain) == 0:
            return None

        first = self.chain[0]
        return DAGChainNode(first.dagNode, first.startPos, 0, 0)

    def getLookupByteValues(self, trim=True):
        """Return all sets of lookup byte values"""
        print "getLookupByteValues:", self, trim
        #lookupBytes = ['']
        offset = 0
        prevOffset = 0
        pos = 0
        lookupBytes = []
        lookupValues = [([], [])]
        while pos < len(self.chain):
            cnode = self.chain[pos]
            dagNode = cnode.dagNode
            hdr = dagNode.hdr
            inst = dagNode.inst
            if isinstance(dagNode, SubHeaderNode):
                startPos = cnode.startPos + dagNode.startPos
                offset = -startPos + prevOffset
                consdepth = startPos + cnode.consumed
                readdepth = startPos + cnode.read
                lengths = hdr.getHeaderLengths()[1]
                if len(lengths) == 1:
                    length = lengths[0]
                else:
                    length = Header.ANY
            elif isinstance(dagNode, HeaderNode):
                startPos = cnode.startPos
                offset = -startPos + prevOffset
                consdepth = startPos + cnode.consumed
                readdepth = startPos + cnode.read
                lengths = dagNode.getLengths()
                if len(lengths) == 1:
                    length = lengths[0]
                else:
                    length = Header.ANY
            elif isinstance(dagNode, PadNode):
                startPos = cnode.startPos + hdr.length()[0]
                offset = -startPos + prevOffset
                consdepth = startPos + cnode.consumed
                readdepth = startPos + cnode.read
                length = dagNode.getTotalLength()
            elif isinstance(dagNode, BarrierNode):
                pos += 1
                continue
            else:
                raise TypeError('Unknown type: %s' % type(dagNode))
            #print dagNode, startPos, offset, consdepth, readdepth, length
            depth = consdepth

            nxtPos = pos + 1
            while nxtPos < len(self.chain):
                nxtCNode = self.chain[nxtPos]
                nxtDAGNode = nxtCNode.dagNode
                nxtHdr = nxtDAGNode.hdr
                nxtInst = nxtDAGNode.inst
                if nxtHdr != hdr or nxtInst != inst:
                    nxtHdr = nxtHdr.name
                    break
                elif isinstance(nxtDAGNode, PadNode):
                    length = nxtDAGNode.getTotalLength()
                    consdepth = length
                    readdepth = length
                    depth = length
                elif isinstance(nxtDAGNode, HeaderNode) or \
                        isinstance(nxtDAGNode, SubHeaderNode):
                    readdepth = consdepth + nxtCNode.read
                    consdepth += nxtCNode.consumed
                    depth += nxtCNode.consumed
                nxtPos += 1

            if nxtPos >= len(self.chain):
                nxtHdr = Header.ANY
                successors = self.successors()
                if len(successors) > 0:
                    successor = successors.pop()
                    if successor.dagNode.hdr != hdr or successor.dagNode.inst != inst:
                        nxtHdr = successor.dagNode.hdr.name
                else:
                    nxtHdr = None

            #prevHdrLookupBytes = hdr.getDecisionComboBytes(startPos, nxtHdr, length)
            #hdrLookupBytes = hdr.getDecisionComboBytes(depth, nxtHdr, length)
            (prevHdrLookupBytes, prevHdrLookupMatch) = hdr.getDecisionComboBytes2(startPos, nxtHdr, length)
            (hdrLookupBytes, hdrLookupValues) = hdr.getDecisionComboBytes2(readdepth, nxtHdr, length)
            #print hdr.name, nxtHdr, startPos, depth, length, prevHdrLookupBytes, hdrLookupBytes
            #print mphlb, mhlb, hlv

            #if trim:
            #    if len(prevHdrLookupBytes) > 0:
            #        drop = len(prevHdrLookupBytes.pop())
            #    else:
            #        drop = 0
            #    hdrLookupBytes = set([lb[drop:] for lb in hdrLookupBytes])

            #if len(hdrLookupBytes) > 0:
            #    newLookupBytes = []
            #    for lb in lookupBytes:
            #        for hlb in hdrLookupBytes:
            #            newLookupBytes.append(lb + hlb)
            #    lookupBytes = newLookupBytes

            prevHdrLookupBytes = map(lambda x : x + offset, prevHdrLookupBytes)
            hdrLookupBytes = map(lambda x : x + offset, hdrLookupBytes)
            lookupBytes.extend(hdrLookupBytes)
            if len(hdrLookupValues) > 0:
                newLookupValues = []
                for lv in lookupValues:
                    for hlv in hdrLookupValues:
                        lvcopy = copy.deepcopy(lv)
                        lvcopy[0].extend(hlv[0])
                        lvcopy[1].extend(hlv[1])
                        newLookupValues.append(lvcopy)
                lookupValues = newLookupValues
            #print lm2

            pos = nxtPos
            prevOffset = offset + consdepth

        # Assemble the lookup bytes
        if len(lookupBytes) > 0:
            lbMin = min(lookupBytes)
            lbMax = max(lookupBytes)
            offset = 0
            if lbMin < 0 and not trim:
                offset = -lbMin
            lbMax += offset

            sortedLookupBytes = sorted(set(lookupBytes))
            newLookupValues = []
            for lv in lookupValues:
                mask = [0 for x in xrange(lbMax + 1)]
                match = [0 for x in xrange(lbMax + 1)]

                for i in xrange(len(lookupBytes)):
                    dest = lookupBytes[i]
                    dest += offset
                    if offset >= 0:
                        mask[dest] |= lv[0][i]
                        match[dest] |= lv[1][i]

                matchStr = ""
                for src in sortedLookupBytes:
                    src += offset
                    if offset >= 0:
                        matchStr += "%02x%02x" % (mask[src], match[src])
                newLookupValues.append(matchStr)
                
            lookupValues = newLookupValues
        else:
            lookupValues = ['']

        print "getLookupByteValues (DONE):", self, lookupValues
        #print lb2, lm2



        return lookupValues

# Basic test code
if __name__ == '__main__':
    hdr = Header.Header('TestHeader')
    hdr.addField('test', 32)
    chainNode1 = DAGChainNode(HeaderNode(hdr, 1), 1, 3, 3)
    chainNode2 = DAGChainNode(HeaderNode(hdr, 2), 0, 2, 2)

    chain = DAGChain()
    chain.add(chainNode1)
    chain.add(chainNode2)

    print chain
