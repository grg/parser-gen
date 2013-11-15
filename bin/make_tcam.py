#!/usr/bin/env python

from Header import Header, ANY
from Field import Field
import string
import sys
import re
import copy
import math
from HeaderLib import exploreHeader, readHeaders
from DAGHeaderNode import HeaderNode
from DAGSubHeaderNode import SubHeaderNode
from DAGBarrierNode import BarrierNode
from DAGPadNode import PadNode
from DAGChain import DAGChain
from DAGChainNode import DAGChainNode
from OptNode import OptNode
import time
#from LookupChain import LookupChain
import argparse
import cProfile
from random import randint, random
from OptContext import OptContext

dataRate = 10
clkFreq = 1

# Bits per cycle
globalBPC = dataRate / clkFreq

# Variables identified by Glen

lookups = 2
lookupWidth = 2
maxSkip = 128
minSkip = 0
firstLookupAtZero = True
windowSize = 30
extract = False
extractBytes = 0
maxHdrs = 3

showResultSets = False
multParent = True
multParentRetry = True
parallelEdge = True
instMerge = True
ternMatchOnState = True
debug = False
tcamMaxState = 255
printTCAM = True
saveTCAM = False

hfile = 'headers.txt'

ccontext = OptContext()

resultsRet = 0
resultsTerminate = 14

extractPos = {}
extractOffset = {}
extractMaxDest = 0


def opt(context, cnode, bpc):
    """
    Attempt to find the optimal TCAM allocation for a DAG rooted at a given
    node and starting at position pos. (Assume decision bytes prior to pos have
    been processed.)

    Return:
      [(lookupClusters, worstByteCount, worstCyc)]
    """
    #print "opt: N=%s, BPC=%d, C=%s" % (cnode, bpc, context)

    optPair = (cnode, bpc)

    if optPair in context.results:
        return context.results[optPair]

    context.count += 1

    startCN = cnode
    while isinstance(startCN.dagNode, BarrierNode) and len(startCN.dagNode.nxt) == 1 and list(startCN.dagNode.nxt)[0] is not None:
        startCN = DAGChainNode(list(startCN.dagNode.nxt)[0], 0, 0, 0)

    minEdgeCnt = None
    minCluster = None
    coveredClusters = findClustersAndCovers(context, startCN)
    for cluster in sorted(coveredClusters):
        consumed = cluster.totalConsumed()
        edgeCount = findEdgeCount(context, cluster)
        baseEdge = edgeCount
        #fringe = findFringe(context, cluster)
        doneBeatsBPC = True
        coverBeatsMinSkip = True
        #print cluster
        #print "  Consumed: %d" % consumed
        #print "  Edge count: %d" % edgeCount
        for cover in coveredClusters[cluster]:
            #print "  Cover: %s   Cons: %d b / %d B   BPC: %d" % (cover, cover.totalConsumed(), cover.totalConsumed() * 8, bpc)
            if cover.done:
                if cover.totalConsumed() * 8 < bpc:
                    doneBeatsBPC = False
            elif cover.totalConsumed() < minSkip:
                coverBeatsMinSkip = False
        if not doneBeatsBPC or not coverBeatsMinSkip:
            continue
        fringe = findFringeWithMinCoverLen(context, cluster, coveredClusters)
        #print "  Fringe:"
        #for f in fringe:
        #    print "    %s: %d" % (f, fringe[f])
        fringeBeatsBPC = True
        for f in fringe:
            if fringeBeatsBPC:
                targetBPC = globalBPC + (bpc - fringe[f] * 8)
                #print "    opt: %s, %d" % (f, targetBPC)
                (fEdgeCnt, fCluster) = opt(context, f, targetBPC)
                if not fEdgeCnt:
                    fringeBeatsBPC = False
                else:
                    edgeCount += fEdgeCnt

        if fringeBeatsBPC:
            if not minEdgeCnt or edgeCount < minEdgeCnt:
                minEdgeCnt = edgeCount
                minCluster = cluster
        #print "CHOICE: (%s, %s): %d   %s %d" % (cnode, bpc, edgeCount, cluster, baseEdge)

    #print "MIN: (%s, %s): %d   %s %d" % (cnode, bpc, minEdgeCnt, minCluster, minBaseEdge)

    context.results[optPair] = (minEdgeCnt, minCluster)

    return context.results[optPair]

def findOptNodes(context, cnode, bpc):
    optPair = (cnode, bpc)

    (minEdgeCnt, minCluster) = context.results[(cnode, bpc)]
    if not minCluster:
        return (max(windowSize, maxSkip) + 1, 1)

    startCN = cnode
    while isinstance(startCN.dagNode, BarrierNode) and len(startCN.dagNode.nxt) == 1 and list(startCN.dagNode.nxt)[0] is not None:
        startCN = DAGChainNode(list(startCN.dagNode.nxt)[0], 0, 0, 0)

    worstBits = 8 * windowSize
    worstCyc = 1
    worstBPC = worstBits * 1.0 / worstCyc
    coveredClusters = findClustersAndCovers(context, startCN)

    optChoice = (cnode, minCluster)
    if optChoice not in context.optNodes:

        # Attempt to identify the worst BPC for this chain
        for cover in coveredClusters[minCluster]:
            if cover.done:
                coverBits = cover.totalConsumed() * 8
                if coverBits < worstBPC:
                    worstBits = coverBits
                    worstCyc = 1
                    worstBPC = worstBits * 1.0 / worstCyc

        context.optNodes.add(optChoice)
        fringe = findFringeWithMinCoverLen(context, minCluster)
        for f in fringe:
            targetBPC = globalBPC + (bpc - fringe[f] * 8)
            (subBits, subCyc) = findOptNodes(context, f, targetBPC)
            subBits += fringe[f] * 8
            subCyc += 1
            subBPC = subBits * 1.0 / subCyc
            if subBPC < worstBPC:
                worstBits = subBits
                worstCyc = subCyc
                worstBPC = subBPC

    context.worstBits = worstBits
    context.worstCyc = worstCyc
    context.worstBPC = worstBPC

    context.bestClusters = set()
    context.bestEdgeCount = 0
    for (cnode, cluster) in context.optNodes:
        if cluster not in context.bestClusters:
            context.bestEdgeCount += findEdgeCount(context, cluster)
        context.bestClusters.add(cluster)

    return (worstBits, worstCyc)

def optFound(context):
    return len(context.optNodes) > 0

def printBestOpt(context, printEdges=True):
    if len(context.optNodes) > 0:
        print "opt-algorithm best edge count: %d" % context.bestEdgeCount
        print "opt-algorithm worst bits-per-cycle: %1.3f   (%1.3f bytes-per-cycle)" % \
                (context.worstBPC, context.worstBPC / 8)
        if printEdges:
            print "opt-algorithm optimal clusters:"
            for cluster in sorted(context.bestClusters):
                print "  %s   (%d edges)" % (cluster, findEdgeCount(context, cluster))
                #cnode = DAGChainNode(cluster.chain[0].dagNode, cluster.chain[0].startPos, 0, 0)
                #coveredClusters = findClustersAndCovers(cnode)
                #fringe = findFringeWithMinCoverLen(cluster, coveredClusters)
                #print "    ",
                #for f in sorted(fringe):
                #    print f,
                #print ""
    else:
        print "opt-algorithm could not find optimal that met required BPC (%d) or minimum skip amount (%d)" % (globalBPC, minSkip)

def printEntries(context):
    if len(context.optNodes) > 0:
        print "opt-algorithm TCAM entries:"
        for cluster in sorted(context.bestClusters):
            print "  %s   (%d edges)" % (cluster, findEdgeCount(context, cluster))
            cnode = cluster.chain[0]
            coveredClusters = findClustersAndCovers(context, cnode)
            for entry in sorted(coveredClusters[cluster]):
                print "     %s" % entry

def numStatesNeeded(chain):
    """Identify the number of states needed by a chain"""
    if len(chain.chain) == 0:
        return 0

    firstCNode = chain.chain[0]
    firstDAGNode = firstCNode.dagNode
    firstHdr = firstDAGNode.hdr
    firstInst = firstDAGNode.inst
    firstPos = firstCNode.startPos
    firstLen = ANY
    if isinstance(firstDAGNode, BarrierNode):
        firstPos += firstDAGNode.barrierLoc
    elif isinstance(firstDAGNode, PadNode):
        firstPos += firstHdr.length()[0]
        firstLen = firstDAGNode.getTotalLength()
    elif isinstance(firstDAGNode, SubHeaderNode):
        firstPos += firstDAGNode.startPos

    return firstHdr.getDecisionCombos(firstPos, ANY, firstLen)

def calcNumStatesNeeded(context):
    context.numStatesNeeded = {}
    for cluster in context.bestClusters:
        context.numStatesNeeded[cluster] = numStatesNeeded(cluster)

def isNextNodeInSameHdr(context, chain):
    """
    Check whether any next node is within the same header
    """
    if len(chain.chain) == 0:
        return False

    firstHdr = chain.chain[0].dagNode.hdr
    firstInst = chain.chain[0].dagNode.inst

    fringe = findFringe(context, chain)
    for f in fringe:
        fringeHdr = f.dagNode.hdr
        fringeInst = f.dagNode.inst

        if fringeHdr == firstHdr and fringeInst == firstInst:
            return True

    return False

def isSuccessorInSameHdr(context, chain, successor):
    """
    Check whether any next node is within the same header
    """
    if len(chain.chain) == 0:
        return False

    if successor is None:
        return False

    firstHdr = chain.chain[0].dagNode.hdr
    firstInst = chain.chain[0].dagNode.inst

    successorHdr = successor.dagNode.hdr
    successorInst = successor.dagNode.inst

    return successorHdr == firstHdr and successorInst == firstInst

def findPrecedingMatchBytes(chain):
    firstCNode = chain.chain[0]
    firstDAGNode = firstCNode.dagNode
    firstHdr = firstDAGNode.hdr
    firstInst = firstDAGNode.inst
    firstPos = firstCNode.startPos
    firstLen = ANY
    if isinstance(firstDAGNode, BarrierNode):
        pos = 1
        while isinstance (firstDAGNode, BarrierNode) and pos < len(chain.chain):
            firstDAGNode = chain.chain[pos].dagNode
            pos += 1

    if isinstance(firstDAGNode, PadNode):
        firstPos += firstHdr.length()[0]
        firstLen = firstDAGNode.getTotalLength()
    elif isinstance(firstDAGNode, SubHeaderNode):
        firstPos += firstDAGNode.startPos

    matchBytes = firstHdr.getDecisionComboBytes(firstPos, ANY, firstLen)

    if len(matchBytes) == 0:
        matchBytes.append('')

    return matchBytes

def cmpMaskMatchArray(a, b):
    aMaskArray = a[0]
    aMatchArray = a[1]
    bMaskArray = b[0]
    bMatchArray = b[1]

    aLen = len(aMaskArray)
    bLen = len(bMaskArray)
    minLen = min(aLen, bLen)
    for i in xrange(minLen):
        aMask = aMaskArray[i]
        bMask = bMaskArray[i]
        aMatch = aMatchArray[i]
        bMatch = bMatchArray[i]

        aMaskBits = 0
        bMaskBits = 0
        for bit in xrange(8):
            if aMask & (2 ** bit) != 0:
                aMaskBits += 1
            if bMask & (2 ** bit) != 0:
                bMaskBits += 1

        if aMaskBits > bMaskBits:
            return -1
        elif bMaskBits > aMaskBits:
            return 1

        if aMatch > bMatch:
            return 1
        elif bMatch > aMatch:
            return -1

    if aLen > bLen:
        return 1
    elif bLen > aLen:
        return -1

    return 0

def cmpMatchBytesStr(a, b):
    aByteArray = map(ord, a.decode('hex'))
    aMaskArray = aByteArray[0::2]
    aMatchArray = aByteArray[1::2]

    bByteArray = map(ord, b.decode('hex'))
    bMaskArray = bByteArray[0::2]
    bMatchArray = bByteArray[1::2]

    return cmpMaskMatchArray((aMaskArray, aMatchArray), (bMaskArray, bMatchArray))

def cmpEntry(a, b):
    c = cmpMaskMatchArray(a, b)
    if c != 0:
        return c
    
    c = cmp(a[3], b[3])
    if c != 0:
        return c

    return cmp(a[2], b[2])

def allocateState(context, chain, availableState):
    ingress = chain.ingress()
    statesNeeded = numStatesNeeded(chain)
    pow2StatesNeeded = int(2 ** math.ceil(math.log(statesNeeded, 2)))
    startPos = 0
    endPos = 0
    while endPos <= len(availableState):
        arraySpan = endPos - startPos + 1
        span = availableState[endPos] - availableState[startPos] + 1
        if availableState[startPos] % pow2StatesNeeded != 0:
            startPos += 1
            endPos = startPos
        elif arraySpan == statesNeeded and \
                span == statesNeeded:
            matchBytes = sorted(findPrecedingMatchBytes(chain), cmpMatchBytesStr)
            if len(matchBytes) > 0 and len(matchBytes[0]) > 0:
                mask = map(ord, matchBytes[0].decode('hex'))[0::2]
                for mb in matchBytes[1:]:
                    newMask = map(ord, mb.decode('hex'))[0::2]
                    for i in xrange(len(mask)):
                        mask[i] |= newMask[i]
            else:
                mask = []

            context.allocation[ingress] = availableState[startPos]
            context.matchBytes[ingress] = matchBytes
            context.matchMask[ingress] = mask
            print "CHAIN:", chain, matchBytes
            del availableState[startPos:endPos + 1]
            return True
        elif span != arraySpan:
            startPos = endPos
        elif span < statesNeeded:
            endPos += 1

    return False

def calcExtractLocs(hdrList, hdrs):
    # Create a sorted list of header names

    global extractMaxDest, extractPos, extractOffset

    firstHdr = hdrList[0]
    firstHdrName = firstHdr.name

    # Walk through all possible sequences of headers to identify the maximum count of each header type
    maxCount = {}
    stack = []
    stack.append((firstHdrName, 0, []))
    while len(stack) > 0:
        (hdrName, nxtIndex, seq) = stack.pop()
        hdr = hdrs[hdrName]

        # Verify we aren't exceeding the reference count
        if hdr.refCount:
            maxCnt = hdr.refCount.maxVal
            maxName = hdr.refCount.name
            cnt = 1
            for prevHdrName in seq:
                prevHdr = hdrs[prevHdrName]
                if prevHdr.refCount and prevHdr.refCount.name == maxName:
                    cnt += 1
            if cnt > maxCnt:
                continue

        endSeq = False
        if hdr.nextHeader and type(hdr.nextHeader) == tuple:
            if nxtIndex < len(hdr.nextHeader[1]):
                stack.append((hdrName, nxtIndex + 1, seq))
                nxtHdrName = hdr.nextHeader[1][nxtIndex][1]

                if nxtHdrName:
                    nxtHdr = hdrs[nxtHdrName]
                    newSeq = copy.copy(seq)
                    newSeq.append(hdrName)
                    stack.append((nxtHdrName, 0, newSeq))
                else:
                    endSeq = True
        else:
            endSeq = True

        if endSeq:
            newSeq = copy.copy(seq)
            newSeq.append(hdrName)

            hdrCounts = {}
            for h in newSeq:
                if h not in hdrCounts:
                    hdrCounts[h] = 0
                hdrCounts[h] += 1

            for h in hdrCounts:
                cnt = hdrCounts[h]
                if h not in maxCount:
                    maxCount[h] = 0
                if cnt > maxCount[h]:
                    maxCount[h] = cnt
    


    # Create a sorted list of header names, sorted by order of appearance
    # (breadth-first search)
    seenHdrs = set()
    hdrNames = []
    pending = []

    pending.append(firstHdrName)
    seenHdrs.add(firstHdrName)

    while len(pending) > 0:
        hdrName = pending.pop(0)
        hdr = hdrs[hdrName]
        hdrNames.append(hdrName)

        nxtHdrs = []
        if hdr.nextHeader and type(hdr.nextHeader) == tuple:
            for (a, nxtHdrName) in hdr.nextHeader[1]:
                if nxtHdrName:
                    if nxtHdrName not in nxtHdrs:
                        nxtHdrs.append(nxtHdrName)
            nxtHdrs.sort()
            for nxtHdr in nxtHdrs:
                if nxtHdr not in seenHdrs:
                    pending.append(nxtHdr)
                    seenHdrs.add(nxtHdr)

    # Create tables of header extract locations/offsets
    currPos = 1
    extractPos = {}
    extractOffset = {}
    for hdrName in hdrNames:
        hdr = hdrs[hdrName]
        numExtBytes = len(hdr.getExtractBytes())
        numInst = maxCount[hdrName]

        extractPos[hdrName] = currPos
        extractOffset[hdrName] = numExtBytes

        currPos += numExtBytes * numInst

    extractMaxDest = currPos


def printTCAMEntry(context, chain):
    """Print a single TCAM entry"""
    # Work out where to put each lookup byte
    lookupBytes = chain.getLookupBytes()
    lookupMap = []
    lookup = 0
    lookupStart = 0
    if firstLookupAtZero:
        while len(lookupBytes) > 0 and lookupBytes[0] - lookupStart < lookupWidth:
            lookupPos = lookupBytes.pop(0)
            lookupMap.append(lookup * lookupWidth + lookupPos - lookupStart)
        lookup += 1
    while len(lookupBytes) > 0:
        lookupStart = lookupBytes[0]
        while len(lookupBytes) > 0 and lookupBytes[0] - lookupStart < lookupWidth:
            lookupPos = lookupBytes.pop(0)
            lookupMap.append(lookup * lookupWidth + lookupPos - lookupStart)
        lookup += 1

    stateFmtStr = '%%%dd/%%%dd' % (stateWidth10, stateWidth10)
    srcWidth = int(math.ceil(math.log10(windowSize)))
    extractDstWidth = int(math.ceil(math.log10(context.maxDest)))
    extractFmtStr = '(%%%dd, %%%dd)' % (srcWidth, extractDstWidth)
    hdrNumWidth = int(math.ceil(math.log10(context.hdrCount)))
    maxDestWidth = int(math.ceil(math.log10(context.maxDest)))
    hdrStartFmtStr = '(%%%dd, %%%dd, %%%dd)' % (srcWidth, hdrNumWidth, maxDestWidth)
    advWidth = int(math.ceil(math.log10(maxSkip)))
    advFmtStr = '%%%dd' % (advWidth)

    lookupBytes = chain.getLookupBytes()

    ingress = chain.ingress()
    firstCNode = chain.chain[0]
    firstDAGNode = firstCNode.dagNode
    firstHdr = firstDAGNode.hdr
    firstInst = firstDAGNode.inst
    firstPos = firstCNode.startPos
    firstLen = ANY
    if isinstance(firstDAGNode, BarrierNode):
        firstPos += firstDAGNode.startPos
    elif isinstance(firstDAGNode, PadNode):
        firstPos += firstHdr.length()[0]
        firstLen = firstDAGNode.getTotalLength()

    baseState = context.allocation[ingress]
    baseMatchBytes = context.matchBytes[ingress]
    if len(baseMatchBytes) > 0:
        baseMatchByteLen = len(baseMatchBytes[0])
    else:
        baseMatchByteLen = 0

    # Walk through each cover and assign entries
    cnode = chain.chain[0]
    covers = findClustersAndCovers(context, cnode)[chain]
    entries = []
    for cover in sorted(covers):
        print chain, cover

        lookupValues = cover.getLookupByteValues(False)
        if len(lookupValues) == 0:
            lookupValues.add('')

        successors = cover.successors()
        if len(successors) > 0:
            successor = successors.pop()
            while isinstance(successor.dagNode, BarrierNode):
                successorDAGNode = list(successor.dagNode.nxt)[0]
                successor = DAGChainNode(successorDAGNode, 0, 0, 0)
            if instMerge:
                successor = firstInstCNode(successor)
            #print cover, successor
            successorState = context.allocation[successor]
            successorMatchBytes = context.matchBytes[successor]
            successorMatchMask = context.matchMask[successor]
            if len(successorMatchBytes) > 0:
                successorMatchByteLen = len(successorMatchBytes[0])
            else:
                successorMatchByteLen = 0

            nxtLookupBytes = []
            for nxtChain in context.bestClusters:
                firstCNode = nxtChain.chain[0]
                if firstCNode.dagNode == successor.dagNode and firstCNode.startPos == successor.startPos:
                    nxtLookupBytes = nxtChain.getLookupBytes()
                    break
            nxtLookupStarts = []
            if firstLookupAtZero:
                nxtLookupStarts.append(0)
                while len(nxtLookupBytes) > 0 and nxtLookupBytes[0] < lookupWidth:
                    nxtLookupBytes.pop(0)
            while len(nxtLookupBytes) > 0:
                nxtLookupStarts.append(nxtLookupBytes.pop(0))
                while len(nxtLookupBytes) > 0 and nxtLookupBytes[0] < lookupWidth + nxtLookupStarts[-1]:
                    nxtLookupBytes.pop(0)
        else:
            successor = None
            successorState = tcamMaxState
            successorMatchBytes = ['']
            successorMatchMask = []
            successorMatchByteLen = 0
            nxtLookupStarts = []
        while len(nxtLookupStarts) < lookups:
            nxtLookupStarts.append(0)

        skip = cover.totalConsumed()
        for lookupValueFull in sorted(lookupValues, cmpMatchBytesStr, reverse=True):
            stateInBytes = lookupValueFull[0:baseMatchByteLen]
            if successorMatchByteLen > 0:
                stateOutBytes = lookupValueFull[-successorMatchByteLen:]
            else:
                stateOutBytes = ''

            #print "SOBi:", stateOutBytes
            stateOutBytes = map(ord, stateOutBytes.decode('hex'))
            for i in xrange(len(successorMatchMask)):
                stateOutBytes[i*2] &= successorMatchMask[i]
                stateOutBytes[i*2+1] &= successorMatchMask[i]
            stateOutBytes = ''.join(["%02x" % x for x in stateOutBytes])

            #print "SOBo:", stateOutBytes

            prevMatchBytes = context.matchBytes[chain.ingress()]
            if len(prevMatchBytes) == 0:
                print "YUP"
                prevMatchBytes = set('')
            print "cover: %s   lookupValue:%s   stateInBytes: %s   prevMatchBytes: %s" % (cover, lookupValueFull, stateInBytes, prevMatchBytes)
            stateInBytesAsArray = map(ord, stateInBytes.decode('hex'))
            stateInBytesMask = stateInBytesAsArray[0::2]
            stateInBytesMatch = stateInBytesAsArray[1::2]
            for pmb in prevMatchBytes:
                pmbAsArray = map(ord, pmb.decode('hex'))
                pmbMask = pmbAsArray[0::2]
                pmbMatch = pmbAsArray[1::2]

                # Check if the masks are compatible
                compat = True
                for i in xrange(len(pmbMask)):
                    compat &= (stateInBytesMask[i] & pmbMask[i]) ^ stateInBytesMask[i] == 0
                    compat &= (stateInBytesMask[i] & pmbMatch[i]) == (stateInBytesMask[i] & stateInBytesMatch[i])

                print "COMPAT:", compat, stateInBytes, pmb

                if not compat:
                    continue

                lookupValue = lookupValueFull[baseMatchByteLen:]
                print "cover: %s   lookupValue:%s   stateInBytes (pmb): %s   stateOutBytes: %s" % (cover, lookupValue, pmb, stateOutBytes)
                stateIn = baseMatchBytes.index(pmb) + baseState
                stateOut = successorMatchBytes.index(stateOutBytes) + successorState

                #print cover, cover.getLookupByteValues(False)
                mask = [0 for x in xrange(lookups * lookupWidth)]
                match = [0 for x in xrange(lookups * lookupWidth)]
                lookupValueByteArray = map(ord, lookupValue.decode('hex'))
                for i in xrange(len(lookupValueByteArray) / 2):
                    dest = lookupMap[i]
                    mask[dest] = lookupValueByteArray[i*2]
                    match[dest] = lookupValueByteArray[i*2+1]

                stateMask = []
                stateMatch = []
                for i in xrange(stateBytes):
                    maskVal = 255
                    if i == stateBytes - 1:
                        maskVal = 2 ** (stateBits % 8) - 1
                        if maskVal == 0:
                            maskVal = 255
                    matchVal = (stateIn >> (i * 8)) & 0xff
                    stateMask.insert(0, maskVal)
                    stateMatch.insert(0, matchVal)
                stateMask.extend(mask)
                stateMatch.extend(match)

                #print cover, "  State-in:", stateIn, "  State-out:", stateOut, "  Lookup-values:", lookupValue, "  Span:", skip, "  Mask:", stateMask, "  Match:", stateMatch

                wildcardMatch = False
                if len(lookupValue) == 0 and isSuccessorInSameHdr(context, cover, successor):
                    wildcardMatch = True
                    numBaseMatchBytes = len(baseMatchBytes)
                    baseMaskSize = int(math.ceil(math.log(numBaseMatchBytes, 2)))
                    nxtState = stateOut & (2 ** stateBits - 1 - (2 ** baseMaskSize - 1))
                    nxtStateMask = (2 ** stateBits - 1 - (2 ** baseMaskSize - 1))
                    pos = stateBytes - 1
                    while baseMaskSize > 0:
                        bits = baseMaskSize
                        if bits > 8:
                            bits = 8
                        stateMask[pos] &= ~(2 ** bits - 1)
                        stateMatch[pos] &= ~(2 ** bits - 1)
                        baseMaskSize -= 8
                else:
                    nxtState = stateOut
                    nxtStateMask = 2 ** stateBits - 1

                #print cover, "  State-in:", stateIn, "  State-out:", stateOut, "  Lookup-values:", lookupValue, "  Span:", skip, "  Mask:", stateMask, "  Match:", stateMatch

                mask = stateMask
                match = stateMatch

                entryStr = 'Match: ([%s]' % (', '.join(['%02x' % val for val in mask]))
                entryStr += ', [%s])' % (', '.join(['%02x' % val for val in match]))
                entryStr += '   Next-State: ' + stateFmtStr % (nxtState, nxtStateMask)
                entryStr += '   Adv: ' + advFmtStr % skip
                entryStr += '   Next-Lookup: [%s]' % (', '.join('%d' % val for val in nxtLookupStarts))
                if extract:
                    coverExtractBytes = []
                    base = 0
                    #print "EXT:", cover
                    for node in cover.chain:
                        #nodeExtractBytes = node.dagNode.hdr.getExtractBytes()
                        startPos = node.startPos
                        if isinstance(node.dagNode, SubHeaderNode):
                            startPos += node.dagNode.startPos
                        elif isinstance(node.dagNode, BarrierNode):
                            startPos += node.dagNode.barrierLoc
                        elif isinstance(node.dagNode, PadNode):
                            base += node.consumed
                            continue
                        endPos = startPos + node.consumed
                        nodeName = '%s-%d' % (node.dagNode.hdr.name, node.dagNode.inst)
                        nodeExtractBytes = copy.copy(context.fieldPos[nodeName])
                        while len(nodeExtractBytes) > 0 and nodeExtractBytes[0][0] < startPos:
                            nodeExtractBytes.pop(0)
                        while len(nodeExtractBytes) > 0 and nodeExtractBytes[-1][0] >= endPos:
                            nodeExtractBytes.pop()
                        nodeExtractBytes = [(src - startPos + base, dest) for (src, dest) in nodeExtractBytes]
                        coverExtractBytes.extend(nodeExtractBytes)
                        base += node.consumed
                        #print "EXT:", nodeName, nodeExtractBytes
                    #entryStr += '   Extract: [%s]' % (', '.join('(%d, %d)' % (src, dst) for (src, dst) in coverExtractBytes))
                    coverExtractBytes.extend([(0, 0) for x in xrange(extractBytes - len(coverExtractBytes))])
                    entryStr += '   Extract: [%s]' % (', '.join(extractFmtStr % (src, dst) for (src, dst) in coverExtractBytes))
                    #print "EXT:", coverExtractBytes

                    #extractBytes = cover.getExtractBytes()
                    #entryStr += '   Extract: [%s]' % (', '.join('%d' % val for val in extractBytes))
                starts = []
                pos = 0
                for node in cover.chain:
                    dagNode = node.dagNode
                    nodeStr = '%s-%d' % (dagNode.hdr.name, dagNode.inst)
                    if isinstance(dagNode, HeaderNode):
                        if node.startPos == 0:
                            starts.append((context.foundHdrNum[nodeStr], pos, len(context.fieldPos[nodeStr])))
                    elif isinstance(dagNode, SubHeaderNode):
                        if dagNode.startPos == 0 and node.startPos == 0:
                            starts.append((context.foundHdrNum[nodeStr], pos, len(context.fieldPos[nodeStr])))
                    pos += node.consumed
                #entryStr += '   Hdr-Starts: [%s]' % (', '.join('(%d, %d)' % (hdr, pos) for (hdr, pos) in starts))
                starts.extend([(0, 0, 0) for x in xrange(maxHdrs - len(starts))])
                entryStr += '   Hdr-Starts: [%s]' % (', '.join(hdrStartFmtStr % (pos, hdr, length) for (hdr, pos, length) in starts))

                #print "  %s   # Match: %s   Nxt-State: %d" % (entryStr, cover, stateOut)
                commentStr = "Match: %s   Nxt-State: %d" % (cover, stateOut)
                entries.append((mask, match, entryStr, commentStr))


                #print cover, "  State-in:", stateIn, "  State-out:", stateOut, "  Lookup-values:", lookupValue, "  Span:", skip, "  Mask:", mask, "  Match:", match

                if wildcardMatch:
                    break
            if wildcardMatch:
                break

    for (mask, match, entryStr, commentStr) in sorted(entries, cmpEntry):
        if printTCAM:
            print "  %s     # %s" % (entryStr, commentStr)
        if saveTCAM:
            tcamFile.write("  %s     # %s\n" % (entryStr, commentStr))


def firstInstCNode(cnode):
    """
    Return a DAGChainNode that points to a first instance
    """
    if cnode.dagNode.inst == 1:
        return cnode
    else:
        dagNode = cnode.dagNode
        if type(dagNode) == HeaderNode:
            newDAGNode = HeaderNode(dagNode.hdr, 1)
        elif type(dagNode) == SubHeaderNode:
            newDAGNode = SubHeaderNode(dagNode.hdr, 1, dagNode.startPos, dagNode.endPos)
        elif type(dagNode) == PadNode:
            newDAGNode = PadNode(dagNode.hdr, 1)
        elif type(dagNode) == BarrierNode:
            newDAGNode = BarrierNode(dagNode.hdr, 1, dagNode.barrierLoc)
        else:
            raise TypeError('Unknown type of DAG node: %s (%s)' % (type(dagNode), dagNode))
        newCNode = DAGChainNode(newDAGNode, cnode.startPos, cnode.consumed, cnode.read)
        return newCNode

def printTCAMEntries(context):
    global tcamFile

    if len(context.optNodes) > 0:
        sortDAG(context)

        calcNumStatesNeeded(context)


        firstNode = context.dagOrderList[0]
        firstCluster = None

        # Identify the state requirements of each cluster
        states = 0
        clustersByState = {}
        for cluster in context.bestClusters:
            cnode = cluster.chain[0]
            statesNeeded = numStatesNeeded(cluster)
            states += statesNeeded
            if statesNeeded not in clustersByState:
                clustersByState[statesNeeded] = []
            clustersByState[statesNeeded].append(cluster)

            if cnode.dagNode == firstNode and cnode.startPos == 0:
                firstCluster = cluster
        #for cluster in sorted(context.bestClusters):
        #    print cluster, isNextNodeInSameHdr(context, cluster)

        # Allocate state to each cluster
        availableState = range(tcamMaxState)
        context.allocation = {}
        context.matchBytes = {}
        context.matchMask = {}
        allocateState(context, firstCluster, availableState)
        for states in reversed(sorted(clustersByState)):
            for cluster in clustersByState[states]:
                if cluster == firstCluster:
                    continue
                if not allocateState(context, cluster, availableState):
                    print "Insufficient state available for table"
                    return
        #for cluster in sorted(context.bestClusters):
        #    #print cluster, context.allocation[cluster.ingress()]
        #    print cluster.ingress(), context.allocation[cluster.ingress()]

        
        dstFileName = hfile.rstrip('.txt')
        dstFileName += '-w%03d-ms%03d-l%02d-lw%02d-flaz%d-e%d' % (windowSize, maxSkip, lookups, lookupWidth, firstLookupAtZero, extract)
        if extract:
            dstFileName += '-eb%02d' % extractBytes
        else:
            dstFileName += '-mh%02d' % maxHdrs
        dstFileName += '.tcam'
        if saveTCAM:
            tcamFile = open(dstFileName, 'w')
        if printTCAM:
            print "tcam-entries:"

        firstNode = context.dagOrderList[0]
        for firstChain in context.bestClusters:
            firstCNode = firstChain.chain[0]
            if firstCNode.dagNode == firstNode and firstCNode.startPos == 0:
                firstLookupBytes = firstChain.getLookupBytes()
                break
        firstLookupStarts = []
        if firstLookupAtZero:
            firstLookupStarts.append(0)
            while len(firstLookupBytes) > 0 and firstLookupBytes[0] < lookupWidth:
                firstLookupBytes.pop(0)
        while len(firstLookupBytes) > 0:
            firstLookupStarts.append(firstLookupBytes.pop(0))
            while len(firstLookupBytes) > 0 and firstLookupBytes[0] < lookupWidth + firstLookupStarts[-1]:
                firstLookupBytes.pop(0)
        while len(firstLookupStarts) < lookups:
            firstLookupStarts.append(0)
        firstLookupStr = 'First-Lookup: [%s]' % (', '.join('%d' % val for val in firstLookupStarts))
        if printTCAM:
            print "  %s" % firstLookupStr
        if saveTCAM:
            tcamFile.write("  %s\n" % firstLookupStr)

        for states in sorted(clustersByState):
            for cluster in sorted(clustersByState[states]):
                printTCAMEntry(context, cluster)
        if saveTCAM:
            tcamFile.close()


def exploreGraph(context, cnode):
    """
    Attempt to find the optimal TCAM allocation for a DAG rooted at a given
    chain node and starting at position pos. (Assume decision bytes prior to pos have
    been processed.)

    Return:
      [(lookupClusters, worstByteCount, worstCyc)]
    """
    #print "exploreGraph: N=%s" % (cnode)

    # Check if we've already done the calculation for the requested node
    if cnode in context.exploreResults:
        return context.exploreResults[cnode]

    #print "exploreGraph: N=%s" % (cnode)

    worstByteCount = max(maxSkip, windowSize) + 1
    worstCyc = 1
    worstBPC = (worstByteCount * 1.0) / worstCyc

    results = []

    startCN = cnode
    while isinstance(startCN.dagNode, BarrierNode) and len(startCN.dagNode.nxt) == 1 and list(startCN.dagNode.nxt)[0] is not None:
        startCN = DAGChainNode(list(startCN.dagNode.nxt)[0], 0, 0, 0)

    clusters = findClustersAndCovers(context, startCN)
    clusterCount = []
    #print "ICC", clusterCount
    #for cluster in clusters:
    #    shortestLen = findShortestLen(cluster, clusters)
    #for cluster in sorted(clusters):
    #    print "Cluster: %s" % (cluster)
    #    print "  Succ:",
    #    for succ in sorted(cluster.successors()):
    #        print succ,
    #    print ""
    #    fringe = findFringe(context, cluster, clusters)
    #    print "  Fringe:",
    #    for f in sorted(fringe):
    #        print f,
    #    print ""
    #print ""

    combs = 0
    counts = 0
    cn = 0
    for cluster in sorted(clusters):
        cn += 1
        fringeMult = 1


        #print "%s: Cluster: %s" % (cnode, cluster)
        #print "Cluster: %s" % (cluster)
        #continue

        edges = findEdgeCount(context, cluster, clusters)
        #print "  Edges: %d" % edges
        shortestLen = findShortestLen(context, cluster, clusters)

        #print "Shortest:", shortestLen

        fringe = findFringe(context, cluster, clusters)
        #print "%s: Fringe:" % cluster,
        #for f in fringe:
        #    print "    %s" % (f),
        #print ""
        #print "numFringe:", len(fringe)
        for fringeNode in fringe:
            #print fringeNode
            #print "Fringe:", fringeNode
            #print "%s: Fringe: %s" % (cnode, fringeNode)
            newResultSets = []

            fringeOpt = exploreGraph(context, fringeNode)
            (fringeClusterCount, fringeArrangements) = context.exploreResultCount[fringeNode]
            fringeMult *= fringeClusterCount
        #print ""
        #print "%s: Cluster: %s" % (cnode, cluster)
        #print "  FringeMult:", fringeMult
        childTotals = 0
        for fringeNode in fringe:
            #print "Fringe:", fringeNode
            newResultSets = []

            fringeOpt = exploreGraph(context, fringeNode)
            (fringeClusterCount, fringeArrangements) = context.exploreResultCount[fringeNode]
            #print "    chilldTotals inc:", fringeArrangements * fringeMult / fringeClusterCount
            childTotals +=  fringeArrangements * fringeMult / fringeClusterCount

        #print "  childTotals:", childTotals
        combs += childTotals + fringeMult
        #print "  combsInc:", childTotals + fringeMult
        counts += fringeMult
        #print "  FringeMult:", fringeMult, counts

        if len(fringe) > 0:
            newByteCount = max(maxSkip, windowSize) + 1
            newCyc = 1
        else:
            newByteCount = 0
            newCyc = 0
        resultSets = [(set(), newByteCount, newCyc)]
        #if len(fringe) > 0:
        #    clusterCount.append(0)
        #else:
        #    clusterCount.append(1)
        for fringeNode in fringe:
            #print "Fringe:", fringeNode
            newResultSets = []

            fringeOpt = exploreGraph(context, fringeNode)
            #fringeCounts = context.exploreResultCount[fringeNode]
            #for fc in fringeCounts:
            #    clusterCount[-1] += fc + 1
            #print fringeNode, clusterCount
            #print "FO:", fringeOpt
            #print "RS-Size:", len(resultSets)
            for (currClusters, currByteCount, currCyc) in resultSets:
                for (fClusters, fByteCount, fCyc) in fringeOpt:
                    newClusters = copy.copy(currClusters)
                    newByteCount = currByteCount
                    newCyc = currCyc

                    newClusters.update(fClusters)
                    newBPC = newByteCount * 1.0 / newCyc
                    fBPC = fByteCount * 1.0 / fCyc
                    if fBPC < newBPC:
                        newByteCount = fByteCount
                        newCyc = fCyc

                    newResultSets.append((newClusters, newByteCount, newCyc))
            #resultSets = newResultSets

        lenSum = 0
        for (newClusters, newByteCount, newCyc) in resultSets:
            newClusters.add(cluster)
            newByteCount += shortestLen
            newCyc += 1
            lenSum += len(newClusters)
            #if resultsRet == 7 and cn == 1:
            #    print "Clusters:"
            #    for c in sorted(newClusters):
            #        print "   ", c


        #print "  [",
        #for cluster in sorted(newClusters):
        #    print cluster,
        #print "] %d %d %1.5f" % (newByteCount, newCyc, newByteCount * 1.0 / newCyc)

            #print "  Fringe clusters:", node.dagNode, node.startPos

            results.append((newClusters, newByteCount, newCyc))
        #print "RS-Size:", len(resultSets), lenSum

    if debug:
        print "exploreGraph: N=%s" % (cnode)
        print "  Results:"
        for (clusters, byteCount, cyc) in results:
            #print "  Result:"
            for cluster in sorted(clusters):
                print "    %s" % cluster
        print ""
        resultCount = 0
        for (clusters, byteCount, cyc) in results:
            resultCount += len(clusters)
        print "  ResultCount:", resultCount
        print "  Counts: %d   Combs: %d" % (counts, combs)
        print ""

    #global resultsRet
    #resultsRet += 1
    #if resultsRet == resultsTerminate:
    #    sys.exit(1)

    context.exploreResults[cnode] = results
    #context.exploreResults[cnode] = []
    #context.exploreResultCount[cnode] = clusterCount
    #context.exploreResultCount[cnode] = (len(clusters), combs)
    context.exploreResultCount[cnode] = (counts, combs)
    #print "exploreGraph: N=%s, C=%d" % (cnode, clusterCount)


    return results

def extendChain(chain,
        lookups=None, lookupWidth=None,
        maxSkip=None, firstLookupAtZero=None,
        windowSize=None, extract=None, extractBytes=None,
        maxHdrs=None):
    """
    Extend a chain staying within the boundaries of the various parameters
    """
    # Populate defaults
    if lookups is None:
        lookups = globals()['lookups']
    if lookupWidth is None:
        lookupWidth = globals()['lookupWidth']
    if maxSkip is None:
        maxSkip = globals()['maxSkip']
    if firstLookupAtZero is None:
        firstLookupAtZero = globals()['firstLookupAtZero']
    if windowSize is None:
        windowSize = globals()['windowSize']
    if extract is None:
        extract = globals()['extract']
    if extractBytes is None:
        extractBytes = globals()['extractBytes']
    if maxHdrs is None:
        maxHdrs = globals()['maxHdrs']

    #print "extendChain: C=%s, LU=%d, LUW=%d, MS=%d, FLUAZ=%s, WS=%d, E=%s, EB=%d MH=%d" % (
    #        chain, lookups, lookupWidth, maxSkip,
    #        firstLookupAtZero, windowSize, extract, extractBytes,
    #        maxHdrs)

    if chain.done or len(chain.chain) == 0:
        return False

    origChainConsumed = chain.totalConsumed()
    origChainRead = chain.totalRead()

    # Temporarily extend the last node
    lastCN = chain.chain[-1]
    origLastConsumed = lastCN.consumed
    origLastRead = lastCN.read
    lastCN.read += lastCN.unread()
    lastCN.consumed += lastCN.unconsumed()

    # Identify the first lookup byte that can't be processed
    # (due to insufficient lookups or sitting outside the window)
    allLB = chain.getLookupBytes()
    #print "Chain: %s   AllLB: %s" % (chain, allLB)
    if firstLookupAtZero:
        lookups -= 1
        while len(allLB) > 0 and allLB[0] < lookupWidth and allLB[0] < windowSize:
            allLB.pop(0)
    while len(allLB) > 0 and lookups > 0 and allLB[0] < windowSize:
        lookupStart = allLB.pop(0)
        lookups -= 1
        while len(allLB) > 0 and allLB[0] < lookupStart + lookupWidth and allLB[0] < windowSize:
            allLB.pop(0)
    # Drop any lookup bytes that may already be consumed in case we've got all
    # zeros for lookups etc
    while len(allLB) > 0 and allLB[0] < origChainRead:
        allLB.pop(0)
    #print "Chain: %s   AllLB: %s" % (chain, allLB)

    # Identify the first extract byte that can't be processed
    # (due to insufficient extracts or sitting outside the window)
    allEB = []
    if extract:
        allEB = chain.getExtractBytes()
        #print "Chain: %s   AllEB: %s" % (chain, allEB)
        while len(allEB) > 0 and extractBytes > 0 and allEB[0] < windowSize:
            allEB.pop(0)
            extractBytes -= 1
        # Drop any lookup bytes that may already be consumed in case we've
        # got all zeros for lookups etc
        while len(allEB) > 0 and allEB[0] < origChainConsumed:
            allEB.pop(0)

    # Extend the read/consumed counters
    maxSize = max(maxSkip, windowSize)
    allLB.append(maxSize)
    allEB.append(maxSize)
    newConsumed = min(chain.totalConsumed(), maxSize, allLB[0], allEB[0])
    newRead = min(chain.totalRead(), maxSize, allLB[0], allEB[0])
    lastCN.consumed = origLastConsumed + newConsumed - origChainConsumed
    lastCN.read = origLastRead + newRead - origChainRead

    #print "Exp:", chain, newConsumed, newRead, allLB, allEB

    chain.calcPatterns()

    return lastCN.unconsumed() == 0 and lastCN.unread() == 0

def extendChainNode(cnode,
        unusedLookupBytes, lookups, lookupWidth, maxSkip, firstLookupAtZero,
        windowSize, extract, extractBytes, maxHdrs):
    """
    Extend a chain node staying within the boundaries of the various parameters
    """
    #print "extendChainNode: N=%s, ULB=%d, LU=%d, LUW=%d, MS=%d, FLUAZ=%s, WS=%d, E=%s, EB=%d MH=%d" % (
    #        cnode, unusedLookupBytes, lookups, lookupWidth, maxSkip,
    #        firstLookupAtZero, windowSize, extract, extractBytes,
    #        maxHdrs)

    # Identify the first lookup byte that can't be processed
    # (due to insufficient lookups or sitting outside the window)
    unprocLB = cnode.unprocLookupBytes()
    if unusedLookupBytes:
        while len(unprocLB) > 0 and unprocLB[0] < unusedLookupBytes and unprocLB[0] < windowSize:
            unprocLB.pop(0)
    while len(unprocLB) > 0 and lookups > 0 and unprocLB[0] < windowSize:
        lookupStart = unprocLB.pop(0)
        lookups -= 1
        while len(unprocLB) > 0 and unprocLB[0] < lookupStart + lookupWidth and unprocLB[0] < windowSize:
            unprocLB.pop(0)

    # Identify the first extract byte that can't be processed
    # (due to insufficient extracts or sitting outside the window)
    unprocEB = []
    if extract:
        unprocEB = cnode.unprocExtractBytes()
        while len(unprocEB) > 0 and extractBytes > 0 and unprocEB[0] < windowSize:
            unprocEB.pop(0)
            extractBytes -= 1

    # Extend the read/consumed counters
    maxExtend = max(maxSkip, windowSize)
    unprocLB.append(maxExtend)
    unprocEB.append(maxExtend)
    unconsumedBytes = cnode.unconsumed()
    unreadBytes = cnode.unread()
    extendConsumed = min(unconsumedBytes, maxExtend, unprocLB[0], unprocEB[0])
    extendRead = min(unreadBytes, maxExtend, unprocLB[0], unprocEB[0])
    cnode.consumed += extendConsumed
    cnode.read += extendRead

    #print "Exp:", cnode

    return extendConsumed == unconsumedBytes and extendRead == unreadBytes

lookupDests = {}


def findClustersAndCovers(context, cnode):
    """
    Identify all possible clusters from cnode reachable
    within a single cycle.
    """
    #clusterStr = "%s--%d" % (node.shortStr(), pos)
    if cnode in context.coveredClustersFromNode:
        return context.coveredClustersFromNode[cnode]

    #print "findClustersAndCover: N=%s" % (cnode)

    baseChain = DAGChain()
    baseChain.add(copy.copy(cnode))
    
    chains = findMaxLenChains(baseChain, lookups, lookupWidth, maxSkip,
            firstLookupAtZero, windowSize, extract, extractBytes, maxHdrs)
    #print "findClustersAndCover: N=%s" % (cnode)
    #for chain in sorted(chains):
    #    print "CHAIN:", chain
    #sys.exit(1)

    coveredClusters = findCoveredClusters(chains)
    #print "findClustersAndCover: N=%s" % (cnode)
    #for chain in sorted(coveredClusters):
    #    print "CHAIN (post):", chain
    #    for cover in sorted(coveredClusters[chain]):
    #        print "  ", cover
    #print ""
    context.coveredClustersFromNode[cnode] = coveredClusters
    context.clustersFromNode[cnode] = sorted(coveredClusters.keys())

    return coveredClusters

def findMaxLenChains(baseChain,
        lookups=None, lookupWidth=None,
        maxSkip=None, firstLookupAtZero=None,
        windowSize=None, extract=None, extractBytes=None,
        maxHdrs=None):
    """
    Identify all maximum length chains reachable from a node starting at
    position pos reachable within a single cycle.
    """
    if lookups is None:
        lookups = globals()['lookups']
    if lookupWidth is None:
        lookupWidth = globals()['lookupWidth']
    if maxSkip is None:
        maxSkip = globals()['maxSkip']
    if firstLookupAtZero is None:
        firstLookupAtZero = globals()['firstLookupAtZero']
    if windowSize is None:
        windowSize = globals()['windowSize']
    if extract is None:
        extract = globals()['extract']
    if extractBytes is None:
        extractBytes = globals()['extractBytes']
    if maxHdrs is None:
        maxHdrs = globals()['maxHdrs']

    #print "findMaxLenChains: BC=%s, LU=%d, LUW=%d, MS=%d, FLUAZ=%s, WS=%d, E=%s, EB=%d MH=%d" % (
    #        baseChain, lookups, lookupWidth, maxSkip,
    #        firstLookupAtZero, windowSize, extract, extractBytes, maxHdrs)

    chains = [baseChain]
    chainsDone = []
    while len(chains) > 0:
        # Attempt to expand each chain
        newChains = []
        for chain in chains:
            #print "Base:", chain
            lastCN = chain.lastNode()
            if lastCN:
                if len(chain.chain) > 1 and isinstance(lastCN.dagNode, BarrierNode):
                    chainsDone.append(chain)
                    continue

                # Identify the resources we have left
                lookupsUsed, unusedLookupBytes = chain.lookupsUsed(lookupWidth, firstLookupAtZero)

                lookupsRem = lookups - lookupsUsed
                skipRem = maxSkip - chain.totalConsumed()
                windowRem = windowSize - chain.totalConsumed()

                extractRem = extractBytes
                hdrsRem = maxHdrs
                if extract:
                    extractRem -= chain.totalExtracted()
                else:
                    hdrsRem -= chain.hdrsFound()

                # Adjust values if we've exhausted the window
                if windowRem <= 0:
                    windowRem = 0
                    lookupsRem = 0
                    unusedLookupBytes = 0
                    extractRem = 0

                #print "  Last:", lastCN, lookupsRem, skipRem, windowRem, extractRem, hdrsRem

                # Attempt to extend the header chain
                if lastCN.unconsumed() > 0 or lastCN.unread() > 0:
                    extended = extendChain(chain,
                            lookups, lookupWidth, maxSkip, firstLookupAtZero,
                            windowSize, extract, extractBytes, maxHdrs)
                    #extended = extendChainNode(lastCN,
                    #        unusedLookupBytes, lookupsRem, lookupWidth,
                    #        skipRem, firstLookupAtZero, windowRem, extract,
                    #        extractRem, hdrsRem)
                    if extended:
                        newChains.append(chain)
                    else:
                        chainsDone.append(chain)
                else:
                    for nxtNode in lastCN.dagNode.nxt:
                        if nxtNode:
                            newChain = chain.copy()
                            newChain.add(nxtNode, 0, 0, 0)
                            newCN = newChain.lastNode()
                            if not isinstance(nxtNode, BarrierNode):
                                extended = extendChain(newChain,
                                        lookups, lookupWidth, maxSkip, firstLookupAtZero,
                                        windowSize, extract, extractBytes, maxHdrs)
                                #extended = extendChainNode(newCN,
                                #        unusedLookupBytes, lookupsRem, lookupWidth,
                                #        skipRem, firstLookupAtZero, windowRem, extract,
                                #        extractRem, hdrsRem)
                                if extended:
                                    newChains.append(newChain)
                                else:
                                    #newChain.trim()
                                    chainsDone.append(newChain)
                            else:
                                chainsDone.append(newChain)
                        else:
                            newChain = chain.copy()
                            newChain.add(None)
                            chainsDone.append(newChain)
            else:
                chainsDone.append(chain)
        chains = newChains
        #print ""
    return chainsDone

def findCoveredClusters(clusters):
    """Find all nodes reachable from each cluster"""

    # Identify the lookup chains
    covers = dict()
    (reachableNewHdr, reachableSameHdr) = findReachable(clusters)
    #for chain in sorted(clusters):
    #    print "CLU:", chain
    #for subchain in sorted(reachableNewHdr):
    #    print "RNR:", subchain
    #    for r in sorted(reachableNewHdr[subchain]):
    #        print "      ", r
    #sys.exit(1)
    for cluster in reachableNewHdr.keys():
        #print "CLU:", cluster
        subchains = cluster.findLookupSubchains()

        # Calculate the set of nodes reachable from the current cluster
        targets = set()
        for subchain in sorted(subchains):
            #print "Add:", subchain
            targets.update(reachableNewHdr[subchain])
        targets.update(reachableSameHdr[subchains[-1]])
        #for t in sorted(targets):
        #    print "T:", t

        # Trim the list to ensure only the longest chains ending at each node
        # are included.
        # Be careful about cases where multiple paths lead to the SAME final node!
        #toRemove = []
        toRemove = set()
        longest = {}
        for target in targets:
            #print "TAR:", target
            targetOrig = target
            lastCN = target.lastNode(True)
            lastDagNode = lastCN.dagNode
            origPos = lastCN.startPos + lastCN.read
            #origPos = target.totalSpan()
            while target.nodeCount() > 0:
                lastCN = target.lastNode(True)
                #lastCN = target.lastNode()
                if target.done:
                    lastDagNode = None
                else:
                    lastDagNode = lastCN.dagNode
                lastPos = lastCN.startPos + lastCN.read
                targetMinusOne = target.copy()
                if target.done:
                    targetMinusOne.done = False
                else:
                    targetMinusOne.chainPop()
                targetMinusOne.patterns=0
                pathStr = str(targetMinusOne) + '  --  ' + str(lastDagNode)
                #print "PATHSTR:", lastPos, origPos, target, pathStr

                if pathStr in longest:
                    #print "Seen it"
                    (otherPos, otherChain, otherOrig) = longest[pathStr]
                    #print otherPos, lastPos, origPos, otherOrig
                    #if otherPos < origPos:
                    if otherPos < lastPos:
                        toRemove.add(otherChain)
                        #print "TO-REM:", otherChain
                        longest[pathStr] = (lastPos, targetOrig, target == targetOrig)
                        #longest[pathStr] = (lastPos, target, target == targetOrig)
                        #if target != targetOrig:
                        #    toRemove.add(target)
                        #    print "TO-REM:", target
                    else:
                        if target == targetOrig:
                            toRemove.add(targetOrig)
                            #print "TO-REM:", targetOrig
                        #toRemove.add(target)
                        pass
                else:
                    longest[pathStr] = (lastPos, targetOrig, target == targetOrig)
                    #longest[pathStr] = (lastPos, target, target == targetOrig)
                    #if target != targetOrig:
                    #    toRemove.add(target)
                    #    print "TO-REM:", target

                target = targetMinusOne
            #print ""

        #for target in sorted(targets):
        #    print "TAR:", target
        for removeChain in toRemove:
            #print "REM:", removeChain
            targets.remove(removeChain)

        covers[cluster] = targets
        #sys.exit(1)

    return covers


def findReachable(chains,
        lookups=None, lookupWidth=None,
        maxSkip=None, firstLookupAtZero=None,
        windowSize=None, extract=None, extractBytes=None,
        maxHdrs=None):
    """
    Identify all headers reachable from a given chain.

    Reachable headers fork off after each decision point (for chains with
    decision points) or at header boundaries (when there is no decision point)
    """

    if lookups is None:
        lookups = globals()['lookups']
    if lookupWidth is None:
        lookupWidth = globals()['lookupWidth']
    if maxSkip is None:
        maxSkip = globals()['maxSkip']
    if firstLookupAtZero is None:
        firstLookupAtZero = globals()['firstLookupAtZero']
    if windowSize is None:
        windowSize = globals()['windowSize']
    if extract is None:
        extract = globals()['extract']
    if extractBytes is None:
        extractBytes = globals()['extractBytes']
    if maxHdrs is None:
        maxHdrs = globals()['maxHdrs']

    reachableNewHdr = dict()
    reachableSameHdr = dict()
    #for chain in chains:
    for chain in sorted(chains):
        #print "CH:", chain
        lookupSubchains = chain.findLookupSubchains()
        #for lc in lookupSubchains:
        for lc in sorted(lookupSubchains):
            #print "   LC:", lc
            #continue
            if lc not in reachableNewHdr:

                chainBase = lc.copy()
                extendedChains = findMaxLenChains(chainBase,
                        0, 0, maxSkip, False, windowSize, extract, extractBytes, maxHdrs)

                extendedChainsNewHdr = set()
                extendedChainsSameHdr = set()
                for ec in extendedChains:
                    if ec.chain[-1].dagNode != lc.chain[-1].dagNode or \
                            ec.done != lc.done:
                        extendedChainsNewHdr.add(ec)
                    else:
                        extendedChainsSameHdr.add(ec)

                reachableNewHdr[lc] = extendedChainsNewHdr
                reachableSameHdr[lc] = extendedChainsSameHdr

    return reachableNewHdr, reachableSameHdr

def findEdgeCount(context, cluster, clusters=None):
    """Count the number of edges from a cluster"""
    if clusters:
        if cluster not in context.edgesPerCluster:
            edgeCount = 0
            for c in clusters[cluster]:
                #print "findEdgeCount: %s: %s: %d" % (cluster, c, c.patterns)
                edgeCount += c.patterns
            context.edgesPerCluster[cluster] = edgeCount
            #context.edgesPerCluster[cluster] = len(clusters[cluster])

    return context.edgesPerCluster[cluster]

def findFringe(context, cluster, clusters=None):
    """Find the 'fringe' nodes from a cluster"""
    if clusters:
        if cluster not in context.fringeByCluster:
            #print "findFringe:", cluster
            #for chain in sorted(clusters[cluster]):
            #    print "   ", chain,
            #    for c in chain.successors():
            #        print c,
            #    print ""
            fringe = set()
            for chain in clusters[cluster]:
                #print "  ", chain
                fringe.update(chain.successors())

            context.fringeByCluster[cluster] = sorted(fringe)
            #print "Fringe is:",
            #for f in fringe:
            #    print f,
            #print ""

    return context.fringeByCluster[cluster]

def findFringeWithMinCoverLen(context, cluster, coveredClusters=None):
    """
    Find the 'fringe' nodes from a cluster and the minimum cover width to
    reach the fringe
    """
    #print "findFringeWithMinCoverLen:", cluster
    if cluster not in context.fringeWithMinCoverLen:
        fringe = {}
        for cover in coveredClusters[cluster]:
            if not cover.done:
                coverLen = cover.totalConsumed()
                for s in cover.successors():
                    if s in fringe:
                        if coverLen < fringe[s]:
                            fringe[s] = coverLen
                    else:
                        fringe[s] = coverLen
        context.fringeWithMinCoverLen[cluster] = fringe

    return context.fringeWithMinCoverLen[cluster]


def findShortestLen(context, cluster, clusters):
    """Find the shortest chain length"""
    if cluster not in context.shortestLenByCluster:
        chains = clusters[cluster]

        if len(chains) == 0:
            return 0

        chainLen = max(maxSkip, windowSize) + 1
        for chain in chains:
            newChainLen = chain.totalConsumed()
            if newChainLen < chainLen:
                chainLen = newChainLen

        context.shortestLenByCluster[cluster] = chainLen

    return context.shortestLenByCluster[cluster]

def buildDAG(headerList, headers):
    seenPath = set()
    #finalPaths = set()

    dagNodes = {}

    def exploreHdrChain(hdr, nxtHdr, headers, hdrInfo, chain = None):
        prevInst = 0
        path = ""
        if chain:
            for hdrNode in chain:
                if hdrNode.hdr == hdr:
                    prevInst += 1
                path += "%s ->" % (hdrNode.hdr.name)
        else:
            chain = []
        refCount = hdr.getRefCount()
        if refCount is not None:
            prevInst = refCount.refCount - 1

        dagNodeStr = "%s:%d" % (hdr.name, prevInst + 1)
        path += "%s" % (hdr.name)
        
        if dagNodeStr in dagNodes:
            dagNode = dagNodes[dagNodeStr]
        else:
            dagNode = HeaderNode(hdr, prevInst + 1)
            dagNodes[dagNodeStr] = dagNode
        if len(chain) > 0:
            chain[-1].nxt.add(dagNode)
        else:
            buildDAG.headNode = dagNode
        chain = copy.copy(chain)
        chain.append(dagNode)

        if nxtHdr:
            nxtHdrName = nxtHdr.name
        else:
            nxtHdrName = "--"
        #print path
        pathNxt = path + " -> " + nxtHdrName
        if pathNxt not in seenPath:
            seenPath.add(pathNxt)
            if nxtHdr:
                ehcFunc = lambda ehcHdr, ehcNxtHdr, ehcHeaders, ehcInfo : \
                        exploreHdrChain(ehcHdr, ehcNxtHdr, ehcHeaders, ehcInfo, chain)
                exploreHeader(nxtHdr, headers, ehcFunc)
            else:
                dagNode.nxt.add(None)
                #finalPaths.add(path)
    
    
    #headHdr = Header(HEAD_NODE)
    #headNode = HeaderNode(headHdr, 0)
    buildDAG.headNode = None
    exploreHeader(headerList[0], headers, exploreHdrChain)

    addPad(dagNodes)

    return buildDAG.headNode, dagNodes

def addPad(dagNodes):
    """Add pad nodes to a DAG"""
    newNodes = {}
    for key in dagNodes.keys():
        node = dagNodes[key]
        lengths = node.getLengths()
        if len(lengths) > 1: # or lengths[0] > maxSkip:
            #minLen = min(lengths)
            #oldNxt = copy.copy(node.nxt)
            #for l in lengths:
            #    if l - minLen > 0:
            #        padNode = PadNode(node.hdr, node.inst, l - minLen)
            #        padNode.nxt = copy.copy(oldNxt)
            #        padNodeStr = '%s:%d:%02d' % (padNode.getName(), padNode.inst, padNode.length)
            #        newNodes[padNodeStr] = padNode
            #        node.nxt.add(padNode)
            minLen = min(lengths)
            oldNxt = node.nxt
            node.nxt = set()
            for l in lengths:
                padNode = PadNode(node.hdr, node.inst, l - minLen)
                padNode.nxt = copy.copy(oldNxt)
                padNodeStr = '%s:%d:%02d' % (padNode.getName(), padNode.inst, padNode.length)
                newNodes[padNodeStr] = padNode
                node.nxt.add(padNode)
    dagNodes.update(newNodes)


walkHistory = []

bestEdgeCount = -1
bestClusters = set
bestWorstBytes = max(windowSize, maxSkip) + 1
bestWorstCyc = 1
#bestWorstBPC = bestWorstBytes * 1.0 / bestWorstCyc
bestWorstBPC = 0
exploreCount = 0
exploreReportAmt = 100000
estReportTime = 10
countEstimated = False

def printParams():
    print "TCAM entry evaluation"
    print "====================="
    print ""
    print "Header file: %s" % hfile
    print "Data rate: %1.3f Gb/s" % dataRate
    print "Clock frequency: %1.3f GHz" % clkFreq
    print "Required processing rate: %1.3f bpc" % globalBPC
    print "Lookups: %d" % lookups
    print "Lookup width: %d bytes" % lookupWidth
    print "Maximum skip: %d bytes" % maxSkip
    print "Minimum skip: %d bytes" % minSkip
    print "First lookup at zero: %s" % firstLookupAtZero
    print "Window size: %d bytes" % windowSize
    print "Perform extraction: %s " % extract
    if extract:
        print "Extract: %d bytes" % extractBytes
    else:
        print "Maximum headers/cycle: %d" % maxHdrs
    print ""
    print "Evaluations:"
    print "  Multiple-parent merge: %s" % multParent
    print "  Multiple-parent merge second pass: %s" % multParentRetry
    print "  Parallel edge barriers: %s" % parallelEdge
    print "  Combine clusters different by instance number: %s" % instMerge
    print "  Ternary state matching: %s" % ternMatchOnState
    print ""
    print "Output:"
    print "  Print TCAM entries: %s" % printTCAM
    print "  Save TCAM entries: %s" % saveTCAM
    print "  TCAM maximum state: %d" % tcamMaxState

def sortDAG(context):
    # Identify the parents
    parents = {}
    parents[context.dag] = []
    context.dagParents = {}
    context.dagParents[context.dag] = []
    context.dagOrder = {}
    context.dagOrderList = []
    pending = set([context.dag])
    seen = set()
    while len(pending) > 0:
        curr = pending
        pending = set()
        for node in curr:
            for nxtNode in node.nxt:
                if nxtNode is not None:
                    if nxtNode not in context.dagParents:
                        parents[nxtNode] = []
                        context.dagParents[nxtNode] = []
                    parents[nxtNode].append(node)
                    context.dagParents[nxtNode].append(node)
                    if nxtNode not in seen:
                        pending.add(nxtNode)
                        seen.add(nxtNode)

    # Assign numbers
    depth = 0
    offset = 0
    pending = set([context.dag])
    while len(pending) > 0:
        curr = pending
        pending = set()
        for node in sorted(curr):
            context.dagOrder[node] = depth
            context.dagOrderList.append(node)
            node.setTopo(depth, offset)
            depth += 1
            offset += node.getLength()
            for nxtNode in node.nxt:
                if nxtNode is not None:
                    parentList = parents[nxtNode]
                    parentList.remove(node)
                    if len(parentList) == 0:
                        pending.add(nxtNode)

def cloneDAG(head):
    """Clone a DAG, duplicating all DAG nodes"""
    newHead = None

    clonedNodes = {}

    # Clone the nodes
    pending = [head]
    seen = set(pending)
    while len(pending) > 0:
        node = pending.pop(0)
        #print "Cloning %s..." % node.shortStr()
        nodeCopy = copy.copy(node)
        newNxt = []
        for nxt in nodeCopy.nxt:
            #print "  Nxt: %s" % nxt
            if nxt is None:
                newNxt.append(None)
            else:
                newNxt.append(nxt.shortStr())
                if nxt not in seen:
                    pending.append(nxt)
                    seen.add(nxt)
        nodeCopy.nxt = newNxt
        clonedNodes[nodeCopy.shortStr()] = nodeCopy
        if newHead is None:
            newHead = nodeCopy

    # Fix the next pointers
    pending = [newHead]
    seen = set()
    while len(pending) > 0:
        node = pending.pop(0)
        #print "Fixing nxt for %s..." % node.shortStr()
        newNxt = set()
        for nxt in node.nxt:
            #print "  Nxt: %s" % nxt
            if nxt is None:
                newNxt.add(None)
            else:
                nxtNode = clonedNodes[nxt]
                newNxt.add(nxtNode)
                if nxtNode not in seen:
                    pending.append(nxtNode)
                    seen.add(nxtNode)
        node.nxt = newNxt

    return newHead

def findDAGNodes(dag, nodes):
    # Find the requested nodes
    nodeCount = len(nodes)
    isFound = [False for node in xrange(nodeCount)]
    foundNodes = [None for node in xrange(nodeCount)]
    toFind = nodeCount

    pending = [dag]
    seen = set(pending)
    while len(pending) > 0 and toFind:
        node = pending.pop(0)
        for i in xrange(nodeCount):
            if not isFound[i]:
                if nodes[i] == node:
                    foundNodes[i] = node
                    isFound[i] = True
                    toFind -= 1
        for nxt in node.nxt:
            if nxt and nxt not in seen:
                pending.append(nxt)
                seen.add(nxt)

    if toFind == 0:
        return foundNodes
    else:
        return None

def insertDAGBarrier(head, targetNode, pos, nxtNode=None):
    """
    Insert a barrier in a DAG
    
    Place the barrier within (or at the end of) targetNode at position pos.

    If pos extends beyond the boundary of the node, add it at the end of the node.
    Will search out pads as well.

    A nxtNode can be specified when the pos is at the end of a node.
    """
    # Handle pos being beyond the end of the node
    if pos > targetNode.getLength():
        nxtIsPad = len(targetNode.nxt) > 0
        for nxt in targetNode.nxt:
            nxtIsPad = nxtIsPad and isinstance(nxt, PadNode)
        if nxtIsPad:
            for nxt in targetNode.nxt:
                insertDAGBarrier(head, nxt, pos - targetNode.getLength(), nxtNode)
            return
        else:
            pos = targetNode.getLength()

    # Find the targetNode and the nxtNode in the graph
    nodes = [targetNode]
    if nxtNode:
        nodes.append(nxtNode)
    nodes = findDAGNodes(head, nodes)
    if nodes is None:
        return

    targetNode = nodes[0]
    if nxtNode:
        nxtNodes = nodes[1]

    # Handle when the pos is the end of the target node
    if pos == targetNode.getLength() and pos > 0:

        if nxtNode:
            barrierNode = BarrierNode(nxtNode.hdr, nxtNode.inst, 0)
            nodes = findDAGNodes(head, [barrierNode])
            if nodes is not None:
                barrierNode = nodes[0]
            else:
                barrierNode.nxt.add(nxtNode)
            for nxt in targetNode.nxt:
                if nxt == nxtNode:
                    targetNode.nxt.add(barrierNode)
                    targetNode.nxt.remove(nxt)
                    break
        else:
            barrierNode = BarrierNode(targetNode.hdr, targetNode.inst, pos)

            # Attempt to ensure that we aren't inserting a duplicate barrier
            nxtIsBarrier = False
            barrierCnt = 0
            for nxt in targetNode.nxt:
                if nxt and isinstance(nxt, BarrierNode) and nxt == barrierNode:
                    nxtIsBarrier = True
                    barrierCnt += 1
            if nxtIsBarrier and barrierCnt == len(targetNode.nxt):
                return

            barrierNode = BarrierNode(targetNode.hdr, targetNode.inst, pos)
            barrierNode.nxt = targetNode.nxt
            targetNode.nxt = set([barrierNode])

        return

    # If we make it to here then pos is within the header
    barrierNode = BarrierNode(targetNode.hdr, targetNode.inst, pos)

    if pos == 0:
        headNode = barrierNode
        tailNode = targetNode
    else:
        headNode = SubHeaderNode(targetNode.hdr, targetNode.inst, 0, pos-1)
        headNode.nxt.add(barrierNode)
        tailNode = SubHeaderNode(targetNode.hdr, targetNode.inst, pos, targetNode.getLength()-1)
        tailNode.nxt = copy.copy(targetNode.nxt)
    barrierNode.nxt.add(tailNode)

    #tailNode = barrierNode
    #tailNode.nxt.add(None)

    # Perform the replacement
    pending = [head]
    seen = set(pending)
    while len(pending) > 0:
        node = pending.pop(0)
        newNxt = set()
        for nxt in node.nxt:
            if nxt is None:
                newNxt.add(nxt)
            elif nxt == targetNode:
                newNxt.add(headNode)
            else:
                newNxt.add(nxt)
                if nxt not in seen:
                    pending.append(nxt)
                    seen.add(nxt)
        node.nxt = newNxt

def printDAG(dag):
    """Print a DAG"""

    pending = [dag]
    seen = set(pending)
    while len(pending) > 0:
        node = pending.pop(0)
        for nxt in sorted(node.nxt):
            if nxt and nxt not in seen:
                seen.add(nxt)
                pending.append(nxt)

    # Work out the length of the longest shortStr
    maxLen = 0
    for node in seen:
        strLen = len(node.shortStr())
        if strLen > maxLen:
            maxLen = strLen

    formatStr = '%%%ds   %%s' % maxLen
    print formatStr % ('NODE', 'CHILDREN')
    for node in sorted(seen):
        childStrs = []
        for nxt in sorted(node.nxt):
            if nxt:
                childStrs.append(nxt.shortStr())
            else:
                childStrs.append('--')
        children = ', '.join(childStrs)
        print formatStr % (node.shortStr(), children)
        #nodeStr = '%s -> [' % node.shortStr()
        #firstChild = True
        #for nxt in sorted(node.nxt):
        #    if not firstChild:
        #        nodeStr += ', '
        #    if nxt:
        #        nodeStr += nxt.shortStr()
        #    else:
        #        nodeStr += str(nxt)
        #    firstChild = False
        #nodeStr += ']'
        #print nodeStr
    print ""

def printDAGParents(context):
    # Work out the length of the longest shortStr
    maxLen = 0
    for node in ccontext.dagOrderList:
        strLen = len(node.shortStr())
        if strLen > maxLen:
            maxLen = strLen

    formatStr = '%%%ds   %%s' % maxLen
    print formatStr % ('NODE', 'PARENT')
    for node in ccontext.dagOrderList:
        if len(ccontext.dagParents[node]) > 0:
            parentStr = ', '.join([parent.shortStr() for parent in sorted(ccontext.dagParents[node])])
        else:
            parentStr = '--'
        print formatStr % (node.shortStr(), parentStr)
        #for parent in ccontext.dagParents[node]:
        #    print parent.shortStr(),
        #print ""
def runExplore(context):
    baseNode = DAGChainNode(context.dag, 0, 0, 0)

    startTime = time.clock()
    results = exploreGraph(context, baseNode)
    endTime = time.clock()

    timeDelta = endTime - startTime

    return timeDelta

def runOpt(context):
    baseNode = DAGChainNode(context.dag, 0, 0, 0)

    startTime = time.clock()
    opt(context, baseNode, globalBPC)
    endTime = time.clock()

    (context.worstBits, context.worstCyc) = findOptNodes(context, baseNode, globalBPC)

    timeDelta = endTime - startTime

    return timeDelta

def runExploreAndOpt(context,
        showExpTime=False, showOptTime=False,
        showResultSet=False,
        printPathsToExplore=False, printResults=True,
        printEdges=True):

    expTime = runExplore(context)

    if showExpTime:
        print "Graph exploration time: %1.03fs" % expTime

    if printPathsToExplore:
        baseNode = DAGChainNode(context.dag, 0, 0, 0)
        (pathCount, combs) = context.exploreResultCount[baseNode]
        if pathCount > 1000000:
            exp = 0
            while pathCount > 1000000:
                exp += 3
                pathCount /= 1000
            pathCount /= 1000.0
            print "Total paths to explore: %1.3f x 10^%d" % (pathCount, exp)
        else:
            print "Total paths to explore: %d" % pathCount

    optTime = runOpt(context)

    if showOptTime:
        print "opt-algorithm runtime: %1.03fs" % optTime
        print "opt-algorithm evaluation steps: %d" % context.count

    if printResults:
        printBestOpt(context, printEdges)

    print "\n"

def combineClustersByInst(context):
    """
    Attempt to combine multiple clusters that only differ by instance number

    Requires counter logic in parser
    """
    if not hasattr(context, 'bestClusters'):
        return

    clustersByHdrPos = {}
    newEdgeCount = 0
    for cluster in sorted(context.bestClusters):
        if len(cluster.chain) > 0:
            zeroCluster = chainInstOne(context, cluster)
            zeroClusterEdgeCount = findEdgeCount(context, zeroCluster)

            #print cluster, zeroCluster
            firstCN = cluster.chain[0]
            firstHdrName = firstCN.dagNode.getName()
            firstPos = firstCN.startPos
            key = (firstHdrName, firstPos)

            if key not in clustersByHdrPos:
                clustersByHdrPos[key] = [zeroCluster]
                newEdgeCount += zeroClusterEdgeCount
            else:
                dupFound = False
                pos = 0
                while not dupFound and pos < len(clustersByHdrPos[key]):
                    superChain = compareChainIgnoreInst(zeroCluster, clustersByHdrPos[key][pos])
                    if superChain:
                        if superChain == zeroCluster:
                            prevEdgeCount = findEdgeCount(context, clustersByHdrPos[key][pos])
                            clustersByHdrPos[key][pos] = zeroCluster
                            newEdgeCount += zeroClusterEdgeCount - prevEdgeCount
                        dupFound = True
                    pos += 1

    context.bestClusters = set()
    for clusters in clustersByHdrPos.values():
        context.bestClusters.update(clusters)
    context.bestEdgeCount = newEdgeCount
    #for (cnode, cluster) in context.optNodes:
    #    if cluster not in context.bestClusters:
    #        context.bestEdgeCount += findEdgeCount(context, cluster)
    #    context.bestClusters.add(cluster)

def chainInstOne(context, chain):
    """
    Attempt to find an equivalent cluster indexed at position one
    """
    if len(chain.chain) == 0:
        return chain

    firstCN = chain.chain[0]
    dagNode = firstCN.dagNode
    if dagNode.inst == 1:
        return chain

    inst = dagNode.inst - 1
    bestCand = None
    while inst > 0:
        dagNodeStr = "%s:%d" % (dagNode.getName(), inst)
        if dagNodeStr in context.dagNodes:
            newDAGNode = context.dagNodes[dagNodeStr]
            newCN = DAGChainNode(newDAGNode, firstCN.startPos, 0, 0)
            clusters = findClustersAndCovers(context, newCN)

            for candidate in clusters:
                candidate = compareChainIgnoreInst(candidate, chain, True)
                if candidate:
                    bestCand = candidate

        inst -= 1

    if bestCand is not None:
        return bestCand

    return chain

def compareChainIgnoreInst(chainA, chainB, matchExactly=False):
    """
    Compare two chains to see if one is a subset of the other,
    ignoring instance numbers

    Return: superset chain if one chain is a superset of the other
            None if the two chains are independent
    """
    aLen = len(chainA.chain)
    bLen = len(chainB.chain)
    if aLen == 0 and bLen == 0:
        return chainA

    for pos in xrange(min(len(chainA.chain), len(chainB.chain))):
        aNode = chainA.chain[pos]
        aHdrName = aNode.dagNode.getName()
        aPos = aNode.startPos

        bNode = chainB.chain[pos]
        bHdrName = bNode.dagNode.getName()
        bPos = bNode.startPos

        if aHdrName != bHdrName or aPos != bPos:
            return None

    if matchExactly and (
        aLen != bLen or
        chainA.chain[-1].consumed != chainB.chain[-1].consumed or
        chainA.chain[-1].read != chainB.chain[-1].read):
        return None

    if aLen > bLen:
        return chainA
    elif bLen > aLen:
        return chainB
    elif chainA.chain[-1].consumed > chainB.chain[-1].consumed:
        return chainA
    elif chainB.chain[-1].consumed > chainA.chain[-1].consumed:
        return chainB
    elif chainA.chain[-1].read > chainB.chain[-1].read:
        return chainA
    elif chainB.chain[-1].read > chainA.chain[-1].read:
        return chainB

    return chainA

def updateDAGNodes(context):
    context.dagNodes = {}
    pending = [context.dag]
    seen = set(pending)
    while len(pending) > 0:
        node = pending.pop(0)
        dagNodeStr = "%s:%d" % (node.getName(), node.inst)
        context.dagNodes[dagNodeStr] = node
        
        for nxt in node.nxt:
            if nxt is not None and nxt not in seen:
                pending.append(nxt)
                seen.add(nxt)

def tryDAGBarrier(context, node):
    print "Procesing node '%s'   (Parents: %d)..." % (node, parentCount)
    posList = node.getDecisionBytes()

    # Trim temporarily
    if len(posList) > 0 and posList[-1] >= node.getLength():
        while len(posList) > 0 and posList[-1] >= node.getLength():
            posList.pop()
    elif node.getLength() not in posList:
        posList.append(node.getLength())
    posList.reverse()
    if 0 not in posList:
        posList.append(0)
    print "Barrier locations:", posList

    roundBestContext = context
    roundBestLoc = None
    for barPos in posList:
    #for barPos in [0,4,20,40]:
        print "Considering %d..." % barPos

        newContext = OptContext()
        newContext.dag = cloneDAG(context.dag)
        insertDAGBarrier(newContext.dag, node, barPos)
        updateDAGNodes(newContext)

        runExploreAndOpt(newContext,
                showExpTime=False, showOptTime=False,
                printPathsToExplore=False, printEdges=False)

        if optFound(newContext) and \
                newContext.bestEdgeCount < roundBestContext.bestEdgeCount:
            roundBestContext = newContext
            roundBestLoc = barPos

    return roundBestContext, roundBestLoc

def tryDAGParallelEdgeBarrier(context, node, nxtNode):
    edges = node.hdr.getDecisionCombos(
            nodeLen, nxtNode.hdr.name, nodeLen)
    print "Processing edges between %s and %s (edges: %d)" % (node, nxtNode, edges)

    posList = nxtNode.getDecisionBytes()
    # Trim temporarily
    if len(posList) > 0 and posList[-1] >= node.getLength():
        while len(posList) > 0 and posList[-1] > nxtNode.getLength():
            print "Popping:", posList.pop()
    elif nxtNode.getLength() not in posList:
        posList.append(nxtNode.getLength())
    posList.reverse()
    if 0 not in posList:
        posList.append(0)
    print "Barrier locations:", posList
    roundBestContext = context
    roundBestLoc = None
    for barPos in posList:
        print "Considering %d..." % barPos

        newContext = OptContext()
        newContext.dag = cloneDAG(context.dag)
        if barPos >= 0:
            insertDAGBarrier(newContext.dag, nxtNode, barPos)
        else:
            insertDAGBarrier(newContext.dag, node, nodeLen, nxtNode)
        updateDAGNodes(newContext)
        #printDAG(newContext.dag)

        runExploreAndOpt(newContext,
                showExpTime=False, showOptTime=False,
                printPathsToExplore=False, printEdges=False)
        #printEntries(newContext)

        if optFound(newContext) and \
                newContext.bestEdgeCount < roundBestContext.bestEdgeCount:
            roundBestContext = newContext
            roundBestLoc = barPos

    return roundBestContext, roundBestLoc

def allocateResultVectorEntries(context):
    global fieldResultSize

    # Verify that we have a dag order list
    if len(context.dagOrderList) == 0:
        sortDAG(context)

    hdrs = []
    hdrNum = 1
    fieldPos = 1
    context.foundHdrNum = {}
    context.fieldPos = {}
    for node in context.dagOrderList:
        if isinstance(node, HeaderNode) or \
                (isinstance(node, SubHeaderNode) and node.startPos == 0):
            nodeName = '%s-%d' % (node.hdr.name, node.inst)
            context.foundHdrNum[nodeName] = hdrNum
            hdrs.append(nodeName)

            extractBytes = node.hdr.getExtractBytes()
            fieldPos = extractPos[node.hdr.name] + (node.inst - 1) * extractOffset[node.hdr.name]
            extractBytes = zip(extractBytes, xrange(fieldPos, fieldPos + len(extractBytes)))
            context.fieldPos[nodeName] = extractBytes

            hdrNum += 1
            #fieldPos += len(extractBytes)

    context.hdrCount = hdrNum + 1
    context.maxDest = fieldPos + 1

    print "header-numbers:"
    for node in hdrs:
        print "  %s: %d" % (node, context.foundHdrNum[node])

if __name__ == "__main__":
    my_parser = argparse.ArgumentParser('Read headers from a given file')
    my_parser.add_argument('hdr_file', type=str,
            help='Parse graph description file')

    # General parameters
    my_parser.add_argument('--data-rate', type=int,
            default=dataRate, help='Data rate (Gbps)')
    my_parser.add_argument('--clk-freq', type=float,
            default=clkFreq, help='Clock frequency (Gbps)')
    my_parser.add_argument('--lookups', type=int,
            default=lookups, help='Lookups per cycle')
    my_parser.add_argument('--lookup-width', type=int,
            default=lookupWidth, help='Lookup width (bytes)')
    my_parser.add_argument('--max-skip', type=int,
            default=maxSkip, help='Maximum jump amount per cycle')
    my_parser.add_argument('--min-skip', type=int,
            default=minSkip, help='Minimum jump amount per cycle')
    my_parser.add_argument('--first-lookup-at-zero', action='store_true',
            default=firstLookupAtZero, help='First lookup is at zero offset')
    my_parser.add_argument('--no-first-lookup-at-zero', action='store_false',
            default=firstLookupAtZero, dest='first_lookup_at_zero',
            help='First lookup is NOT at zero offset')
    my_parser.add_argument('--window', type=int,
            default=windowSize, help='Size of window from which to extract lookups')
    my_parser.add_argument('--extract', action='store_true',
            default=extract, help='Table outputs fields to extract')
    my_parser.add_argument('--no-extract', action='store_false',
            default=extract, dest='extract',
            help='Table outputs header types only')
    my_parser.add_argument('--extract-bytes', type=int,
            default=extractBytes, help='Number of bytes to extract (extract only)')
    my_parser.add_argument('--max-hdrs', type=int,
            default=maxHdrs, help='Maximum number of headers identifiable in a single cycle (no-extract only)')

    my_parser.add_argument('--mult-parent', action='store_true',
            default=multParent, help='Try merging at nodes with mulitple parents')
    my_parser.add_argument('--no-mult-parent', action='store_false',
            default=multParent, dest='mult_parent',
            help='Do NOT try merging at nodes with multiple parents')

    my_parser.add_argument('--mult-parent-retry', action='store_true',
            default=multParentRetry, help='Retry failed multiple parent merging')
    my_parser.add_argument('--no-mult-parent-retry', action='store_false',
            default=multParentRetry, dest='mult_parent_retry',
            help='Do NOT retry failed multiple parent merging')

    my_parser.add_argument('--parallel-edge', action='store_true',
            default=parallelEdge, help='Try forcing barriers for parallel multiple edges')
    my_parser.add_argument('--no-parallel-edge', action='store_false',
            default=parallelEdge, dest='parallel_edge',
            help='Do NOT try forcing barriers for parallel multiple edges')

    my_parser.add_argument('--inst-merge', action='store_true',
            default=instMerge, help='Combine clusters that differ only by instance number')
    my_parser.add_argument('--no-inst-merge', action='store_false',
            default=instMerge, dest='inst_merge',
            help='Do NOT combine clusters that differ only by instance number')

    my_parser.add_argument('--ternary-state-match', action='store_true',
            default=ternMatchOnState, help='Perform ternary matching on state variables (single entry within a header)')
    my_parser.add_argument('--no-ternary-state-match', action='store_false',
            default=ternMatchOnState, dest='ternary_state_match',
            help='Do NOT perform ternary matching on state variables')

    my_parser.add_argument('--print-tcam', action='store_true',
            default=printTCAM, help='Print TCAM entries')
    my_parser.add_argument('--no-print-tcam', action='store_false',
            default=printTCAM, dest='print_tcam',
            help='Do NOT print TCAM entries')

    my_parser.add_argument('--save-tcam', action='store_true',
            default=saveTCAM, help='Print TCAM entries')
    my_parser.add_argument('--no-save-tcam', action='store_false',
            default=saveTCAM, dest='save_tcam',
            help='Do NOT save TCAM entries')

    my_parser.add_argument('--max-tcam-state', type=int,
            default=tcamMaxState, help='Maximum TCAM state value')

    my_parser.add_argument('--show-result-sets', action='store_true',
            default=showResultSets, help='Show all potential result sets ')

    args = my_parser.parse_args()


    hfile = args.hdr_file

    dataRate = args.data_rate
    clkFreq = args.clk_freq
    globalBPC = dataRate / clkFreq
    lookups = args.lookups
    lookupWidth = args.lookup_width
    maxSkip = args.max_skip
    minSkip = args.min_skip
    firstLookupAtZero = args.first_lookup_at_zero
    windowSize = args.window
    extract = args.extract
    extractBytes = args.extract_bytes
    extract = args.extract
    maxHdrs = args.max_hdrs

    multParent = args.mult_parent
    multParentRetry = args.mult_parent_retry
    parallelEdge = args.parallel_edge
    instMerge = args.inst_merge
    ternMatchOnState = args.ternary_state_match
    DAGChain.singleEntryForInternal = ternMatchOnState
    printTCAM = args.print_tcam
    saveTCAM = args.save_tcam

    showResultSets = args.show_result_sets

    tcamMaxState = args.max_tcam_state
    stateBits = int(math.ceil(math.log(tcamMaxState, 2)))
    stateBytes = int(math.ceil(stateBits / 8.0))
    stateWidth10 = int(math.ceil(math.log10(tcamMaxState)))

    if minSkip > 0:
        print "Stopping: minSkip parameter does not work correctly. Specified: %d" % minSkip
        sys.exit(1)

    printParams()

    (headerList, headers) = readHeaders(hfile)
    firstHdr = headerList[0]

    # Optimize the headers
    for hdr in headerList:
        if hdr.nextHeader is not None and type(hdr.nextHeader) == tuple:
            hdr.optNextHeader()

    ccontext.dag, ccontext.dagNodes = buildDAG(headerList, headers)
    printDAG(ccontext.dag)
    calcExtractLocs(headerList, headers)
    headNode = ccontext.dag

    sortDAG(ccontext)
    printDAGParents(ccontext)

    print "\n\n\n"
    print "Initial run:"
    runExploreAndOpt(ccontext, showExpTime=True, showOptTime=True, printPathsToExplore=True, printEdges=True)
    bestContext = ccontext

    if multParent:
        print "\nAttempting optimization of nodes with multiple parents...\n"

        barrierLocs = []
        unprocHeaders = []
        for node in reversed(ccontext.dagOrderList):
            parentCount = len(ccontext.dagParents[node])
            if parentCount <= 1:
                continue

            bestContext, bestLoc = tryDAGBarrier(bestContext, node)
            if bestLoc is not None:
                barrierLocs.append((node, bestLoc))
            else:
                unprocHeaders.append(node)

        if multParentRetry and len(unprocHeaders) > 0:
            print "Retrying:",
            for h in unprocHeaders:
                print h,
            print ""
            for node in unprocHeaders:
                bestContext, bestLoc = tryDAGBarrier(bestContext, node)
                if bestLoc is not None:
                    barrierLocs.append((node, bestLoc))

        print ""
        print "Multi-parent optimization results:"
        print "----------------------------------"
        printBestOpt(bestContext)
        print ""
        if len(barrierLocs) > 0:
            print "Barrier locations:"
            for (node, loc) in barrierLocs:
                print "  %s: %d" % (node, loc)
        else:
            print "Multi-parent barriers do not improve results"
        print "\n\n\n"

    if parallelEdge:
        print "\nAttempting optimization of parallel edges to downstream nodes...\n"
        # Insert barriers between adjacent nodes when multiple patterns cause the transition
        parallelLocs = []
        for node in ccontext.dagOrderList:
            if isinstance(node, HeaderNode):
                nodeLen = node.getLength()
                for nxt in node.nxt:
                    if nxt and node.hdr.getDecisionCombos(
                            nodeLen, nxt.hdr.name, nodeLen) > 1:
                        bestContext, bestLoc = tryDAGParallelEdgeBarrier(bestContext, node, nxt)
                        if bestLoc is not None:
                            parallelLocs.append((node, nxt, bestLoc))

        print ""
        print "Parallel edge optimization results:"
        print "-----------------------------------"
        printBestOpt(bestContext)
        print ""
        if len(parallelLocs) > 0:
            print "Parallel edge barrier locations:"
            for (node, nxtNode, loc) in parallelLocs:
                print "  %s -> %s: %d" % (node, nxtNode, loc)
        else:
            print "Parallel edge barriers do not improve results"
        print "\n\n\n"


    if instMerge:
        print "\nAttempting combining clusters that differ only by instance number...\n"
        combineClustersByInst(bestContext)

    print ""
    print "Final optimization results:"
    print "---------------------------"
    printDAG(bestContext.dag)
    printBestOpt(bestContext)
    print ""
    printEntries(bestContext)
    if printTCAM or saveTCAM:
        print ""
        allocateResultVectorEntries(bestContext)
        print ""
        printTCAMEntries(bestContext)
