#!/usr/bin/env python

from Header import Header
from HeaderLib import exploreHeader, readHeaders, DONE
from Field import Field
import string
import sys
import re
import argparse
from pyparsing import *
import json

hdrStarts = {}
hdrSeqs = []
hdrSeqStrs = set()
seenHdrs = set()
orderedHdrs = []
simpleHdrs = {}

def resetVars():
    global hdrStarts, hdrSeqs, hdrSeqStrs, seenHdrs, orderedHdrs, simpleHdrs
    hdrStarts = {}
    hdrSeqs = []
    hdrSeqStrs = set()
    seenHdrs = set()
    orderedHdrs = []
    simpleHdrs = {}

def matchToVerilog(match):
    """Concert a match array pair to a verilog string"""
    matchStr = "%d'b" % (len(match[0]) * 8)
    for j in xrange(len(match[0])):
        maskByte = match[0][j]
        dataByte = match[1][j]
        for k in xrange(8):
            if maskByte & (0x80 >> k):
                matchStr += "%d" % ((dataByte & (0x80 >> k)) >> (7 - k))
            else:
                matchStr += "?"
    return matchStr

def processHdrStart(hdr, nxtHdr, headers, hdrInfo, pos = 0, hdrSeq = []):
    if hdr.name not in hdrStarts:
        hdrStarts[hdr.name] = set()
    hdrStarts[hdr.name].add(pos)

    hdrSeq = hdrSeq + [hdr.name, pos, hdrInfo.length]
    if nxtHdr:
        phsFunc = lambda phsHdr, phsNxtHdr, phsHeaders, phsInfo : \
                processHdrStart(phsHdr, phsNxtHdr, phsHeaders, phsInfo,
                        pos + hdrInfo.length, hdrSeq)
        exploreHeader(nxtHdr, headers, phsFunc)
    else:
        hdrSeqStr = ""
        for i in xrange(0, len(hdrSeq), 3):
            hdrSeqStr += "%s:%d:%d==" % (hdrSeq[i], hdrSeq[i+1], hdrSeq[i+2])
        if hdrSeqStr not in hdrSeqStrs:
            hdrSeqStrs.add(hdrSeqStr)
            hdrSeqs.append(hdrSeq)

def generateJSONObject(hdr, headers, hdrInfo):
    """Porocess a single header in preparation for JSONification"""
    #if hdr.name == 'ieee802-1q':
    #    print hdr
    #    print hdr.refCount
    #    print type(hdr.refCount.name)
    #    print hdr.refCount.name
    #    print len(hdr.refCount.name)
    #    sys.exit(1)
    if hdr.name not in seenHdrs:
        seenHdrs.add(hdr.name)
        orderedHdrs.append(hdr.name)

        #print hdr.name

        hdrDict = {}
        if hdrInfo.lenIsVar:
            hdrDict['len'] = 0
            hdrDict['lenBytes'] = hdrInfo.lenBytes
            lenMap = {}
            for i in xrange(len(hdrInfo.lengths)):
                lenMatch = matchToVerilog(hdrInfo.lenMatch[i])
                length = hdrInfo.lengths[i]
                lenMap[lenMatch] = length
            hdrDict['lenMap'] = lenMap
        else:
            hdrDict['len'] = hdrInfo.lengths[0]

        if len(hdrInfo.nxtHdrs) > 0:
            hdrDict['nxtHdrBytes'] = hdrInfo.nxtHdrBytes
            nxtHdrMap = {}
            for i in xrange(len(hdrInfo.nxtHdrs)):
                nxtHdrMatch = matchToVerilog(hdrInfo.nxtHdrMatch[i])
                nxtHdr = hdrInfo.nxtHdrs[i]
                if nxtHdr:
                    nxtHdrMap[nxtHdrMatch] = nxtHdr.name
                else:
                    nxtHdrMap[nxtHdrMatch] = "DONE"
            hdrDict['nxtHdrMap'] = nxtHdrMap

        hdrDict['starts'] = sorted(hdrStarts[hdr.name])

        hdrDict['extractFields'] = []
        for fieldName in hdr.extractFields:
            field = hdr.fields[fieldName]
            hdrDict['extractFields'].append(fieldName)
            hdrDict['extractFields'].append(field.pos)
            hdrDict['extractFields'].append(field.width)

        if hdrInfo.defNxtHdrVal is not None:
            hdrDict['defNxtHdrVal'] = matchToVerilog(hdrInfo.defNxtHdrVal)

        if hdr.refCount:
            hdrDict['refCountName'] = hdr.refCount.name

        simpleHdrs[hdr.name] = hdrDict

        for nxtHdr in hdrInfo.nxtHdrs:
            if nxtHdr:
                exploreHeader(nxtHdr, headers, generateJSONObject, False)

def generateJSON(headerList, headers, dstFile = None):
    """Print out a JSON version of the header list"""
    resetVars()

    exploreHeader(headerList[0], headers, processHdrStart, True)
    exploreHeader(headerList[0], headers, generateJSONObject, False)

    if dstFile is None:
        print json.dumps(indent = 4, obj =
                {
                    'firstHdr': orderedHdrs[0],
                    'headers': simpleHdrs,
                    'hdrSeqs': hdrSeqs,
                })
    else:
        f = open(dstFile, 'w')
        json.dump(indent = 4, fp = f, obj =
                {
                    'firstHdr': orderedHdrs[0],
                    'headers': simpleHdrs,
                    'hdrSeqs': hdrSeqs,
                })
        f.close()

hfile = 'headers.txt'
if __name__ == "__main__":
    my_parser = argparse.ArgumentParser('Read headers from a given file')
    my_parser.add_argument('hdr_file', type=str,
            help='Parse graph description file')
    my_parser.add_argument('json_file', type=str,
            help='Output JSON file')
    args = my_parser.parse_args()

    #if args.json_file:
    #    json_file = args.json_file
    #else:
    #    json_file = args.hdr_file
    #    if json_file.endswith(".txt"):
    #        json_file = json_file[:-4] + ".json"
    #print args.hdr_file, json_file

    jsonFile = args.json_file
    if jsonFile == '-':
        jsonFile = None

    (headerList, headers) = readHeaders(args.hdr_file)
    generateJSON(headerList, headers, jsonFile)
