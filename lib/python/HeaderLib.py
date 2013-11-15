#!/usr/bin/env python
#
# Library of header functions
#
# Useful functions:
#   getHeaderLengths: return the header lengths of a given header
#

from Header import Header
from Field import Field
from RefCount import RefCount
from pyparsing import *
import sys
import re

START = 0
DONE = 255

# Variables used in reading headers
shouldMergeFixedTransitions = True
shouldTrimNonReachable = True
wantWildcard = True
headerBNF = None


class HeaderInfo:
    """Simple class for returning header info"""
    def __init__(self, length, lenIsVar, matchBytes, match):
        self.length = length
        self.lenIsVar = lenIsVar
        self.matchBytes = matchBytes
        self.match = match

class HeaderInfoAll:
    """Simple class for returning header info for all length/match combos"""
    def __init__(self, lenIsVar, lenBytes, lenMatch, lengths, nxtHdrBytes, nxtHdrMatch, nxtHdrs, defNxtHdrVal):
        self.lenIsVar = lenIsVar
        self.lenBytes = lenBytes
        self.lenMatch = lenMatch
        self.lengths = lengths
        self.nxtHdrBytes = nxtHdrBytes
        self.nxtHdrMatch = nxtHdrMatch
        self.nxtHdrs = nxtHdrs
        self.defNxtHdrVal = defNxtHdrVal

hdrLengths = {}
def getHeaderLengths(hdr):
    """Get the header length(s) and the set of fields influencing the length"""
    global hdrLengths

    if hdr.name not in hdrLengths:

        # Work out the length of the current field -- may be variable!
        (minLen, maxLen, lenIsVariable, lenFields) = hdr.length()
        if lenIsVariable:
            (lenFields, lenFieldVals) = hdr.getLengthVarValues()

            # Calculate the set of valid lengths
            lengths = []
            for lenFieldValSet in lenFieldVals:
                valMap = {}
                for i in xrange(len(lenFields)):
                    valMap[lenFields[i]] = lenFieldValSet[i]
                lengths.append(hdr.doLengthCalc(valMap))

            # Get the list of bytes to extract for length field identification
            (lenFieldBytes, lenFieldPos, lenTotalWidth) = hdr.getFieldByteLocs(lenFields)
            lenContent = hdr.getFieldByteContents(lenFields, lenFieldVals, lenFieldBytes, lenFieldPos)
        else:
            (lenFields, lenFieldVals) = ([], [])
            lengths = [minLen]
            (lenFieldBytes, lenFieldPos) = ([], [])
            lenContent = []

        hdrLengths[hdr.name] =  (lenIsVariable, lengths, lenFields, lenFieldBytes, lenFieldPos, lenContent)

    return hdrLengths[hdr.name]

def getDecisionPos(hdr):
    """Get the first and last decision byte for the header"""
    decBytes = getAllDecisionBytes(hdr)

    if len(decBytes) > 0:
        return decBytes[0], decBytes[-1]
    else:
        return -1, -1

def getAllDecisionBytes(hdr):
    """Get all decision bytes for the header"""
    (lenIsVariable, lengths, lenFields, lenFieldBytes, lenFieldPos, lenContent) = getHeaderLengths(hdr)
    if hdr.nextHeader and isinstance(hdr.nextHeader, tuple):
        # Get the fields that are used in the lookup and the set of
        # values->next headers that they can take on
        fields = hdr.nextHeader[0]
        fieldMap = hdr.nextHeader[1]
        (fieldBytes, fieldPos, totalWidth) = hdr.getFieldByteLocs(fields) 

        # If we have a variable-length header, calculate a merger of the decision
        # fields with the length fields
        if lenIsVariable:
            (mergedFieldBytes, map1, map2) = hdr.mergeLocs(fieldBytes, lenFieldBytes)
        else:
            mergedFieldBytes = fieldBytes

        return mergedFieldBytes
    elif lenIsVariable:
        return lenFieldBytes
    else:
        return []

def exploreHeader(hdr, headers, callback, callPerLenNxtHdr = True):
    """Explore the headers reachable from the given header"""
    # Get the set of lengths and associated variables for this header
    (lenIsVariable, lengths, lenFields, lenFieldBytes, lenFieldPos, lenContent) = getHeaderLengths(hdr)
    lengths = [x / 8 for x in lengths]

    hdr.incRefCount()
    if hdr.refCountExceedsLimit():
        hdr.decRefCount()
        return
    
    # Walk through the combination of lengths and next headers
    # (assume that headers have been merged)
    if hdr.nextHeader:
        # Case: one of N next headers

        # Get the fields that are used in the lookup and the set of
        # values->next headers that they can take on
        fields = hdr.nextHeader[0]
        fieldMap = hdr.nextHeader[1]
        (fieldBytes, fieldPos, totalWidth) = hdr.getFieldByteLocs(fields) 

        if callPerLenNxtHdr:
            # If we have a variable-length header, calculate a merger of the decision
            # fields with the length fields
            if lenIsVariable:
                (mergedFieldBytes, map1, map2) = hdr.mergeLocs(fieldBytes, lenFieldBytes)
            else:
                mergedFieldBytes = fieldBytes

            # Walk through all of the decision value combinations
            for (vals, nextHeaderName) in fieldMap:
                if nextHeaderName:
                    nxtHdr = headers[nextHeaderName]
                else:
                    nxtHdr = None

                fieldMatch = hdr.getFieldByteContentsSingle(fields, vals, fieldBytes, fieldPos)

                # Walk through all of the possible lengths
                for i in xrange(len(lengths)):
                    length = lengths[i]
                    if lenIsVariable:
                        match = hdr.mergeContents(mergedFieldBytes, map1, map2, fieldMatch, lenContent[i])
                    else:
                        match = fieldMatch

                    hdrInfo = HeaderInfo(length, lenIsVariable, mergedFieldBytes, match)
                    callback(hdr, nxtHdr, headers, hdrInfo)
        else:
            nxtHdrMatch = []
            nxtHdrs = []

            # Walk through all of the decision value combinations
            for (vals, nextHeaderName) in fieldMap:
                if nextHeaderName:
                    nxtHdr = headers[nextHeaderName]
                else:
                    nxtHdr = None

                fieldMatch = hdr.getFieldByteContentsSingle(fields, vals, fieldBytes, fieldPos)

                nxtHdrMatch.append(fieldMatch)
                nxtHdrs.append(nxtHdr)

            defNxtHdrVal = None
            if hdr.getDefNxtHdrVal() is not None:
                defNxtHdrVal = hdr.getFieldByteContentsSingle(fields, hdr.getDefNxtHdrVal(), fieldBytes, fieldPos)

            hdrInfo = HeaderInfoAll(lenIsVariable, lenFieldBytes, lenContent,
                    lengths, fieldBytes, nxtHdrMatch, nxtHdrs, defNxtHdrVal)
            callback(hdr, headers, hdrInfo)

    else:
        # Case: no next header
        if callPerLenNxtHdr:
            for i in xrange(len(lengths)):
                length = lengths[i]
                if lenIsVariable:
                    match = lenContent[i]
                else:
                    match = None

                hdrInfo = HeaderInfo(length, lenIsVariable, lenFieldBytes, match)
                callback(hdr, None, headers, hdrInfo)
        else:
            hdrInfo = HeaderInfoAll(lenIsVariable, lenFieldBytes, lenContent,
                    lengths, None, None, [], None)
            callback(hdr, headers, hdrInfo)

    hdr.decRefCount()

if __name__ == "__main__":
    import argparse
    from read_headers import readHeaders
    import copy

    hfile = 'headers.txt'

    seenPath = set()
    finalPaths = set()
    def exploreHdrChain(hdr, nxtHdr, headers, hdrInfo, chain = None):
        if chain:
            path = chain + " -> "
        else:
            path = ""
        path += "%s (%d)" % (hdr.name, hdrInfo.length)
        if nxtHdr:
            nxtHdrName = nxtHdr.name
        else:
            nxtHdrName = "--"
        pathNxt = path + " -> " + nxtHdrName
        if pathNxt not in seenPath:
            seenPath.add(pathNxt)
            if nxtHdr:
                ehcFunc = lambda ehcHdr, ehcNxtHdr, ehcHeaders, ehcInfo : \
                        exploreHdrChain(ehcHdr, ehcNxtHdr, ehcHeaders, ehcInfo, path)
                exploreHeader(nxtHdr, headers, ehcFunc)
            else:
                finalPaths.add(path)
    
    
    def exampleHeaderWalk(hdr, headers):
        exploreHeader(headerList[0], headers, exploreHdrChain)

        for path in sorted(finalPaths):
            print path


    seenHdrs = set()
    def exploreHdrChainMerged(hdr, headers, hdrInfo):
        if not hdr.name in seenHdrs:
            print hdr.name
            seenHdrs.add(hdr.name)

            for nxtHdr in hdrInfo.nxtHdrs:
                if nxtHdr:
                    exploreHeader(nxtHdr, headers, exploreHdrChainMerged, False)


    def exampleHeaderWalkMerged(hdr, headers):
        exploreHeader(headerList[0], headers, exploreHdrChainMerged, False)


    my_parser = argparse.ArgumentParser('Read headers from a given file')
    my_parser.add_argument('--hdr_file', metavar='hdr_file', type=str,
            default=hfile, help='Header description file')
    args = my_parser.parse_args()

    (headerList, headers) = readHeaders(args.hdr_file)
    print "List of headers (from walk)"
    print "==========================="
    exampleHeaderWalkMerged(headerList[0], headers)
    print "\n"

    print "Header paths"
    print "============"
    exampleHeaderWalk(headerList[0], headers)

def crackKey(hdr, key, fields):
    numFields = len(fields)
    mask = []
    data = []

    # Identify the widths of the various fields
    widths = hdr.getFieldWidths(fields)
    widths.reverse()

    if isinstance(key, int):
        for width in widths:
            fieldMask = 2 ** width - 1
            fieldData = key & fieldMask
            mask.append(fieldMask)
            data.append(fieldData)
            key >>= width
    else:
        for width in widths:
            fieldMask = key[-width:]
            fieldMask = fieldMask.replace('0', '1').replace('x', '0')
            fieldMask = int(fieldMask, 2)

            fieldData = key[-width:]
            fieldData = fieldData.replace('x', '0')
            fieldData = int(fieldData, 2)

            mask.append(fieldMask)
            data.append(fieldData)
            key = key[:-width]

    mask.reverse()
    data.reverse()

    return (mask, data)

def getHeaderBNF():
    global headerBNF

    if not headerBNF:
        LBRACE, RBRACE, COLON, COMMA, EQ, LPAREN, RPAREN = map(Suppress, "{}:,=()")
        identifier = Word(alphas,alphanums+'_'+'-')
        identifierOrStar = Or(Literal('*'), identifier)
        integer = Word(nums)
        integerOrStar = Literal('*') ^ integer
        hexval = Combine(Literal('0x') + Word(hexnums))
        binary = Combine(Literal('b') + Word('01'))
        ternary = Combine(Literal('b') + Word('01x'))

        intOrHex = hexval ^ integer
        intOrHexOrBin= hexval ^ integer ^ binary
        intOrHexOrTern= hexval ^ integer ^ ternary

        extract = CaselessLiteral('extract')
        shiftOp = Literal('<<') ^ Literal('>>')
        addSubOp = Literal('+') ^ Literal('-')
        mulOp = Literal('*')
        identifierOrInteger = integer ^ identifier

        #expression = Forward()
        #atom = identifierOrInteger | (LPAREN + expression + RPAREN)
        #mulExp = Group(atom + OneOrMore(mulOp +
        #    atom)) | atom
        #addSubExp = Group(mulExp + OneOrMore(addSubOp + mulExp)) ^ mulExp
        #expression << (Group(addSubExp + OneOrMore(shiftOp + addSubExp)) ^ addSubExp)
        expression = Forward()
        atom = identifierOrInteger | (LPAREN + expression + RPAREN)
        mulExp = (atom + OneOrMore(mulOp +
            atom)) | atom
        addSubExp = (mulExp + OneOrMore(addSubOp + mulExp)) ^ mulExp
        expression << ((addSubExp + OneOrMore(shiftOp + addSubExp)) ^ addSubExp)

        fieldValue = Group(identifier('name') + COLON + integerOrStar('value') + Optional(COLON + extract('extract')))
        fieldList = fieldValue + ZeroOrMore(COMMA + fieldValue) + Optional(COMMA)
        fields = Suppress('fields') + LBRACE + fieldList + RBRACE
        pseudofields = Suppress('pseudo-fields') + LBRACE + fieldList + RBRACE


        mapKeys = Group(intOrHexOrTern + ZeroOrMore(COMMA + intOrHexOrTern))
        mapEntry = Group(mapKeys('keys') + COLON + identifier('value'))
        mapList = mapEntry + ZeroOrMore(COMMA + mapEntry) + Optional(COMMA)

        from_header = Group(identifier + ZeroOrMore(COMMA + identifier))
        map_start = Suppress('map') + LPAREN + from_header('from_header') + RPAREN
        map_body = LBRACE + Optional(mapList('maplist')) + RBRACE
        mapTable = Group(map_start + map_body)
        next_header = Suppress('next_header') + EQ + (mapTable('mapping') ^ identifier('field'))
        next_header_def = Suppress('next_header_def') + EQ + identifier

        max_var = Suppress('max_var') + EQ + identifier
        max_val = Suppress('max') + EQ + integer

        #length = Suppress('length') + EQ + expression
        length = Suppress('length') + EQ + Group(expression)

        max_length = Suppress('max_length') + EQ + integer

        singleHeader = Group(identifier('hdr') + LBRACE + Optional(fields('fields')) + Optional(pseudofields('pseudofields')) + Optional(next_header('next_header')) + Optional(next_header_def('next_header_def')) + Optional(max_var('maxvar')) + Optional(max_val('maxval')) + Optional(length('hdr_len')) + Optional(max_length('max_len')) + RBRACE)

        headerBNF = ZeroOrMore(singleHeader)

        comment = Literal('#') + Optional(restOfLine)

        headerBNF.ignore(comment)

    return headerBNF

def readHeaders(filename):
    """Read all of the headers from a file"""

    fh = open(filename)
    data = fh.read()
    fh.close()
    #data = filename

    parser = getHeaderBNF()

    intRE = re.compile(r'^\d+$')
    opRE = re.compile(r'^[+\-*]|<<|>>$')

    refCounts = {}
    headerList = []
    headers = {}
    for item in parser.parseString(data, True):
        #print item

        #print item.hdr
        if item.hdr not in headers:
            hdr = Header(item.hdr)
            headerList.append(hdr)
            headers[item.hdr] = hdr
            if item.fields != '':
                for fieldData in item.fields:
                    (name, width) = fieldData[0:2]
                    if width == '*':
                        width = None
                    else:
                        width = int(width)
                    hdr.addField(name, width)
                    if len(fieldData) == 3:
                        hdr.addExtractField(name)

            if item.pseudofields != '':
                for (name, width) in item.pseudofields:
                    if width == '*':
                        width = None
                    else:
                        width = int(width)
                    hdr.addPseudofield(name, width)

            if item.next_header != '':
                if item.next_header.field != '':
                    hdr.setNextHeader(str(item.next_header.field))
                else:
                    from_fields = item.next_header.mapping.from_header.asList()
                    rangeCount = 0
                    widths = hdr.getFieldWidths(from_fields)
                    hdrMap = {}
                    for (keys, value) in item.next_header.mapping.maplist:
                        for key in keys:
                            if key.find('0x') == 0:
                                key = int(key, 16)
                            elif key.find('b') == 0:
                                pass
                            else:
                                key = int(key)

                            (mask, data) = crackKey(hdr, key, from_fields)

                            hdrMap[key] = ((mask, data), value)

                            # Approximate counting of the number of entries covered
                            wildcards = 0
                            for i in xrange(len(widths)):
                                fieldWidth = widths[i]
                                fieldMask = mask[i]
                                wildcards += sum([(~fieldMask >> shift) & 1 for shift in xrange(fieldWidth)])
                            rangeCount += 2 ** wildcards

                    # Attempt to sort the header map by key
                    keys = hdrMap.keys()
                    keys.sort()
                    hdrList = []
                    for key in keys:
                        hdrList.append(hdrMap[key])

                    # Add a wildcard entry if we haven't covered all inputs
                    if rangeCount < 2 ** sum(widths) and wantWildcard:
                        hdrList.append((([0] * len(from_fields), [0] * len(from_fields)), None))

                    hdr.setNextHeader((from_fields, hdrList))

            if item.next_header_def != '':
                if item.next_header == '' or item.next_header.field != '':
                    raise Exception("next_header_def value specified but next_header is not specified/is not a map")
                from_fields = item.next_header.mapping.from_header.asList()
                (mask, data) = crackKey(hdr, item.next_header_def[0], from_fields)
                hdr.setDefNxtHdrVal((mask, data))

            if item.maxvar != '':
                item.maxvar = str(item.maxvar)
                if item.maxvar not in refCounts:
                    refCounts[item.maxvar] = RefCount(item.maxvar)
                hdr.setRefCount(refCounts[item.maxvar])

            if item.maxval != '':
                if not hdr.getRefCount():
                    if item.hdr not in refCounts:
                        refCounts[item.hdr] = RefCount(item.hdr)
                    hdr.setRefCount(refCounts[item.hdr])
                hdr.getRefCount().setMax(int(item.maxval[0]))

            if item.hdr_len != '':
                newExp = []
                if type(item.hdr_len[0]) == str:
                    newExp.append(item.hdr_len[0])
                else:
                    #print item.hdr_len[0]
                    #print type(item.hdr_len[0])
                    #print len(item.hdr_len[0])
                    #print len(item.hdr_len.asList())
                    for exp in item.hdr_len[0]:
                        if intRE.match(exp):
                            exp = int(exp)
                            newExp.append(exp)
                        elif opRE.match(exp):
                            newExp.append(exp)
                        else:
                            newExp.append(exp)
                hdr.setCalcLength(newExp)
                
                #hdr.setMax(int(item.maxval[0]))

            if item.max_len != '':
                hdr.setMaxLength(int(item.max_len[0]))
                
            #print hdr
        else:
            print "Error: header '%s' seen multiple times" % item.hdr
            sys.exit(-1)
        #print hdr

    # Merge fixed/deterministic transitions if requested
    if shouldMergeFixedTransitions:
        headerList = mergeTransitions(headerList, headers)

    # Work out which headers are reachable
    if shouldTrimNonReachable:
        headerList = trimNonReachable(headerList, headers)

    return (headerList, headers)

def mergeTransitions(headerList, headers):
    mergedHeaderList = []
    remapHdrs = {}
    for hdr in headerList:
        origHdrName = hdr.name
        if hdr.nextHeader:
            pos = 1
            while isinstance(hdr.nextHeader, str):
                nextHdr = headers[hdr.nextHeader]

                mergedHeader = Header('%s+%s' % (hdr.name, nextHdr.name))
                for field in hdr.fieldList:
                    mergedHeader.addField(field.name, field.width)
                for field in nextHdr.fieldList:
                    mergedHeader.addField('%d-%s' % (pos, field.name), field.width)

                for field in hdr.pseudofieldList:
                    mergedHeader.addPseudofield(field.name, field.width)
                for field in nextHdr.pseudofieldList:
                    mergedHeader.addPseudofield('%d-%s' % (pos, field.name), field.width)

                mergedHeader.nextHeader = nextHdr.nextHeader
                if isinstance(mergedHeader.nextHeader, tuple):
                    mergedHeader.nextHeader = copy.deepcopy(mergedHeader.nextHeader)
                    from_fields = mergedHeader.nextHeader[0]
                    for i in xrange(len(from_fields)):
                        from_fields[i] = '%d-%s' % (pos, from_fields[i])
                    #print from_fields

                if hdr.getRefCount() and nextHdr.getRefCount():
                    print "Headers being merged both have a reference count. Exiting..."
                    sys.exit(-1)
                elif hdr.getRefCount():
                    mergedHeader.setRefCount(hdr.getRefCount())
                else:
                    mergedHeader.setRefCount(nextHdr.getRefCount())

                if hdr.calcLength:
                    print "First header being merged has a calculated length. Unable to merge..."
                    sys.exit(-1)
                elif nextHdr.calcLength:
                    newExp = [hdr.length()[0], '+'] + nextHdr.calcLength
                    for i in xrange(len(newExp)):
                        if isinstance(newExp[i], str) and not opRE.match(newExp[i]):
                            newExp[i] = '%d-%s' % (pos, newExp[i])
                    mergedHeader.setCalcLength(newExp)

                if hdr.maxLength:
                    print "First header being merged has a max length. Unable to merge..."
                    sys.exit(-1)
                elif nextHdr.maxLength:
                    mergedHeader.setMaxLength(hdr.length()[0] + nextHdr.maxLength)

                remapHdrs[origHdrName] = mergedHeader.name
                hdr = mergedHeader
                headers[hdr.name] = hdr
                pos += 1

        mergedHeaderList.append(hdr)

    # Update references to the new headers
    for hdr in mergedHeaderList:
        if hdr.nextHeader:
            if isinstance(hdr.nextHeader, tuple):
                hdrList = hdr.nextHeader[1]
                for i in xrange(len(hdrList)):
                    nxtHdrName = hdrList[i][1]
                    if nxtHdrName in remapHdrs:
                        hdrList[i] = (hdrList[i][0], remapHdrs[nxtHdrName])

    return mergedHeaderList

def trimNonReachable(headerList, headers):
    reachable = []
    reachableNames = set()
    reachableHeaderList = []
    if len(headerList) > 0:
        reachable.append(headerList[0])
    while len(reachable) > 0:
        nextReachable = []
        while len(reachable) > 0:
            hdr = reachable.pop(0)
            if hdr.name not in reachableNames:
                reachableHeaderList.append(hdr)
                reachableNames.add(hdr.name)
                if hdr.nextHeader:
                    if isinstance(hdr.nextHeader, str):
                        nxtHdr = headers[hdr.nextHeader]
                        if nxtHdr not in reachableNames:
                            nextReachable.append(nxtHdr)
                    elif isinstance(hdr.nextHeader, tuple):
                        for (maskData, nxtHdrName) in hdr.nextHeader[1]:
                            if nxtHdrName and nxtHdrName not in reachableNames:
                                nextReachable.append(headers[nxtHdrName])

        reachable = nextReachable

    return reachableHeaderList

def setMergeFixedTransitions(merge):
    global shouldMergeFixedTransitions
    shouldMergeFixedTransitions = merge

def setTrimNonReachable(trim):
    global shouldTrimNonReachable
    shouldTrimNonReachable = trim

def setWantWildcard(want):
    global wantWildcard
    wantWildcard = want


# Test code
hfile = 'headers.txt'
if __name__ == "__main__":
    my_parser = argparse.ArgumentParser('Read headers from a given file')
    my_parser.add_argument('--hdr_file', metavar='hdr_file', type=str,
            default=hfile, help='Header description file')
    args = my_parser.parse_args()

    (headerList, headers) = readHeaders(args.hdr_file)

    print "Headers:",
    for hdr in headerList:
        print hdr.name,
    print ""

    for hdr in headerList:
        (l, m, o, oList) = hdr.length()

        if not o:
            print "%s:%d" % (hdr.name, l)
        else:
            if m:
                print "%s:%d+* (max: %d)" % (hdr.name, l, m)
                print "  Header fields: %s" % str(hdr.getLengthVarValues())
            else:
                print "%s:%d+*" % (hdr.name, l)
