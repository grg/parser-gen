#!/usr/bin/env python

"""A header found within a packet"""

import Field
import RefCount
import copy
import re

DEFAULT = 'DEFAULT'
MAX_HDR_LEN = 512
ANY = 'ANY'

opre = re.compile(r'^[+\-*/]$')

class Header():
    """Header found within a packet"""
    
    def __init__(self, name):
        self.name = name
        self.fields = {}
        self.fieldList = []
        self.pseudofields = {}
        self.pseudofieldList = []
        self.extractFields = []
        self.nextHeader = None
        self.calcLength = None
        self.maxLength = None
        self.refCount = None
        self.lastPos = 0
        self.defNxtHdrVal = None
        self._decCombos = {}
        self._length = None
        
        self.hdrLengths = None

    def addField(self, name, width = None):
        if name not in self.fields:
            field = Field.Field(name, width, self.lastPos)
            self.fields[name] = field
            self.fieldList.append(field)
            if width:
                self.lastPos += width
        else:
            pass

    def addExtractField(self, name):
        if name not in self.extractFields:
            self.extractFields.append(name)
        else:
            pass

    def addPseudofield(self, name, width = None):
        if name not in self.pseudofields:
            field = Field.Field(name, width)
            self.pseudofields[name] = field
            self.pseudofieldList.append(field)
        else:
            pass

    def setRefCount(self, refCount):
        self.refCount = refCount

    def getRefCount(self):
        return self.refCount

    def setDefNxtHdrVal(self, defNxtHdrVal):
        self.defNxtHdrVal = defNxtHdrVal

    def getDefNxtHdrVal(self):
        return self.defNxtHdrVal

    def setNextHeader(self, nextHeader):
        self.nextHeader = nextHeader

    def setMaxLength(self, maxLength):
        self.maxLength = maxLength

    def length(self):
        if self._length is None:
            minLen = 0
            maxLen = None
            optional = False
            optionalFields = []
            for field in self.fieldList:
                if field.width >= 0:
                    minLen += field.width
                else:
                    optional = True
                    optionalFields.append(field)
            if self.maxLength:
                maxLen = self.maxLength

            self._length = (minLen, maxLen, optional, optionalFields)

        return self._length

    def setCalcLength(self, calcLength):
        self.calcLength = calcLength

    def doLengthCalc(self, valMap):
        if not self.calcLength:
            return self.length()[0]

        calcStr = ''
        for tok in self.calcLength:
            if isinstance(tok, int):
                calcStr += str(tok)
            elif opre.match(tok):
                calcStr += tok
            else:
                calcStr += str(valMap[tok])

        return eval(calcStr)

    def __str__(self):
        fieldList = ''
        for field in self.fieldList:
            if fieldList != '':
                fieldList += ', '
            if field.width >= 0:
                fieldList += '%s:%d' % (field.name, field.width)
            else:
                fieldList += '%s:?' % (field.name)
        fieldList = 'fields={%s}' % fieldList

        length = ''
        (hLen, hMaxLen, hOpt, hOptList) = self.length()
        if hOpt:
            expr = ''
            if self.calcLength:
                for exp in self.calcLength:
                    if expr != '':
                        expr += ' '
                    expr += str(exp)
            else:
                expr = '?'

            if hMaxLen:
                expr += ', max = %d' % hMaxLen

            length = 'length=%d+* (%s)' % (hLen, expr)
        else:
            length = 'length=%d' % hLen

        optional = ''
        if self.nextHeader:
            if type(self.nextHeader) == str:
                optional += ", nextHeader='%s'" % self.nextHeader
            elif isinstance(self.nextHeader, Field.Field):
                optional += ", nextHeader='%s'" % self.nextHeader.name
            elif type(self.nextHeader) == tuple:
                name = ''
                if isinstance(self.nextHeader[0], str):
                    name = self.nextHeader[0]
                elif isinstance(self.nextHeader[0], list):
                    for hdr in self.nextHeader[0]:
                        if name != '':
                            name += ', '
                        name += hdr
                else:
                    name = self.nextHeader[0].name
                optional += ", nextHeader='%s'->%s" % (name,
                        str(self.nextHeader[1]))

        if self.refCount:
            optional += ', refCount=%s' % self.refCount

        return "Header: [name='%s', %s, %s%s]" % (self.name, fieldList,
                length, optional)

    def incRefCount(self):
        if self.refCount and isinstance(self.refCount, RefCount.RefCount):
            self.refCount.inc()

    def decRefCount(self):
        if self.refCount and isinstance(self.refCount, RefCount.RefCount):
            self.refCount.dec()

    def refCountAtLimit(self):
        if self.refCount and isinstance(self.refCount, RefCount.RefCount):
            return self.refCount.atLimit()
        else:
            return False

    def refCountExceedsLimit(self):
        if self.refCount and isinstance(self.refCount, RefCount.RefCount):
            return self.refCount.exceedsLimit()
        else:
            return False

    def decisionLoc(self):
        """Get the last byte used in the lookup decision process"""
        if not self.nextHeader or isinstance(self.nextHeader, str):
            return 0

        if isinstance(self.nextHeader, tuple):
            decPos = 0
            fields = self.nextHeader[0]

            # Identify the bytes that need to be extracted
            pos = 0
            for field in self.fieldList + self.pseudofieldList:
                if field.width > 0:
                    if field.name in fields:
                        endPos = pos + field.width
                        decPos = int((endPos - 1) / 8) + 1

                    pos += field.width

            return decPos

        return 0

    def lengthLoc(self):
        """Get the last byte used in the length calculation process"""
        if not self.calcLength:
            return 0

        lengthFields = self.getLengthFields()
        decPos = 0

        # Identify the bytes that need to be extracted
        pos = 0
        for field in self.fieldList + self.pseudofieldList:
            if field.width > 0:
                if field.name in lengthFields:
                    endPos = pos + field.width
                    decPos = int((endPos - 1) / 8) + 1

                pos += field.width

        return decPos

    def decisionLengthLoc(self):
        """Get the last byte used in the lookup decision/length process"""
        return max(self.decisionLoc(), self.lengthLoc())

    def check(self):
        """Check validity of header

        Requirements:
          - integer number of bytes
        """
        # Check length
        minLen = self.length()[0]
        if minLen % 8 != 0:
            return False

        return True

    def getFieldByteLocs(self, fields):
        """Get a list of bytes corresponding to a set of fields

        Params:
          fields - list of fields to get bytes of

        Return:
          (fieldBytes, fieldPos, fieldWidth)
            - fieldBytes - list of bytes to extract from header
            - fieldPos   - dict of field name -> (first byte, first bit in extracted header)
            - totalWidth - total width of extracted fields
        """
        fieldBytes = []
        fieldPos = {}
        totalWidth = 0

        # Identify the bytes that need to be extracted
        pos = 0
        for field in self.fieldList + self.pseudofieldList:
            if field.width > 0:
                if field.name in fields:
                    endPos = pos + field.width
                    firstByte = int(pos / 8)
                    lastByte = int((endPos - 1) / 8)

                    if len(fieldBytes) > 0 and fieldBytes[-1] == firstByte:
                        fieldBytes += range(firstByte + 1, lastByte + 1)
                    else:
                        fieldBytes += range(firstByte, lastByte + 1)

                    fieldPos[field.name] = (len(fieldBytes) - (lastByte - firstByte) - 1, pos % 8)

                    totalWidth += field.width
                pos += field.width

        return (fieldBytes, fieldPos, totalWidth)

    def getFieldByteContentsSingle(self, fields, vals, fieldBytes, fieldPos):
        """Get a list of field byte contents for a single set of values

        Params:
          fields     - list of fields to get byte contents of
          vals       - list of values for corresponding bytes or
                       pair of list of values of mask/data bytes
          fieldBytes - list of byte locaitons of the fields
          fieldPos   - map of field name and first byte/bit in *extracted* bytes

        Return:
          (contentMask, contentData)
            - contentMask - mask values of content bytes
            - contentData - data values of content bytes
        """
        numFields = len(fields)
        numBytes = len(fieldBytes)
        numBits = 8 * numBytes
        mask = 0
        data = 0

        for i in xrange(numFields):
            field = fields[i]
            startPos = fieldPos[field]
            if field in self.fields:
                fieldWidth = self.fields[field].width
            else:
                fieldWidth = self.pseudofields[field].width

            lastBit = startPos[0] * 8 + startPos[1] + fieldWidth
            shift = numBits - lastBit

            if isinstance(vals, list):
                newMask = (2 ** fieldWidth) - 1
                newMask <<= shift

                newData = vals[i] << shift
            else:
                newMask = vals[0][i] << shift
                newData = vals[1][i] << shift

            mask |= newMask
            data |= newData

        contentMask = [0] * numBytes
        contentData = [0] * numBytes
        for i in xrange(numBytes):
            contentMask[i] = int((mask >> (8 * (numBytes - i - 1))) & 0xff)
            contentData[i] = int((data >> (8 * (numBytes - i - 1))) & 0xff)

        return (contentMask, contentData)

    def getFieldByteContents(self, fields, fieldVals, fieldBytes, fieldPos):
        """Get a list of field byte contents for a single set of values

        Params:
          fields     - list of fields to get byte contents of
          fieldVals  - list of list of values for corresponding fields
          fieldBytes - list of byte locaitons of the fields
          fieldPos   - map of field name and first byte/bit in *extracted* bytes

        Return:
          [(contentMask, contentData)]
            - contentMask - mask values of content bytes
            - contentData - data values of content bytes
        """
        content = []
        for val in fieldVals:
            content.append(self.getFieldByteContentsSingle(fields, val, fieldBytes, fieldPos))

        return content

    def getFieldByteLocsContents(self, fields, fieldVals):
        """Get a list of bytes locations and a list of contents

        Params:
          fields - list of fields to get bytes of
          fieldVals  - list of list of values for corresponding fields

        Return:
          (fieldBytes, [(contentMask, contentData)]
            - fieldBytes - list of bytes to extract from header
            - contentMask - mask values of content bytes
            - contentData - data values of content bytes
        """
        (fieldBytes, fieldPos, totalWidth) = self.getFieldByteLocs(fields)
        content = self.getFieldByteContents(fields, fieldVals, fieldBytes, fieldPos)

        return (fieldBytes, content)

    def mergeLocs(self, fieldBytes1, fieldBytes2):
        """Merge two sets of extract byte locations

        Params:
          fieldBytes1 - extract byte locations 1
          fieldBytes2 - extract byte locations 2

        Return:
          (mergedFieldBytes, map1, map2)
            - mergedFieldBytes - merged list of bytes to extract from header
            - map1 - mapping of positions from fieldBytes1
            - map2 - mapping of positions from fieldBytes2
        """

        # Create a list of merged field bytes
        mergedFieldBytes = set(fieldBytes1)
        mergedFieldBytes.update(fieldBytes2)
        mergedFieldBytes = list(mergedFieldBytes)
        mergedFieldBytes.sort()
        numBytes = len(mergedFieldBytes)

        # Work out the mapping from the two set of locs/contents
        map1 = []
        map2 = []
        for field in fieldBytes1:
            for i in xrange(numBytes):
                if field == mergedFieldBytes[i]:
                    map1.append(i)
                    break
        for field in fieldBytes2:
            for i in xrange(numBytes):
                if field == mergedFieldBytes[i]:
                    map2.append(i)
                    break

        return (mergedFieldBytes, map1, map2)

    def mergeContents(self, mergedFieldBytes, map1, map2, content1, content2):
        """Merge two sets of content

        Params:
          mergedFieldBytes - list of merged extract locations
          map1             - mapping from content1 to merged fields
          map2             - mapping from content2 to merged fields
          content1         - content bytes to be merged
          content2         - content bytes to be merged

        Return:
          mergedContent - list of merged content mask/data bytes
        """
        numBytes = len(mergedFieldBytes)

        mask = [0] * numBytes
        data = [0] * numBytes

        for i in xrange(len(content1[0])):
            mask[map1[i]] |= content1[0][i]
            data[map1[i]] |= content1[1][i]

        for i in xrange(len(content2[0])):
            mask[map2[i]] |= content2[0][i]
            data[map2[i]] |= content2[1][i]

        return (mask, data)

    def mergeLocsContents(self, locsContent1, locsContent2):
        """Merge two sets of bytes/locations

        Params:
          locsContent1 - (extractBytes, content) tuple
          locsContent2 - (extractBytes, content) tuple

        Return:
          (fieldBytes, [(contentMask, contentData)]
            - fieldBytes - list of bytes to extract from header
            - contentMask - mask values of content bytes
            - contentData - data values of content bytes
        """

        (mergedFieldBytes, map1, map2) = self.mergeLocs(locsContent1[0], locsContent2[0])

        mergedContent = []
        for content1 in locsContent1[1]:
            for content2 in locsContent2[1]:
                mergedContent.append(self.mergeLocs(mergedFieldBytes, map1, map2, content1, content2))

        return (mergedFieldBytes, mergedContent)


    def getLookupBytes_i(self):
        """Internal function to get list of lookup bytes

        Return:
          (lookupBytes, fieldPos, lookupWidth)
            - lookupBytes - list of bytes to lookup from header
            - fieldPos    - dict of field name -> (first byte, first bit in extracted header)
            - lookupWidth - total width of lookup fields
        """
        if not self.nextHeader or isinstance(self.nextHeader, str):
            return ([], {}, 0)

        if isinstance(self.nextHeader, tuple):
            fields = self.nextHeader[0]
            return self.getFieldByteLocs(fields)

    def getLookupBytes(self):
        """Get a list of the lookup bytes

        Return: array of lookup bytes to extract from header
        """
        return self.getLookupBytes_i()[0]

    def getLengthBytes_i(self):
        """Internal function to get list of length bytes

        Return:
          (lengthBytes, fieldPos, lookupWidth)
            - lengthBytes - list of bytes specifying length
            - fieldPos    - dict of field name -> (first byte, first bit in extracted header)
            - lookupWidth - total width of length fields
        """
        lengthFields = self.getLengthFields()

        if len(lengthFields) == 0:
            return ([], {}, 0)

        return self.getFieldByteLocs(lengthFields)

    def getLengthBytes(self):
        """Get a list of the length bytes

        Return: array of length bytes to extract from header
        """
        return self.getLengthBytes_i()[0]

    def getDecisionBytes(self):
        """
        Get a list of all decision (length & lookup) bytes

        Return: array of decision bytes to extract from header
        """
        return self.mergeLocs(self.getLookupBytes(), self.getLengthBytes())

    def getLengthFields(self):
        """Get a list of variables in the length expression"""
        if not self.calcLength:
            return []

        lengthFields = []
        for tok in self.calcLength:
            if isinstance(tok, str) and not opre.match(tok):
                lengthFields.append(tok)

        return lengthFields

    def getLookupFields(self):
        """Get a list of lookup fields"""
        if not self.nextHeader or isinstance(self.nextHeader, str):
            return []

        if isinstance(self.nextHeader, tuple):
            return self.nextHeader[0]

    def getLookupLengthFields(self):
        """Get a merged list of lookup/length fields"""
        fields = self.getLookupFields() + self.getLengthFields()
        for pos in xrange(len(fields)):
            field = fields[pos]
            firstLoc = self.getFieldByteLocs([field])[0][0]
            fields[pos] = (firstLoc, field)
        fields.sort()
        fields = [field for (loc, field) in fields]

        return fields

    def getLengthVarValues(self):
        """Get the list of length variables and valid values"""
        lengthVars = self.getLengthFields()
        validVals = []

        if len(lengthVars) > 0:
            minLen = self.length()[0]
            maxLen = self.maxLength
            if not maxLen:
                maxLen = MAX_HDR_LEN
            vals = [0] * len(lengthVars)
            term = 0
            calcLen = -1
            valMap = {}
            for i in xrange(len(vals)):
                valMap[lengthVars[i]] = 0
            while calcLen < maxLen and term < 255:
                calcLen = self.doLengthCalc(valMap)
                if calcLen >= minLen and calcLen <= maxLen:
                    validVals.append(copy.copy(vals))
                incPos = len(vals) - 1
                while incPos >= 0:
                    if vals[incPos] == 255:
                        vals[incPos] = 0
                        valMap[lengthVars[incPos]] = 0
                        incPos -= 1
                        if incPos == -1:
                            term = 255
                    else:
                        vals[incPos] += 1
                        valMap[lengthVars[incPos]] = vals[incPos]
                        incPos = -1

        return (lengthVars, validVals)

    def getExtractBytes_i(self):
        """Internal function to get list of bytes to extract

        Return:
          (extractBytes, fieldPos, lookupWidth)
            - extractBytes - list of bytes to extract
            - fieldPos    - dict of field name -> (first byte, first bit in extracted header)
            - lookupWidth - total width of length fields
        """
        if len(self.extractFields) == 0:
            return ([], {}, 0)

        return self.getFieldByteLocs(self.extractFields)

    def getExtractBytes(self):
        """Get a list of the bytes to extract

        Return: array of bytes to extract
        """
        return self.getExtractBytes_i()[0]

    def getFieldWidths(self, fields):
        """Get the widths of multiple fields

        Params:
          fields -- list of fields to get widths of

        Return:
          widths -- list of field widths
        """
        widths = []
        for field in fields:
            if field in self.fields:
                fieldWidth = self.fields[field].width
            else:
                fieldWidth = self.pseudofields[field].width
            widths.append(fieldWidth)
        return widths

    def getDecisionCombos(self, depth, desiredNxtHdrName=ANY, desiredLength=ANY):
        """Number of decision byte combinations up until a specified depth"""
        key = (depth, desiredNxtHdrName, desiredLength)
        if key in self._decCombos:
            return self._decCombos[key]

        #print "getDecisionCombos: %s %d %s %s" % (self.name, depth,
        #        desiredNxtHdrName, desiredLength)
        mergedFieldBytes = []
        matches = []

        # Get the set of lengths and associated variables for this header
        (lenIsVariable, lengths, lenFields, lenFieldBytes, lenFieldPos, lenContent) = self.getHeaderLengths()
        lengths = [x / 8 for x in lengths]
        minLength = min(lengths)
        
        if desiredLength < minLength:
            desiredLength = ANY

        # Walk through the combination of lengths and next headers
        # (assume that headers have been merged)
        if self.nextHeader:
            # Case: one of N next headers

            # Get the fields that are used in the lookup and the set of
            # values->next headers that they can take on
            fields = self.nextHeader[0]
            fieldMap = self.nextHeader[1]
            (fieldBytes, fieldPos, totalWidth) = self.getFieldByteLocs(fields) 

            # If we have a variable-length header, calculate a merger of the decision
            # fields with the length fields
            if lenIsVariable:
                (mergedFieldBytes, map1, map2) = self.mergeLocs(fieldBytes, lenFieldBytes)
            else:
                mergedFieldBytes = fieldBytes

            # Walk through all of the decision value combinations
            for (vals, nextHeaderName) in fieldMap:
                if desiredNxtHdrName == ANY or desiredNxtHdrName == nextHeaderName:
                    fieldMatch = self.getFieldByteContentsSingle(fields, vals, fieldBytes, fieldPos)

                    # Walk through all of the possible lengths
                    for i in xrange(len(lengths)):
                        length = lengths[i]
                        if lenIsVariable:
                            match = self.mergeContents(mergedFieldBytes, map1, map2, fieldMatch, lenContent[i])
                        else:
                            match = fieldMatch

                        if desiredLength == ANY or desiredLength == length:
                            matches.append(match)
        else:
            mergedFieldBytes = lenFieldBytes

            # Case: no next header
            for i in xrange(len(lengths)):
                length = lengths[i]
                if lenIsVariable:
                    match = lenContent[i]
                else:
                    match = None
                if desiredLength == ANY or desiredLength == length:
                    matches.append(match)

        matchStrs = set()
        numDecBytes = 0
        for byteNum in mergedFieldBytes:
            if byteNum < depth:
                numDecBytes += 1

        for match in matches:
            matchStr = ""
            for i in xrange(numDecBytes):
                matchStr += "%02x%02x" % (match[0][i], match[1][i])
            matchStrs.add(matchStr)

        #print matchStrs

        numCombos = len(matchStrs)
        if numCombos < 1:
            numCombos = 1

        self._decCombos[key] = numCombos

        return numCombos

    def getDecisionComboBytes(self, depth, desiredNxtHdrName=ANY, desiredLength=ANY):
        """
        Get the various decision combination values up until up until a specified depth"""
        #print "getDecisionComboBytes: %s %d %s %s" % (self.name, depth,
        #        desiredNxtHdrName, desiredLength)
        mergedFieldBytes = []
        matches = []

        # Get the set of lengths and associated variables for this header
        (lenIsVariable, lengths, lenFields, lenFieldBytes, lenFieldPos, lenContent) = self.getHeaderLengths()
        lengths = [x / 8 for x in lengths]
        minLength = min(lengths)
        
        if desiredLength < minLength:
            desiredLength = ANY

        # Walk through the combination of lengths and next headers
        # (assume that headers have been merged)
        if self.nextHeader:
            # Case: one of N next headers

            # Get the fields that are used in the lookup and the set of
            # values->next headers that they can take on
            fields = self.nextHeader[0]
            fieldMap = self.nextHeader[1]
            (fieldBytes, fieldPos, totalWidth) = self.getFieldByteLocs(fields) 

            # If we have a variable-length header, calculate a merger of the decision
            # fields with the length fields
            if lenIsVariable:
                (mergedFieldBytes, map1, map2) = self.mergeLocs(fieldBytes, lenFieldBytes)
            else:
                mergedFieldBytes = fieldBytes

            # Walk through all of the decision value combinations
            for (vals, nextHeaderName) in fieldMap:
                if desiredNxtHdrName == ANY or desiredNxtHdrName == nextHeaderName:
                    fieldMatch = self.getFieldByteContentsSingle(fields, vals, fieldBytes, fieldPos)

                    # Walk through all of the possible lengths
                    for i in xrange(len(lengths)):
                        length = lengths[i]
                        if lenIsVariable:
                            match = self.mergeContents(mergedFieldBytes, map1, map2, fieldMatch, lenContent[i])
                        else:
                            match = fieldMatch

                        if desiredLength == ANY or desiredLength == length:
                            matches.append(match)
        else:
            mergedFieldBytes = lenFieldBytes

            # Case: no next header
            for i in xrange(len(lengths)):
                length = lengths[i]
                if lenIsVariable:
                    match = lenContent[i]
                else:
                    match = None
                if desiredLength == ANY or desiredLength == length:
                    matches.append(match)

        matchStrs = set()
        numDecBytes = 0
        for byteNum in mergedFieldBytes:
            if byteNum < depth:
                numDecBytes += 1

        for match in matches:
            matchStr = ""
            for i in xrange(numDecBytes):
                matchStr += "%02x%02x" % (match[0][i], match[1][i])
            matchStrs.add(matchStr)

        #print matchStrs

        return matchStrs

    def getDecisionComboBytes2(self, depth, desiredNxtHdrName=ANY, desiredLength=ANY):
        """
        Get the various decision combination values up until up until a specified depth"""
        #print "getDecisionComboBytes: %s %d %s %s" % (self.name, depth,
        #        desiredNxtHdrName, desiredLength)
        mergedFieldBytes = []
        matches = []

        # Get the set of lengths and associated variables for this header
        (lenIsVariable, lengths, lenFields, lenFieldBytes, lenFieldPos, lenContent) = self.getHeaderLengths()
        lengths = [x / 8 for x in lengths]
        minLength = min(lengths)
        
        if desiredLength < minLength:
            desiredLength = ANY

        # Walk through the combination of lengths and next headers
        # (assume that headers have been merged)
        if self.nextHeader:
            # Case: one of N next headers

            # Get the fields that are used in the lookup and the set of
            # values->next headers that they can take on
            fields = self.nextHeader[0]
            fieldMap = self.nextHeader[1]
            (fieldBytes, fieldPos, totalWidth) = self.getFieldByteLocs(fields) 

            # If we have a variable-length header, calculate a merger of the decision
            # fields with the length fields
            if lenIsVariable:
                (mergedFieldBytes, map1, map2) = self.mergeLocs(fieldBytes, lenFieldBytes)
            else:
                mergedFieldBytes = fieldBytes

            # Walk through all of the decision value combinations
            for (vals, nextHeaderName) in fieldMap:
                if desiredNxtHdrName == ANY or desiredNxtHdrName == nextHeaderName:
                    fieldMatch = self.getFieldByteContentsSingle(fields, vals, fieldBytes, fieldPos)

                    # Walk through all of the possible lengths
                    for i in xrange(len(lengths)):
                        length = lengths[i]
                        if lenIsVariable:
                            match = self.mergeContents(mergedFieldBytes, map1, map2, fieldMatch, lenContent[i])
                        else:
                            match = fieldMatch

                        if desiredLength == ANY or desiredLength == length:
                            matches.append(match)
        else:
            mergedFieldBytes = lenFieldBytes

            # Case: no next header
            for i in xrange(len(lengths)):
                length = lengths[i]
                if lenIsVariable:
                    match = lenContent[i]
                else:
                    match = None
                if desiredLength == ANY or desiredLength == length:
                    matches.append(match)

        matchStrs = set()
        numDecBytes = 0
        for byteNum in mergedFieldBytes:
            if byteNum < depth:
                numDecBytes += 1

        for match in matches:
            matchStr = ""
            for i in xrange(numDecBytes):
                matchStr += "%02x%02x" % (match[0][i], match[1][i])
            matchStrs.add(matchStr)

        #print matchStrs

        newMatches = []
        for match in matchStrs:
            matchValueByteArray = map(ord, match.decode('hex'))
            newMatches.append((matchValueByteArray[0::2], matchValueByteArray[1::2]))

        return mergedFieldBytes[0:numDecBytes], newMatches
        #return matchStrs

    def getHeaderLengths(self):
        """Get the header length(s) and the set of fields influencing the length"""
        if not self.hdrLengths:
            # Work out the length of the current field -- may be variable!
            (minLen, maxLen, lenIsVariable, lenFields) = self.length()
            if lenIsVariable:
                (lenFields, lenFieldVals) = self.getLengthVarValues()

                # Calculate the set of valid lengths
                lengths = []
                for lenFieldValSet in lenFieldVals:
                    valMap = {}
                    for i in xrange(len(lenFields)):
                        valMap[lenFields[i]] = lenFieldValSet[i]
                    lengths.append(self.doLengthCalc(valMap))

                # Get the list of bytes to extract for length field identification
                (lenFieldBytes, lenFieldPos, lenTotalWidth) = self.getFieldByteLocs(lenFields)
                lenContent = self.getFieldByteContents(lenFields, lenFieldVals, lenFieldBytes, lenFieldPos)
            else:
                (lenFields, lenFieldVals) = ([], [])
                lengths = [minLen]
                (lenFieldBytes, lenFieldPos) = ([], [])
                lenContent = []

            self.hdrLengths =  (lenIsVariable, lengths, lenFields, lenFieldBytes, lenFieldPos, lenContent)

        return self.hdrLengths

    def _greedyMerge(self, fieldWidths, matches):
        """Attempt to merge multiple matches"""
        merged = True
        #print matches
        while len(matches) > 1 and merged:
            merged = False

            # Identify candidates for merging
            candidates = []
            cPos = 0
            cFound = False
            numParts = len(matches[0][0])
            part = 0
            while not cFound and part < numParts:
                bit = fieldWidths[part] - 1
                while not cFound and bit >= 0:
                    #print part, bit
                    for lPos in xrange(len(matches)):
                        for rPos in xrange(lPos + 1, len(matches)):
                            (leftMask, leftMatch) = matches[lPos]
                            (rightMask, rightMatch) = matches[rPos]
                            if leftMask == rightMask:
                                newMask = copy.copy(leftMask)
                                newMask[part] &= (2 ** fieldWidths[part] - 1) ^ (2 ** bit)
                                #print lPos, rPos, part, bit, leftMask, newMask
                                leftMatchNew = [leftMatch[i] & newMask[i] for
                                        i in xrange(numParts)]
                                rightMatchNew = [rightMatch[i] & newMask[i] for
                                        i in xrange(numParts)]
                                #print leftMatch, leftMatchNew, rightMatch, rightMatchNew
                                if leftMatchNew == rightMatchNew:
                                    #if cPos not in candidates:
                                    #    candidates[cPos] = []
                                    candidates.append((lPos, rPos))
                                    cFound = True
                    if not cFound:
                        cPos += 1
                        bit -= 1
                if not cFound:
                    part += 1

            # Attempt to merge the candidates
            #print "Candidates:", candidates
            #print "Matches:", matches
            mergedPositions = set()
            mergedMatches = []
            while len(candidates) > 0:
                (lPos, rPos) = candidates.pop()
                if lPos not in mergedPositions and rPos not in mergedPositions:
                    (leftMask, leftMatch) = matches[lPos]
                    newMask = copy.copy(leftMask)
                    newMask[part] &= (2 ** fieldWidths[part] - 1) ^ (2 ** bit)
                    newMatch = [leftMatch[i] & newMask[i] for i in
                            xrange(numParts)]
                    mergedPositions.add(lPos)
                    mergedPositions.add(rPos)
                    mergedMatches.append((newMask, newMatch))
                    merged = True
            for i in xrange(len(matches)):
                if i not in mergedPositions:
                    mergedMatches.append(matches[i])
            #print "MergedMatches:", mergedMatches

            matches = mergedMatches

        #print ""
        return matches

    def optNextHeader(self):
        """Attempt to reduce the number of next-header entries by merging"""
        if type(self.nextHeader) == tuple:
            fields = self.nextHeader[0]
            allFields = copy.copy(self.fields)
            allFields.update(self.pseudofields)
            fieldWidths = [allFields[field].width for field in fields]
            transitions = self.nextHeader[1]

            # Identify which headers are pointed to by multiple match values
            dsts = {}
            for (match, nxtHdr) in transitions:
                if nxtHdr not in dsts:
                    dsts[nxtHdr] = []
                dsts[nxtHdr].append(match)

            # Attempt to merge headers
            newTransitions = []
            for dst in dsts:
                if len(dsts[dst]) > 1:
                    mergedMatches = self._greedyMerge(fieldWidths, dsts[dst])
                    for match in mergedMatches:
                        newTransitions.append((match, dst))
                else:
                    newTransitions.append((dsts[dst][0], dst))

            # Update the nextHeader field
            self.nextHeader = (fields, newTransitions)


# Basic test code
if __name__ == '__main__':
    hdr = Header('TestHeader')
    print hdr

