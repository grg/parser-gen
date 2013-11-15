#!/usr/bin/env python

"""Field within a header"""

class Field():
    """Single field in a header"""
    
    def __init__(self, name, width = -1, pos = -1):
        self.name = name
        #self.minWidth = width
        #self.maxWidth = width
        self.width = width
        self.pos = pos

    def __str__(self):
        if self.width is not None:
            width = '%d' % self.width
        else:
            width = '*'
        return "Field: [name='%s', width='%s', pos='%d']" % (self.name, width,
                self.pos)

    def setWidth(self, width):
        self.width = width
        #if self.minWidth < 0:
        #    self.minWidth = width
        #if self.maxWidth < 0:
        #    self.maxWidth = width

    #def setMinWidth(self, minWidth):
    #    self.minWidth = minWidth

    #def setMaxWidth(self, maxWidth):
    #    self.maxWidth = maxWidth

# Basic test code
if __name__ == '__main__':
    hdr = Field('TestField')
    print hdr

