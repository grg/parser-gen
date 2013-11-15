#!/usr/bin/env python

"""Reference counter"""

class RefCount():
    """Reference counter"""

    def __init__(self, name, maxVal = -1):
        self.name = name
        self.maxVal = maxVal
        self.refCount = 0

    def __str__(self):
        return "RefCount: [name='%s', max=%d, refCount=%d]" % (self.name, self.maxVal,
                self.refCount)

    def atLimit(self):
        """Is this counter at the limit"""
        return self.maxVal >= 0 and self.refCount >= self.maxVal

    def exceedsLimit(self):
        """Is this counter at the limit"""
        return self.maxVal >= 0 and self.refCount > self.maxVal

    def inc(self):
        self.refCount += 1

    def dec(self):
        self.refCount -= 1

    def setMax(self, maxVal):
        self.maxVal = maxVal

    def __str__(self):
        return "[name=%s, max=%d, count=%d]" % (self.name, self.maxVal, self.refCount)

# Basic test code
if __name__ == '__main__':
    refCount = RefCount('Example')
    print refCount
