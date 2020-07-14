import numpy

def fingerprint(a):
    '''Fingerprint of numpy array'''
    a = numpy.asarray(a)
    return numpy.dot(numpy.cos(numpy.arange(a.size)), a.ravel())
finger = fp = fingerprint
