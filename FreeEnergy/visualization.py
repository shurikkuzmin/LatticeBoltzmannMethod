import numpy
import pylab

density = numpy.loadtxt("phase020000.dat")
pylab.imshow(density)
pylab.show()
