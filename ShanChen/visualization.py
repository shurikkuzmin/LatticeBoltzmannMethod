import numpy
import pylab

density = numpy.loadtxt("density30000.dat")
pylab.imshow(density)
pylab.show()
