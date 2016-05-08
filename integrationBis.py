from math import *
import splint as sp

def deriv(f, step):
    """ Derives a function f 
    f : a function, represented by a tabular where step is the space between two adjacent abscissas 
    """
    
    fprim=[]
    fprim.append( (f[1] - f[0]) / step ) # First special case : derivative of the function on the first point.

    for i in range (1, len(f)-1):
        fprim.append( (f[i+1] - f[i-1]) / step ) # General case : approximation of the derivative on a point with the adjacent points.

    fprim.append( (f[len(f)-1] - f[len(f)-2]) / step ) # Second special case : derivative of the function on the last point.
    
    return fprim


def length(spline, method, step):
    fprim = deriv(spline, step)

    g = map(lambda x: sqrt( 1 + x**2 ), fprim)

    return method(g, step)


def trapezium(f, step):
# trapezoidal rule, approximation method of an integral

    I = 0 

    for i in range(len(f)-1): 
        I += step * (f[i] + f[i+1]) / 2

    return I

def left_rectangle(f, step):
# left rectangle rule, approximation mehod of an integral

    I = 0

    for i in range(len(f)-1):
        I+= step * f[i]

    return I

def test():
    step = 10E-2
    (spline1, spline2) = sp.wings_interpolation("DU84132V.DAT", step)
    print "Trapezoidal rule :"
    print "upper :" , length(spline1, trapezium, step)
    print "lower :" , length(spline2, trapezium, step)
    print " "
    print "Left rectangle rule :"
    print "upper :" , length(spline1, left_rectangle, step)
    print "lower :" , length(spline2, left_rectangle, step)

test()
