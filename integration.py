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


def integration(filename, method, step):
    (f1, f2) = sp.wings_interpolation(filename, step)
    fprim1 = deriv(f1, step)
    fprim2 = deriv(f2, step)

    g1 = map(lambda x: sqrt( 1 + x**2 ), fprim1)
    g2 = map(lambda x: sqrt( 1 + x**2 ), fprim2)

    return method(g1, step) + method(g2, step)


def trapezium(f, step):
# midpoint rule approximation method of an integral

    I = 0 

    for i in range(len(f)-1): 
        I += step * (f[i] + f[i+1]) / 2

    return I

#print integration("DU84132V.DAT", trapezium, 10E-2)
