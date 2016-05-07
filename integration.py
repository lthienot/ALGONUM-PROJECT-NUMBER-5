from math import *
import splint as sp

def deriv(f, step=0.001):
    """ The "deriv" function takes 2 arguments :
        - f: array of real, step is the space between two adjacent abscissas
        - step: real, step for x
        Returns fprim, the derivate of the function f ... <- TO CHECK
    """
    fprim=[]
    fprim.append( (f[1] - f[0]) / step ) # First special case: derivative of the function on the first point.

    for i in range (1, len(f)-1):
        fprim.append( (f[i+1] - f[i-1]) / step ) # General case: approximation of the derivative on a point with the adjacent points.

    fprim.append( (f[len(f)-1] - f[len(f)-2]) / step ) # Second special case : derivative of the function on the last point.
    
    return fprim


def length(spline, method, step=0.001):
    """ The "length" function takes 3 arguments :
        - spline: array of real, contains the ordinate of the points of the given part
        - method: function to use ... <- TO CHECK
        - step: real, step for x
        Returns ... <- TO CHECK
    """
    fprim = deriv(spline, step)
    g = map(lambda x: sqrt( 1 + x**2 ), fprim)
    return method(g, step)



def trapezium(f, step=0.001):
    """ The "trapezium" function takes 2 arguments :
        - f: ... <- TO CHECK
        - step: real, step for x
        Returns I, ... <- TO CHECK
    """
# midpoint rule approximation method of an integral
    I = 0 
    for i in range(len(f)-1): 
        I += step * (f[i] + f[i+1]) / 2
    return I



def test():
    step = 10E-2
    (spline1, spline2) = sp.wings_interpolation("DU84132V.DAT", step)
    print "Upper:" , length(spline1, trapezium, step)
    print "Lower:" , length(spline2, trapezium, step)



if __name__ ==  '__main__':
    test()