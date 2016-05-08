# coding: utf-8

import matplotlib as ma
ma.use('Agg')
import matplotlib.pyplot as plt

from math import *
import splint as sp

def deriv(f, step=0.001):
    """ The "deriv" function takes 2 arguments :
        - f: array of real, the list of ordinate for the f function, with x the abscissa
        - step: real, step for x
        Returns fprim, the derivate of points of the function f
    """
    fprim=[]
    fprim.append( (f[1] - f[0]) / step ) # First special case: derivative of the function on the first point.

    for i in range (1, len(f)-1):
        fprim.append( (f[i+1] - f[i-1]) / step ) # General case: approximation of the derivative on a point with the adjacent points.

    fprim.append( (f[len(f)-1] - f[len(f)-2]) / step ) # Second special case : derivative of the function on the last point.
    
    return fprim


def length(spline, method, step=0.001):
    """ The "length" function takes 3 arguments :
        - spline: array of real, contains the ordinate of the f function
        - method: function to use to compute the length
        - step: real, step for x
        Returns the length of the given part of the plane curve
    """
    fprim = deriv(spline, step)
    g = map(lambda x: sqrt( 1 + x**2 ), fprim)
    return method(g, step)



def left_rect(f, step=0.001):
    """ The "left_rect" function takes 2 arguments :
        - f: array of real, the list of ordinate for the f function
        - step: real, step for x
        Returns I, the approximated integral of f
    """
    I = 0 
    for i in range(len(f)-1): 
        I += step * f[i]
    return I



def right_rect(f, step=0.001):
    """ The "right_rect" function takes 2 arguments :
        - f: array of real, the list of ordinate for the f function
        - step: real, step for x
        Returns I, the approximated integral of f
    """
    I = 0 
    for i in range(len(f)-1): 
        I += step * f[i+1]
    return I



def trapezium(f, step=0.001):
    """ The "trapezium" function takes 2 arguments :
        - f: array of real, the list of ordinate for the f function
        - step: real, step for x
        Returns I, the approximated integral of f
    """
# midpoint rule approximation method of an integral
    I = 0 
    for i in range(len(f)-1): 
        I += step * (f[i] + f[i+1]) / 2
    return I




def simpson(f, step=0.001):
    """ The "simpson" function takes 2 arguments :
        - f: array of real, the list of ordinate for the f function
        - step: real, step for x
        Returns I, the approximated integral of f
    """
    I = 0
    n=len(f)-1

    sumA=0
    for i in range(1, n/2):
        sumA+= f[2*i]
    sumB=0
    for i in range(1, n/2+1):
        sumB+= f[2*i-1]
    
    I= (step/3)*(f[0] + 2*sumA + 4*sumB + f[n])
    return I


def compare_rules(file, epsMin=0.05, epsMax=0.0001):
    """ The "compare_rules" function takes 3 arguments :
        - file: string, the name of the file which contains the data to load
        - epsMin: real, the higher interval to compute lengths with different rules
        - epsMax: real, the lower interval to compute lengths with different rules
        Displays a figure comparing the convergence of the different rules.
    """
    x=[epsMin]

    arrTrapez=[]
    arrSimpson=[]
    arrLRect=[]
    arrRRect=[]
    
    while (epsMin>epsMax):
        (spline1, spline2) = sp.wings_interpolation(file, epsMin)
        arrTrapez+=[length(spline1, trapezium, epsMin)]
        arrSimpson+=[length(spline1, simpson, epsMin)]
        arrLRect+=[length(spline1, left_rect, epsMin)]
        arrRRect+=[length(spline1, right_rect, epsMin)]
        epsMin=epsMin-epsMax
        x+=[epsMin]

    x=x[:(len(x)-1)]
    
    plt.xlim(x[0], x[len(x)-1])
    plt.title("Convergence speed of different rules")

    plt.xlabel("Interval size")
    plt.ylabel("Length")
    plt.plot(x,arrTrapez, linewidth=1.0, label='Trapezium')
    plt.plot(x,arrSimpson, linewidth=1.0, label='Simpson')
    plt.plot(x,arrLRect, linewidth=1.0, label='Left Rectangle')
    plt.plot(x,arrRRect, linewidth=1.0, label='Right Rectangle')
    
    plt.legend(
           scatterpoints=1,
           loc='better')
    
    plt.show()
    #plt.savefig("compare_rules.png")
    
    

def test():
    file = "DU84132V.DAT"
    step = 10E-3        # Step have to be lower than the number of point given by the .dat file in order to keep a good precision (step <= 1/load_foil(file)[0] and step <= 1/load_foil(file)[1])

    (spline1, spline2) = sp.wings_interpolation(file, step)
    print("Trapezium rule:")
    print "- Upper:" , length(spline1, trapezium, step)
    print "- Lower:" , length(spline2, trapezium, step)
    print("Simpson rule:")
    print "- Upper:" , length(spline1, simpson, step)
    print "- Lower:" , length(spline2, simpson, step)
    print("Left rectangle rule:")
    print "- Upper:" , length(spline1, left_rect, step)
    print "- Lower:" , length(spline2, left_rect, step)
    print("Right rectangle rule:")
    print "- Upper:" , length(spline1, right_rect, step)
    print "- Lower:" , length(spline2, right_rect, step)
    
    compare_rules(file)

if __name__ ==  '__main__':
    test()
