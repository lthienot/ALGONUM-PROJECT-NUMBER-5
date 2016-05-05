# coding: utf-8

import numpy as np 
import math as m
from load_foil import load_foil

def spline(x, y, n, yp1, ypn):
    """ The "spline" function takes 5 arguments :
        - x: table of real, size n
        - y: table of real, size n
        - n: integer, size of the table x and y
        - yp1: real
        - ypn: real
        Returns y2, a table of real of size n
    """
    y2 = [None]*n
    i = 2
    p = 0
    qn = 0
    sig = 0
    un = 0
    u = [None]*n
    
    n = n-1 #Arrays are defined in the interval [0,n-1] and not [1,n]
    k = n-1
    
    if (yp1 > 0.99e30):
        y2[0]=0
        u[0]=0
    else:
        y2[0]=-0.5
        u[0]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1)
    for i in range(1, n):
        sig = (x[i]-x[i-1])/(x[i+1]-x[i-1])
        p = sig*y2[i-1]+2.0
        y2[i] = (sig-1.0)/p
        u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1])
        u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p
    if (ypn > 0.99e30):
        qn=0.
        un=0.
    else:
        qn=0.5
        un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]))
    y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0)
    while (k >= 0):
        y2[k]=y2[k]*y2[k+1]+u[k]
        k = k - 1
    return y2


def derivate(ix, iy):
    """ The "derivate" function takes 2 arguments :
        - ix: table of real
        - iy: table of real
        Returns 2 values:
        - a1: real, the derivate of the first point of iy.
        - an: real, the derivate of the last point of iy.
    """
    n = len(ix) - 1
    a1 = (iy[1]-iy[0])/(ix[1]-ix[0])
    an = (iy[n-1]-iy[n])/(ix[n-1]-ix[n])
    return (a1, an)
    

def data_spline(file):
    """ The "data_spline" function takes 1 argument :
        - file: string, the name of the file which contains the data.
        Returns 2 values:
        - splin: table of real, the spline of the top of the the wing
        - splin2: table of real, the spline of the bottom of the wing
    """
    (ex,ey,ix,iy) = load_foil(file)
    ex=int(ex[0])
    ey=int(ey[0])

    (yp1, ypn)=derivate(ix[:ex], iy[:ex])
    #print(yp1, ypn)
    splin=spline(ix[:ex], iy[:ex], ex, yp1, ypn)
    
    (yp1, ypn)=derivate(ix[ex:], iy[ex:])
    #print(yp1, ypn)
    splin2=spline(ix[ex:], iy[ex:], ey, yp1, ypn)

    return splin, splin2



def main():
    print(data_spline("aaa.dat"))


if __name__ ==  '__main__':
    main()
    
