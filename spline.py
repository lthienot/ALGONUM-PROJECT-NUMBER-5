# coding: utf-8

import numpy as np 
import math as m
from load_foil import load_foil

def spline(x, y, n, yp1, ypn, y2):
    """ The "spline" function takes 5 arguments :
        x: table of real
        y: table of real
        n: integer
        yp1: real
        ypn: real
        y2: table of real
        Return 
    """
    i = 2
    k = n-1
    p
    qn
    sig
    un
    u = [None]*n
    if (yp1 > 0.99e30):
        y2[1]=0
        u[1]=0
    else:
        y2[1]=-0.5
        u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1)
    for i in range(2, n):
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
    while (k > 0):
        u2[k]=y2[k]*y2[k+1]+u[k]
        k = k - 1
    return

    
        
def main():
    (ex,ey,ix,iy) = load_foil("aaa.dat")
    spline(ix, iy, (len(ix[2]), #yp1, ypn, y2)

if __name__ ==  '__main__':
    main()
    
