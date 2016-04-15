import math as m
import numpy as np

def spline_interpolation(xa,ya,y2a,n,x,y):
    """Spline_interpolation returns a cubic-spline interpolated value y (given a value of x)
    It takes 6 arguments :
    xa: array of length n of the considered points
    ya: array of length n tabulating the function (ya_i=f(xa_i))
    y2a: array of length n, output of function spline
    n: integer
    x: a value of x
    y: cubic-spline interpolated value y """

    # int k, khi, klo
    # float a, b,h
    a=0.
    b=0.
    h=0.

    k=0
    klo=1
    khi=n
    
    while (khi-klo >1):
        k=(khi+klo)/2
        if (xa[k] > x):
            khi=k
        else:
            klo=k
    h=float(xa[khi])-float(xa[klo])
    if (h == 0.0):
        raise ("Bad xa input to routine splint")
    a=(xa[khi]-x)/h
    b=(x-xa[klo])/h
    return a*ya[klo]+b*ya[khi]+(((a**3)-a)*y2a[khi])*(h**2)/6.0
            
