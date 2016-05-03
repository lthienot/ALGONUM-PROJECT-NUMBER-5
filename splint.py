import math as m
import numpy as np
from spline import *
from load_foil import load_foil
import matplotlib.pyplot as plt

def spline_interpolation(xa,ya,y2a,n,x):
    """Spline_interpolation returns a cubic-spline interpolated value y (given a value of x)
    It takes 6 arguments :
    xa: array of length n of the considered points
    ya: array of length n tabulating the function (ya_i=f(xa_i))
    y2a: array of length n, output of function spline
    n: integer
    x: a value of x
 """

    # int k, khi, klo
    # float a, b,h
    a=0.
    b=0.
    h=0.

    k=0
    klo=0
    khi=n-1
    
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
    assert(khi<n)
    assert(klo<n)
    return a*ya[klo]+b*ya[khi]+(((a**3)-a)*y2a[khi])*(h**2)/6.0

def wings_interpolation(filename,eps):
    #"DU84132V.DAT"
    (ex,ey,ix,iy) = load_foil(filename)    
    (yp1, ypn)=derivate(ix, iy)
    n = int(ex[0])
    
    splin1=spline(ix[:n], iy[:n], n, yp1, ypn)
    splin2=spline(ix[n:], iy[n:], n, yp1, ypn)

    y1=[]
    y2=[]
    x=[0]
    i=0
    
    while (x[i]<1):
        y1+=[spline_interpolation(ix[:n],iy[:n],splin1,n,x[i])] 
        y2+=[spline_interpolation(ix[n:],iy[n:],splin2,n,x[i])]
        x+=[x[i]+eps]
        i+=1

    return y1,y2

#def display_wing(

def main():
    yover,yunder=wings_interpolation("DU84132V.DAT",0.001)

    # plt.ylim(-0.5,0.5)
    # plt.plot(x,yover, linewidth=1.0)
    # plt.plot(x,yunder, linewidth=1.0)
    # plt.show()


if __name__ ==  '__main__':
    main()
    

