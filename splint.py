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

def splint_improved(xTable,nUpper,nLower,yTable,eps):
    """ """

    n=len(xTable)
    assert(nLower+nUpper==n)

    xTableUpper=ix[:nUpper]
    assert(len(xTableUpper)==nUpper)
    xTableLower=ix[nUpper:]
    assert(len(xTableLower)==nLower)
    yTableUpper=yx[:nUpper]
    assert(len(yTableUpper)==nUpper)
    yTableLower=yx[nUpper:]
    assert(len(yTableLower)==nLower)

    (fstDeriv1, fstDerivn)=derivate(xTable,yTable)

    scdDerivTableUpper=spline(xTableUpper,yTableUpper,nUpper,fstDeriv1,fstDerivn)
    scdDerivTableLower=spline(xTableLower,yTableLower,nLower,fstDeriv1,fstDerivn)
    
    yUpper=[]
    yLower=[]
    x=[0]

    a=0.
    b=0.
    DiffOverUnderX=0.
    iUnderX=0
    iOverX=1
    i=0
    while (x[i]<1):

        if (xTableUpper[iOverX]<=x[i]):
            iOverX+=1
            iUnderX+=1
        assert(iOverX<nUpper)

        DiffOverUnderX=float(xTableUpper[iOverX])-float(xTableUpper[iUnderX])
        a=(xTableUpper[iOverX]-x[i])/DiffOverUnderX
        b=(x[i]-xTableUpper[iUnderX])/DiffOverUnderX

        yUpper+=[a*yTableUpper[iUnderX]+b*yTableUpper[iOverX]+(((a**3)-a)*scdDerivTableUpper[iOverX])*(DiffOverUnderX**2)/6.0]

        x+=[x[i]+eps]
        i+=1


    x=x[:len(yUpper)]
    
    a=0.
    b=0.
    DiffOverUnderX=0.
    iOverX=1
    iUnderX=0
    i=0
    while (i<len(x)):

        if (xTableLower[iOverX]<=x[i]):
            iOverX+=1
            iUnderX+=1
        assert(iOverX<nLower)
        
        DiffOverUnderX=float(xTableLower[iOverX])-float(xTableLower[iUnderX])
        a=(xTableLower[iOverX]-x[i])/DiffOverUnderX
        b=(x[i]-xTableLower[iUnderX])/DiffOverUnderX

        yLower+=[a*yTableLower[iUnderX]+b*yTableLower[iOverX]+(((a**3)-a)*scdDerivTableLower[iOverX])*(DiffOverUnderX**2)/6.0]

        i+=1

    return yUpper,yLower

def wings_interpolation2(filename,eps):
    (nx1,nx2,xTable,yTable) = load_foil(filename)
    nUpper=int(nx1[0])
    nLower=int(nx2[0])
    return splint_improved(xTable,yTable,nUpper,nLower,eps)

def wings_interpolation(filename,eps):
    #"DU84132V.DAT"
    (ex,ey,ix,iy) = load_foil(filename)    
    (yp1, ypn)=derivate(ix, iy)
    n = int(ex[0])
    
    splin1=spline(ix[:n], iy[:n], n, yp1, ypn)
    splin2=spline(ix[n:], iy[n:], n, yp1, ypn)

    yUpper=[]
    yLower=[]
    x=[0]
    i=0
    
    while (x[i]<1):
        yUpper+=[spline_interpolation(ix[:n],iy[:n],splin1,n,x[i])] 
        yLower+=[spline_interpolation(ix[n:],iy[n:],splin2,n,x[i])]
        x+=[x[i]+eps]
        i+=1

    return yUpper,yLower

def display_wing(yUpper,yLower,eps):
    x=[0]
    i=0
    while (x[i]<1):
        x+=[x[i]+eps]
        i+=1
    plt.ylim(-0.5,0.5)
    plt.plot(x[:(len(yUpper))],yUpper, linewidth=1.0)
    plt.plot(x[:(len(yLower))],yLower, linewidth=1.0)
    plt.show()

def main():
    yUpper,yLower=wings_interpolation("DU84132V.DAT",0.001)
    display_wing(yUpper,yLower,0.001)

    yUpper2,yLower2=wings_interpolation2("DU84132V.DAT",0.001)
    display_wing(yUpper2,yLower2,0.001)

if __name__ ==  '__main__':
    main()
    

