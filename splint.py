import math as m
import numpy as np
from spline import *
from load_foil import load_foil
import matplotlib.pyplot as plt

def spline_interpolation(xa,ya,y2a,n,x):

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

def splint_aux(i,n,x,xTable,yTable,scdDerivTable,iOverX,iUnderX):
    """Splint_aux is called by splint_improved
    """
    a=0.
    b=0.
    DiffOverUnderX=0.
    
    if (xTable[iOverX]<=x[i]):
        iOverX+=1
        iUnderX+=1
    assert(iOverX<n)
        
    DiffOverUnderX=float(xTable[iOverX])-float(xTable[iUnderX])
    a=(xTable[iOverX]-x[i])/DiffOverUnderX
    b=(x[i]-xTable[iUnderX])/DiffOverUnderX
    
    return a*yTable[iUnderX]+b*yTable[iOverX]+(((a**3)-a)*scdDerivTable[iOverX])*(DiffOverUnderX**2)/6.0 , iUnderX, iOverX

def splint_improved(xTable,yTable,nUpper,nLower,eps):
    """Splint_improved returns tables of cubic-spline interpolated values y, for x from 0 to 1 (step eps), for the upper and lower part
    It takes 5 arguments :
    xTable: array the considered points
    yTable: array tabulating the function (yTable_i=f(xTable_i))
    nUpper: number of given values for the upper part
    nLower: number of given values for the lower part
    eps: step for x
    """

    xTableUpper=xTable[:nUpper] #Initial table must be split into two parts: upper and lower part
    xTableLower=xTable[nUpper:]
    yTableUpper=yTable[:nUpper]
    yTableLower=yTable[nUpper:]

    (fstDeriv1, fstDerivn)=derivate(xTable,yTable) #Calculates first derivative at points 1 and n

    scdDerivTableUpper=spline(xTableUpper,yTableUpper,nUpper,fstDeriv1,fstDerivn) #Computes second derivative
    scdDerivTableLower=spline(xTableLower,yTableLower,nLower,fstDeriv1,fstDerivn)
    
    yUpper=[] #Final interpolation of the upper part
    yLower=[] #Final interpolation of the lower part
    x=[0]

    ### Upper part:
    iUnderX=0
    iOverX=1
    i=0
    while (x[i]<1):
        (y,iUnderX,iOverX)=splint_aux(i,nUpper,x,xTableUpper,yTableUpper,scdDerivTableUpper,iOverX,iUnderX)
        yUpper+=[y]
        x+=[x[i]+eps]
        i+=1
    x=x[:len(yUpper)]

    ### Lower part:
    iOverX=1
    iUnderX=0
    i=0
    while (i<len(x)):
        (y,iUnderX,iOverX)=splint_aux(i,nLower,x,xTableLower,yTableLower,scdDerivTableLower,iOverX,iUnderX)
        yLower+=[y]
        i+=1

    return yUpper,yLower

def wings_interpolation(filename,eps):
    #"DU84132V.DAT"
    (nx1,nx2,xTable,yTable) = load_foil(filename)
    nUpper=int(nx1[0])
    nLower=int(nx2[0])
    return splint_improved(xTable,yTable,nUpper,nLower,eps)

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

    yUpper2,yLower2=wings_interpolation("DU84132V.DAT",0.001)
    display_wing(yUpper2,yLower2,0.001)


if __name__ ==  '__main__':
    main()
    

