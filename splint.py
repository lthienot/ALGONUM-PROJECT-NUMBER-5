# coding: utf-8

import matplotlib as ma
ma.use('Agg')
import matplotlib.pyplot as plt

import math as m
import numpy as np
from spline import data_spline
from load_foil import load_foil

def spline_interpolation(xa,ya,y2a,n,x):
    """ The "spline_interpolation" function takes 5 arguments :
        - xa: 
        - ya: 
        - y2a: 
        - n: 
        - x: 
        Returns ...
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



def splint_aux(n,x,xTable,yTable,scdDerivTable,iOverX,iUnderX):
    """ The "splint_aux" function takes 8 arguments :
        - n: integer, size of xTable, yTable and scdDerivTable (the 3 arrays have the same size)
        - x: real, the abscissa where y is calculated
        - xTable: array of real, contains the abscissa of the points
        - yTable: array of real, contains the ordinate of the points
        - scdDerivTable: array of real, contains the spline corresponding at the xTable and the yTable.
        - iOverX: 
        - iUnderX: 
        Returns 3 values:
        - y: real, the calculated ordinate of the given point x
        - iUnderX
        - iOverX
    """
    a=0.
    b=0.
    DiffOverUnderX=0.
    if (xTable[iOverX]<=x):
        iOverX+=1
        iUnderX+=1
    assert(iOverX<n)
        
    DiffOverUnderX=float(xTable[iOverX])-float(xTable[iUnderX])
    a=(xTable[iOverX]-x)/DiffOverUnderX
    b=(x-xTable[iUnderX])/DiffOverUnderX
    
    return a*yTable[iUnderX]+b*yTable[iOverX]+(((a**3)-a)*scdDerivTable[iOverX])*(DiffOverUnderX**2)/6.0 , iUnderX, iOverX



def wings_interpolation(filename,eps):
    """ The "wings_interpolation" function takes 2 arguments :
        - filename: string, the name of the file to load
        - eps: integer, step for x
        Returns tables of cubic-spline interpolated values y, for x from 0 to 1 (step eps), for the upper and lower part.
    """
    (nUpper,nLower,xTableUpper,xTableLower, yTableUpper, yTableLower, scdDerivTableUpper, scdDerivTableLower) = data_spline(filename)

    yUpper=[] #Final interpolation of the upper part
    yLower=[] #Final interpolation of the lower part
    x=[0]

    ### Upper part:
    iUnderX=0
    iOverX=1
    i=0
    while (x[i]<1):
        (y,iUnderX,iOverX)=splint_aux(nUpper,x[i],xTableUpper,yTableUpper,scdDerivTableUpper,iOverX,iUnderX)
        yUpper+=[y]
        x+=[x[i]+eps]
        i+=1
    x=x[:len(yUpper)]

    ### Lower part:
    iOverX=1
    iUnderX=0
    i=0
    while (i<len(x)):
        (y,iUnderX,iOverX)=splint_aux(nLower,x[i],xTableLower,yTableLower,scdDerivTableLower,iOverX,iUnderX)
        yLower+=[y]
        i+=1
        
    return yUpper,yLower



def display_wing(yUpper,yLower,eps):
    """ The "display_wing" function takes 3 arguments :
        - yUpper: array of real, the list of ordinate for the upper wing
        - yLower: array of real, the list of ordinate for the lower wing
        - eps: integer, step for the 2 arrays
        Displays the figure
    """
    x=[0]
    i=0
    while (x[i]<1):
        x+=[x[i]+eps]
        i+=1
    plt.ylim(-0.5,0.5)
    plt.plot(x[:(len(yUpper))],yUpper, linewidth=1.0)
    plt.plot(x[:(len(yLower))],yLower, linewidth=1.0)
    #plt.show()
    plt.savefig("DU84132V.png")
    
    
    
def main():
    yUpper2,yLower2=wings_interpolation("DU84132V.DAT",0.001)
    display_wing(yUpper2,yLower2,0.001)



if __name__ ==  '__main__':
    main()
    

