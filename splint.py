# coding: utf-8

import matplotlib as ma
ma.use('Agg')
import matplotlib.pyplot as plt

import math as m
import numpy as np
from spline import spline_data, calculate_spline_data
from load_foil import load_foil
from scipy.interpolate import interp1d


def splint_aux(n,x,xTable,yTable,scdDerivTable,iOverX):
    """ The "splint_aux" function takes 8 arguments :
        - n: integer, size of xTable, yTable and scdDerivTable (the 3 arrays have the same size)
        - x: real, the abscissa where y is calculated
        - xTable: array of real, contains the abscissa of the points
        - yTable: array of real, contains the ordinate of the points
        - scdDerivTable: array of real, contains the spline corresponding at the xTable and the yTable.
        - iOverX: integer, current indice to access abscissa and ordinate in xTable and yTable over (and under with iOverX-1) the approximated x value.
        Returns 3 values:
        - y: real, the approximated ordinate to the given point x
        - iOverX: integer, current indice updated (then, x is between xTable[iOverX-1] and xTable[iOverX])
    """
    a=0.
    b=0.
    DiffOverUnderX=0.
    if (xTable[iOverX]<=x):
        iOverX+=1
    assert(iOverX<n)
        
    DiffOverUnderX=float(xTable[iOverX])-float(xTable[iOverX-1])
    a=(xTable[iOverX]-x)/DiffOverUnderX
    b=(x-xTable[iOverX-1])/DiffOverUnderX
    
    return a*yTable[iOverX-1]+b*yTable[iOverX]+(((a**3)-a)*scdDerivTable[iOverX])*(DiffOverUnderX**2)/6.0 , iOverX



def splint_improved(xTable, yTable, nUpper, nLower, eps=0.001):
    """ The "splint_improved" function takes 5 arguments :
        - xTable: array of real, the considered points
        - yTable: array of real, tabulating the function (yTable_i=f(xTable_i))
        - nUpper: integer, number of given values for the upper part
        - nLower: integer, number of given values for the lower part
        - eps: real, step for x
         Returns tables of cubic-spline interpolated values y, for x from 0 to 1 (step eps), for the upper and lower part
    """
    (nUpper,nLower,xTableUpper,xTableLower, yTableUpper, yTableLower, scdDerivTableUpper, scdDerivTableLower) = calculate_spline_data(xTable, yTable, nUpper, nLower)

    yUpper=[] #Final interpolation of the upper part
    yLower=[] #Final interpolation of the lower part
    x=[0]
    
    ### Upper part:
    iOverX=1
    i=0
    while (x[i]<1):
        (y,iOverX)=splint_aux(nUpper,x[i],xTableUpper,yTableUpper,scdDerivTableUpper,iOverX)
        yUpper+=[y]
        x+=[x[i]+eps]
        i+=1
    x=x[:len(yUpper)]

    ### Lower part:
    iOverX=1
    i=0
    while (i<len(x)):
        (y,iOverX)=splint_aux(nLower,x[i],xTableLower,yTableLower,scdDerivTableLower,iOverX)
        yLower+=[y]
        i+=1
        
    return yUpper,yLower



def wings_interpolation(filename,eps=0.001):
    """ The "wings_interpolation" function takes 2 arguments :
        - filename: string, the name of the file to load
        - eps: real, step for x
        Returns tables of cubic-spline interpolated values y, for x from 0 to 1 (step eps), for the upper and lower part.
    """
    (nUpper,nLower,ix,iy) = load_foil(filename)
    nUpper=int(nUpper[0])
    nLower=int(nLower[0])
    
    return splint_improved(ix, iy, nUpper, nLower, eps)



def display_wing(yUpper,yLower,eps=0.001):
    """ The "display_wing" function takes 3 arguments :
        - yUpper: array of real, the list of ordinate for the upper wing
        - yLower: array of real, the list of ordinate for the lower wing
        - eps: real, step for the 2 arrays
        Displays the wing
    """
    x=[0]
    i=0
    while (x[i]<1):
        x+=[x[i]+eps]
        i+=1
    plt.ylim(-0.13,0.33)
    plt.title("Refined airfoil")
    plt.xlabel("airfoil x-coordinates")
    plt.ylabel("airfoil y-coordinates")
    plt.plot(x[:(len(yUpper))],yUpper, linewidth=1.0)
    plt.plot(x[:(len(yLower))],yLower, linewidth=1.0)
    plt.show()
    #plt.savefig("DU84132V.png")



def test_wing(filename, eps=0.001):
    """ The "test_wing" function takes 2 arguments :
        - filename: string, the name of the file to load
        - eps: real, step for x
        Checks that every yTable data are consonant with the yUpper and yLower arrays data
    """
    (nx1,nx2,xTable,yTable) = load_foil(filename)
    nUpper=int(nx1[0])
    nLower=int(nx2[0])
    yUpper,yLower = splint_improved(xTable, yTable, nUpper, nLower, eps)

    x=[0]
    i=0
    while (x[i]<1):
        x+=[x[i]+eps]
        i+=1
    x=x[:(len(yUpper))]
    for i in range (nUpper):
        if (xTable[i] in x):
            assert(yTable[i] in yUpper)
    for i in range (nUpper,nLower,1):
        if (xTable[i] in x):
            assert(yTable[i] in yLower)
    print("Test interpolation : PASSED.")
    
    
    
    
def test_wing2(filename):
    """ The "test_wing2" function takes 1 argument :
        - filename: string, the name of the file to load
        Displays the wing calculated with the interpolation of the scipy library.
    """
    (nUpper,nLower,xTableUpper,xTableLower, yTableUpper, yTableLower, scdDerivTableUpper, scdDerivTableLower) = spline_data(filename)

    up = interp1d(xTableUpper, yTableUpper)
    low = interp1d(xTableLower, yTableLower)
    
    plt.plot(xTableUpper, up(xTableUpper), linewidth=1.0)
    plt.plot(xTableLower, low(xTableUpper), linewidth=1.0)
    plt.show()
    #plt.savefig(filename+"_lib.png")
    

    
def main():
    filename="DU84132V.DAT"
    eps=0.001
    yUpper2,yLower2=wings_interpolation(filename, eps)
    display_wing(yUpper2, yLower2, eps)
    test_wing(filename)
    test_wing2(filename)



if __name__ ==  '__main__':
    main()
    