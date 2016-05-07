# coding: utf-8

import matplotlib as ma
ma.use('Agg')
import matplotlib.pyplot as plt

import math as m
import numpy as np
from splint import wings_interpolation

def height(yUpper,yLower):
    """ The "height" function takes 2 arguments :
        - yUpper: array of real, the list of ordinate for the upper wing
        - yLower: array of real, the list of ordinate for the lower wing
        Returns 2 values:
        - hMin: real, the lowest ordinate in yLower
        - hMax: real, the highest ordinate in yUpper
    """
    hMin=yUpper[0]
    hMax=yLower[0]
    n=len(yUpper)
    for i in range(1,n):
        if yUpper[i]>hMax:
            hMax=yUpper[i]
        if yLower[i]<hMin:
            hMin=yLower[i]
    print(hMin, hMax)
    return hMin,hMax


def f_lambda(yTable,lambd,h):
    """ The "f_lambda" function takes 4 arguments :
        - yTable: array of real, contains the ordinate of the points
        - lambd: real, vertical position of the current curve
        - h: real, 3h is the limit of vertical interval
        Returns f, an array of real of size len(yTable). f is the list of ordinates of one curve.
    """
    f=[]
    for i in range(len(yTable)):
        f+=[ ((1-lambd)*yTable[i]) + (lambd*3*h) ]
    return f


def compute_curves(yTable,lambdEps,h):
    """ The "compute_curves" function takes 4 arguments :
        - yTable: array of real, contains the ordinate of the points
        - lambdEps: real, vertical step between each curve to calculate
        - h: real, 3h is the limit of vertical interval
        Returns curves, an array of 1/lambdEps arrays. Each array represent the list of ordinates of one curve.
    """
    lambd=lambdEps
    curves=[]
    for i in range(int(1./lambdEps)):
        curves+=[f_lambda(yTable,lambd,h)]
        lambd+=lambdEps
    return curves


def display_curves(yUpper,yLower,curves,xEps):
    """ The "compute_curves" function takes 4 arguments :
        - yUpper: array of real, the list of ordinate for the upper wing
        - yLower: array of real, the list of ordinate for the lower wing
        - curves:
        - xEps: real, step for x
        Displays the airflow below and above the wing
    """
    x=[0]
    i=0
    while (x[i]<1):
        x+=[x[i]+xEps]
        i+=1
    x=x[:(len(x)-1)]
    
    plt.ylim(-0.13,0.33)
    plt.title("Laminar flow above and below the wing")
    plt.xlabel("airfoil x-coordinates")
    plt.ylabel("airfoil y-coordinates")
    plt.plot(x,yUpper, linewidth=1.0)
    plt.plot(x,yLower, linewidth=1.0)
    for i in range(len(curves)):
        plt.plot(x,curves[i], linewidth=1.0)
    plt.show()
    #plt.savefig("DU84132V_curves.png")
    
    
    
def main():
    xEps=0.001
    yUpper,yLower=wings_interpolation("DU84132V.DAT",xEps)
    hMin,hMax=height(yUpper,yLower)
    curves=[]
    lambdaEpsUpper=0.1
    lambdaEpsLower=0.1
    curves+=compute_curves(yUpper,lambdaEpsUpper,hMax)
    curves+=compute_curves(yLower,lambdaEpsLower,hMin)
    display_curves(yUpper,yLower,curves,xEps)



if __name__ ==  '__main__':
    main()