import math as m
import numpy as np
import matplotlib.pyplot as plt
from splint import wings_interpolation

def height(yUpper,yLower):
    hMin=yUpper[0]
    hMax=yLower[0]
    n=len(yUpper)
    for i in range(1,n):
        if yUpper[i]>hMax:
            hMax=yUpper[i]
        if yLower[i]<hMin:
            hMin=yLower[i]
    
    return hMin,hMax

def f_lambda(yTable,xEps,lambd,h):
    
    f=[]
    for i in range(len(yTable)):
        f+=[ ((1-lambd)*yTable[i]) + (lambd*3*h) ]

    return f

def compute_curves(y,xEps,lambdEps,h):

    lambd=lambdEps
    curves=[]
    for i in range(int(1./lambdEps)):
        curves+=[f_lambda(y,xEps,lambd,h)]
        lambd+=lambdEps

    return curves

def display_curves(yUpper,yLower,curves,xEps):
    x=[0]
    i=0
    while (x[i]<1):
        x+=[x[i]+xEps]
        i+=1
    x=x[:(len(x)-1)]
    
    plt.ylim(-0.13,0.33)
    plt.plot(x,yUpper, linewidth=1.0)
    plt.plot(x,yLower, linewidth=1.0)
    for i in range(len(curves)):
        plt.plot(x,curves[i], linewidth=1.0)
    plt.show()
    
def main():
    xEps=0.001
    yUpper,yLower=wings_interpolation("DU84132V.DAT",xEps)
    hMin,hMax=height(yUpper,yLower)
    curves=[]
    lambdaEpsUpper=0.1
    lambdaEpsLower=0.1
    curves+=compute_curves(yUpper,xEps,lambdaEpsUpper,hMax)
    curves+=compute_curves(yLower,xEps,lambdaEpsLower,hMin)
    display_curves(yUpper,yLower,curves,xEps)

main()
