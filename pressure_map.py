# coding: utf-8

import matplotlib as ma
ma.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as col
import math as m
import numpy as np
from splint import wings_interpolation
from integration import length
from integration import simpson
from integration import trapezium

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
    """ The "display_curves" function takes 4 arguments :
        - yUpper: array of real, the list of ordinate for the upper wing
        - yLower: array of real, the list of ordinate for the lower wing
        - curves: array of 1/lambdEps arrays. Each array represent the list of ordinates of one curve.
        - xEps: real, step for x
        Displays the airflow below and above the wing
    """
    x=[0]
    i=0
    while (x[i]<1):
        x+=[x[i]+xEps]
        i+=1
    x=x[:(len(x)-1)]
    hMin,hMax=height(yUpper,yLower)
    plt.ylim(3*hMin,3*hMax)
    plt.title("Laminar flow above and below the wing")
    plt.xlabel("airfoil x-coordinates")
    plt.ylabel("airfoil y-coordinates")
    plt.plot(x,yUpper, linewidth=1.0)
    plt.plot(x,yLower, linewidth=1.0)
    for i in range(len(curves)):
        plt.plot(x,curves[i], linewidth=1.0)
    plt.show()
    plt.savefig("DU84132V_curves.png")



def mat_pressure(curves,eps,hmax,hmin):
    """ The "mat_pressure" function takes 4 arguments :
        - curves: array of 1/lambdEps arrays. Each array represent the list of ordinates of one curve.
        - eps: real, step for x
        - hMin: real, the lowest ordinate of the curves
        - hMax: real, the highest ordinate of the curves
        Displays the pressure map
    """
    curves_length=[]
    for i in range(len(curves)):
        curves_length+=[length(curves[i],trapezium,eps)]

    nX=int(1/eps)
    nY=int((hmax-hmin)/(eps/3))
    mat=np.zeros((nY,nX))
    mini=1
    maxi=-1
    for i in range(nX):
        for j in range(nY):
            sum=0
            nbvalues=0
            for k in range(len(curves)-1,0,-1):
                if (curves[k][i] <= hmin+(eps/3)*j and curves[k][i] > hmin+(eps/3)*(j-1)):
                    nbvalues+=1
                    sum+=(curves_length[k])
            if (nbvalues==0):
                nbvalues=1
            mat[j][i]=(sum/nbvalues)**2
            if (mat[j][i]>max):
                maxi=mat[j][i]
            if(mat[j][i]<min):
                mini=mat[j][i]

    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            mat[i][j]=(mat[i][j]-mini)*(256/(maxi-mini))

    c_dict_pressure = [
                '#ffffff', # noir
                '#ff0000', # rouge
                '#00ffff', # jaune
                '#000000'  # blanc
        ]

    c_map_pressure = col.LinearSegmentedColormap.from_list('pressure', c_dict_pressure,  N=256, gamma=1.0)
    cm.register_cmap(cmap=c_map_pressure)

    plt.ylim(0,nY-1)
    plt.xlim(0,nX-1,0.1)
    plt.imshow(mat, interpolation='none', cmap='pressure')
    plt.show
    plt.savefig("mat.png")                   
    
    
    
def main():
    xEps=0.01           # Have to be small enough to keep a good precision with the map pressure. Without computed curves, the map pressure can't create points because they don't exist. 
    yUpper,yLower=wings_interpolation("DU84132V.DAT",xEps)
    hMin,hMax=height(yUpper,yLower)
    curves=[]
    lambdaEpsUpper=0.01
    lambdaEpsLower=0.01
    curves+=compute_curves(yUpper,lambdaEpsUpper,hMax)
    curves+=compute_curves(yLower,lambdaEpsLower,hMin)
    
    display_curves(yUpper,yLower,curves,xEps)
    mat_pressure(curves, xEps, 3*hMax,3*hMin)



if __name__ ==  '__main__':
    main()
