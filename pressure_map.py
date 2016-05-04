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
    for i in range(len(y)):
        f+=[ ((1-lambd)*yTable[i]) + (lambd*3*h) ]

    return f

def compute_curves(yUpper,yLower,xEps,lambdEps,hMax,hMin):

    lambd=lambdEps
    curves=[]
    for i in range(1./lambdEps):
        curves+=[f_lambda(yUpper,xEps,lambd,hMax)]
        curves+=[f_lambda(yLower,xEps,lambd,hMin)]
        lambd+=lambdEps
        
    return curves

def display_curves(yUpper,yLower,curves,xEps)
    x=[0]
    i=0
    while (x[i]<1):
        x+=[x[i]+xEps]
        i+=1
    x=x[:(len(x)-1)]
    
    plt.ylim(-0.5,0.5)
    plt.plot(x,yUpper, linewidth=1.0)
    plt.plot(x,yLower, linewidth=1.0)
    for i in range(len(curves)):
        plt.plot(x,curves[i], linewidth=1.0)
    plt.show()
    
def main():
    xEps=0.001
    yUpper,yLower=wings_interpolation("DU84132V.DAT",xEps)
    hMin,hMax=height(yUpper,yLower)
    lambdEps=0.1
    curves=compute_curves(yUpper,yLower,xEps,lambdEps,hMax,hMin)
    display_curves()

main()
