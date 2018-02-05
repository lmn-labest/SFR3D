import math as m
#import myPlot as plot
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
"""
 matrizUpdate: atualizao da matriz de coeficientes

 x -> (x,y,z)
 p -> (sigma,rho,beta)
"""
def matrizUpdate(x,dt,p):

    sig = p[0]
    rho = p[1]
    beta = p[2]
    
    r0 = 1 + dt*sig
    r1 = -dt*sig
    r2 = 0.0

    r3 = -dt*rho
    r4 = 1 + dt
    r5 = dt*x[0]

    r6 = -dt*x[1]
    r7 = 0.0
    r8 = 1 + dt*beta


    a = np.array([[r0,r1,r2],[r3,r4,r5],[r6,r7,r8]])

    return a

"""
r = b - k*x

b -> vetor de forcas
k -> matriz de coeficeintes
x -> vetor x
"""

def residuo(b,k,x):

    r = 3*[0.0]
    r[0] = b[0] -(k[0][0]*x[0] + k[0][1]*x[1] + k[0][2]*x[2])
    r[1] = b[1] -(k[1][0]*x[0] + k[1][1]*x[1] + k[1][2]*x[2])
    r[2] = b[2] -(k[2][0]*x[0] + k[2][1]*x[1] + k[2][2]*x[2])


    return r

"""
produto interno x*y

"""
def dot(x,y):

    tmp = 0.e0
    for i in range(0,len(x)):
        tmp += x[i]*y[i]

    return tmp

"""
Metodo implicito Euler-newthon-raphson
"""
def implicito(X0,p,dt,nStep):
    tG = []
    xG = []
    yG = []
    zG = []
    t = 0

    x = [X0[0],X0[1],X0[2]]

    tG.append(t)
    xG.append(x[0])
    yG.append(x[1])
    zG.append(x[2])
    for i in range(0,nStep):        
        b     = [x[0],x[1],x[2]]
        modR0 = m.sqrt(dot(b,b))
        for j in range(0,100):
            k = matrizUpdate(x,dt,p)   
   
            r = residuo(b,k,x)
            modR = m.sqrt(dot(r,r))
            if  modR/modR0 < 1.e-11 :
                break
       
            dx = np.linalg.solve(k, r)
            x[0] += dx[0]
            x[1] += dx[1]
            x[2] += dx[2]
            
        t += dt
        tG.append(t)
        xG.append(x[0])
        yG.append(x[1])
        zG.append(x[2])

    return tG, xG, yG, zG

"""
metodo explicito Euler
"""
def explicito(X0,p,dt,nStep):


    sigma = p[0]
    rho   = p[1]
    beta  = p[2]

    tG = []
    xG = []
    yG = []
    zG = []
    t = 0

    x0 = X0[0]
    y0 = X0[1]
    z0 = X0[2]

    tG.append(t)
    xG.append(x0)
    yG.append(y0)
    zG.append(z0)
    for i in range(0,nStep):

        t += dt
        x = x0 + sigma*dt*(y0-x0)
        y = y0 + dt*(x0*(rho-z0) - y0)
        z = z0 + dt*(x0*y0-beta*z0)

        x0 = x
        y0 = y
        z0 = z
        
        tG.append(t)
        xG.append(x)
        yG.append(y)
        zG.append(z)
        

    return tG, xG, yG, zG


def func(x,p):

    sig  = p[0]
    rho  = p[1]
    beta = p[2]

    f1 = sig*(x[1]-x[0])
    f2 = x[0]*(rho-x[2]) - x[1]
    f3 = x[0]*x[1] - beta*x[2]

    return (f1,f2,f3) 


def updateRk4(x,h,k):

   return (x[0] + h*k[0],x[1] + h*k[1],x[2] + h*k[1])

def yNew(x,k1,k2,k3,k4,h):
   

    for i in range(0,3):
        x[i] += (h/6.0)*(k1[i]+2*k2[i]+2*k3[i]+k4[i])

"""
metodo explicito rungaKutta
"""
def rungaKuttaExplicito(x0,p,dt,nStep):

    tG = []
    xG = []
    yG = []
    zG = []
    t = 0

    tG.append(t)
    xG.append(x0[0])
    yG.append(x0[1])
    zG.append(x0[2])

    x = [x0[0],x0[1],x0[2]]

    for i in range(0,nStep):

        t += dt
            
        k1 = func(x,p)

        y = updateRk4(x,0.5*dt,k1)              
        k2 = func(y,p)
              
        y = updateRk4(x,0.5*dt,k2)
        k3 = func(y,p)
              
        y = updateRk4(x,dt,k3)
        k4 = func(y,p) 
              
        yNew(x,k1,k2,k3,k4,dt)
  
        tG.append(t)
        xG.append(x[0])
        yG.append(x[1])
        zG.append(x[2])
        

    return tG, xG, yG, zG

"""
z = x - y
"""
def dif(x,y):
    z = []
    for i in range(0,len(y)):
        z.append(x[i] - y[i])

    return z

def main():

    sigma = 10.0
    rho   = 28
    beta  = 8/3  
    dt    = 1.e-4
    nStep = 1000000

    p  = (sigma,rho,beta)

#    tG , xG, yG, zG =explicito(X0,p,dt,nStep)
#    tG1 , xG1, yG1, zG1 = implicito(X0,p,dt,nStep)
    X0 = (1000.0,1000.0,1000.0)
    tG1 , xG1, yG1, zG1 = rungaKuttaExplicito(X0,p,dt,nStep)
#   eps = 1.e-8
    eps = 1.e-06
    X0 = (1000.0+eps,1000.0+eps,1000.0+eps)
    tG2 , xG2, yG2, zG2 = rungaKuttaExplicito(X0,p,dt,nStep)
 
    d1 = dif(xG1 ,xG2)
#    d2 = dif(xG1,xG2)
    
    f, (ax1,ax2,ax3) = plt.subplots(3,sharex=True)

#   f, ax1 = plt.subplots()

#    rx = (0.0,50.0)
#    ry = (-40,40)

#     ex = r'$\bar t$'
#    ey = r'$\bar x$'

#    f.xlabel(ex,fontsize="xx-large")
    #ax1.ylabel(ey,fontsize="xx-large")
#    f.title(legenda, fontsize="larger") 

#    plt.xlabel(ex,fontsize="xx-large")
#    plt.ylabel(ey,fontsize="xx-large")
#    plt.title("Lorentz", fontsize="larger")   


    ax1.plot(tG1,xG1,label='x1',color='red')
    ax1.set_ylabel("x(t)")

    ax2.plot(tG1,xG2,label='x2',color='blue')
    ax2.set_ylabel(r'$\hat x(t)$')
    
    ax3.plot(tG1,d1,label='x2-x1',color='green')
    ax3.set_xlabel(r'$\bar t$')
    ax3.set_ylabel(r'$\hat x(t)\ - x(t)$')

    file = 'lorentz.png'
    plt.savefig(file,ext="png", close=True, verbose=True
               ,bbox_inches='tight',dpi=100) 
    plt.show()

    print (xG1[-1])

    fig2 = plt.figure()
    ax5 = fig2.gca(projection='3d')
    ax5.grid(True)
    ax5.plot(xG1,yG1,zG1,label="x(t)"         ,color='red')
    ax5.plot(xG2,yG2,zG2,linestyle=':',label=r'$\hat x(t)$' ,color='blue')
    ax5.legend()
    file = 'lorentz_t.png'
    plt.savefig(file,ext="png", close=True, verbose=True
               ,bbox_inches='tight',dpi=100) 
    plt.show()


main()
