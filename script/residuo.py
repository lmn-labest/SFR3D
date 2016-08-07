#!/usr/bin/python

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import sys


def grafico(name,mesh,ifig):

    legenda = ['rU1','rU2','rU3','rMass']

    it    = []
    u1    = []
    u2    = []
    u3    = []
    rp    = []

    prename = name.split('.')

    with open(name,"r") as f:
        data = f.read()


        data = data.split('\n')

        sdata = data[2:len(data)-1]

    
        for row in sdata:
            srow = row.split()    
            it.append(srow[0]) 
            u1.append(srow[1]) 
            u2.append(srow[2]) 
            u3.append(srow[3]) 
            rp.append(srow[4]) 

        plt.figure(ifig)
        plt.semilogy(it,u1, basey= 10,label=legenda[0])
        plt.semilogy(it,u2, basey= 10,label=legenda[1])
        plt.semilogy(it,u3, basey= 10,label=legenda[2])
        plt.semilogy(it,rp, basey= 10,label=legenda[3])
        plt.xlabel('Iterations')
        plt.ylabel('Log(r)')
        plt.title("Simple: "+mesh)
#    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.legend(prop={'size':12})
        plt.grid(True)
    
#       plt.show()
        plt.savefig(prename[0]+'.png',ext="png", close=True, verbose=True
                   ,bbox_inches='tight',dpi=1000)    

def main(argv):
    
    argc = len(argv) 
    if argc < 2:
        print ("%s : input")%(argv[0])
        return    

    name    = argv[1] 
    mesh    = 'Mesh 1'

    grafico(name,mesh,1)
    
if __name__ == "__main__":
    main(sys.argv)
