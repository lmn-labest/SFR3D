#!/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 10:41:53 2016

@author: henrique
"""

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties


def grafico(prename,mesh,ifig):
    sufixo= ['_diag_log_hist1.txt' 
            ,'_diag2_log_hist1.txt'
            ,'_diag4_log_hist1.txt'
            ,'_diag8_log_hist1.txt'
            ,'_diag16_log_hist1.txt'
            ,'_diag32_log_hist1.txt'
            ,'_diag64_log_hist1.txt'
            ]

    legenda = ['Diag'
              ,'Block 2x2'
              ,'Block 4x4'
              ,'Block 8x8'
              ,'Block 16x16'
              ,'Block 32x32'
              ,'Block 64x64']

    x     = []
    y     = []
    names = []
    it    = 0

    for name in sufixo:
        names.append(prename+name)

    for name in names:
        print (name),it
        with open(name,"r") as f:
            data = f.read()


        data = data.split('\n')

        sdata = data[1:len(data)-1]

    
        for row in sdata:
            srow = row.split()    
            x.append(srow[1]) 
            y.append(srow[2]) 

        plt.figure(ifig)
        plt.semilogy(x,y, basey= 10,label=legenda[it])
        plt.xlabel('Iterations')
        plt.ylabel('Log(r)')
        plt.title("CG: "+mesh)
#    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.legend(prop={'size':8})
        plt.grid(True)
    
        it= it + 1
        x.clear()
        y.clear()
 
#plt.show()
        plt.savefig(prename+'.png',ext="png", close=True, verbose=True
                   ,bbox_inches='tight',dpi=1000)    

def main():
    
    prename = 'solo1_mec'
    mesh    = 'Mesh 1'
    grafico(prename,mesh,1)    

    prename = 'solo1_1_mec'
    mesh    = 'Mesh 2'
    grafico(prename,mesh,2)  

    prename = 'solo2_mec'
    mesh    = 'Mesh 3'
    grafico(prename,mesh,3)  

    prename = 'solo2_1_mec'
    mesh    = 'Mesh 4'
    grafico(prename,mesh,4)

    prename = 'solo2_2_mec'
    mesh    = 'Mesh 5'
    grafico(prename,mesh,5)
	
    prename = 'solo2_3_mec'
    mesh    = 'Mesh 6'
    grafico(prename,mesh,6)
	
    prename = 'solo3_mec'
    mesh    = 'Mesh 7'
    grafico(prename,mesh,7)		

if __name__ == "__main__":
    main()
