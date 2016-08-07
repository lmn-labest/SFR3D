#!/usr/bin/python

import sys


def meshConv(name):

       

    prename = name.split('.')[0]
    outName = prename + '_conv.dat'
    fOut    = open(outName,"w")
    fOut.write('cells\n') 
    with open(name,"r") as fIn:
        data = fIn.read()

        data = data.split('\n')
        
#tipo da malha
        cellType = data[0]
        print cellType
        if cellType == 'quad4' :
            nen   = 4
            nFace = 4
        
        sdata = data[1:len(data)-3]
        
        for row in sdata:
          srow = row.split()
# numel
          strV = str(srow[0]) + ' '    
          fOut.write(strV)
# mat
          strV = str(srow[nen+1])     
          fOut.write(strV)
          
# 
          strV =' 3 ' + str(nen) + ' ' + str(nFace) + ' ' 
          fOut.write(strV)

#conectividades
          for no in range(1,nen+1):
            strV = str(srow[no]) + ' '
            fOut.write(strV)
             
          fOut.write('\n') 
            

    fOut.write('endCells\nreturn\n') 
        
    fOut.close()    


def main(argv):
    
    argc = len(argv) 
    if argc < 2:
        print ("%s : input")%(argv[0])
        return    

    name    = argv[1] 

    meshConv(name)
    
if __name__ == "__main__":
    main(sys.argv)
