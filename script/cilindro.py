#!/usr/bin/python
import sys
import math as m
#**********************************************************************
def main(argv):
  progArg = (['raio1 raio2 raio3 nTeta nRaio fileOut'])
#checando os agumentos  
  nArgs = len(argv)
  if nArgs < 7:
    sys.stderr.write("Usage: %s "%argv[0])
    for arg in progArg:
      print arg 
    return 1

# ... atribuindo variaveis
  r1      = float(argv[1])
  r2      = float(argv[2])
  r3      = float(argv[3])
  nTeta   = int(argv[4])
  nRaio   = int(argv[5])
  fileOut = argv[6]
# ... 
# ... coordenadas
  fileName = fileOut + '.dat'
  f = open(fileName,'w')
  f.write("coordinates\n")
  nCoor = 1
  dteta = 0.5*m.pi/(nTeta-1) 

  dRaio = (r2-r1)/(nRaio-1) 
# interno   
  teta = 0.0
  y0   = 0.0
  r    = r1  
  for i in range(0,nRaio):
    for j in range(0,nTeta):
      x    = r*m.cos(teta)
      y    = r*m.sin(teta)
      teta = teta + dteta 
      f.write("%10d %15.6f %15.6f\n"%(nCoor,x,y))
      nCoor = nCoor + 1
    teta   = 0.0 
    r      = r + dRaio
  
  dRaio = (r3-r2)/(nRaio-1) 
# externo   
  teta = 0.0
  y0   = 0.0
  r    = r2 + dRaio
  for i in range(0,nRaio-1):
    for j in range(0,nTeta):
      x    = r*m.cos(teta)
      y    = r*m.sin(teta)
      teta = teta + dteta 
      f.write("%10d %15.6f %15.6f\n"%(nCoor,x,y))
      nCoor = nCoor + 1
    teta   = 0.0 
    r      = r + dRaio
  f.write("end coordinates\n")

# interno elementos
  f.write("quad4\n")
  nElm = 1
  for i in range(0,nRaio-1):
    for j in range(0,nTeta-1):
      no1 = 1+j+i*nTeta
      no2 = no1 + 1
      no4 = (i+1)*nTeta + j + 1 
      no3 = no4 + 1 
      f.write("%9d %6d %6d %6d %6d %2d\n"%(nElm,no1,no4,no3,no2,1))
      nElm  = nElm + 1

# externo elementos
  for i in range(0,nRaio-1):
    for j in range(0,nTeta-1):
      no1 = j+ (i+nRaio-1)*nTeta + 1
      no2 = no1 + 1
      no4 = (i+nRaio)*nTeta + j + 1 
      no3 = no4 + 1 
      f.write("%9d %6d %6d %6d %6d %2d\n"%(nElm,no1,no4,no3,no2,2))
      nElm  = nElm + 1
#
  f.write("end quad4\n")
#**********************************************************************

#condicao de contorno termico
  f.write("constraintemp\n")
  for i in range(0,nTeta):
    f.write("%9d %6d\n"%(i+1,1))
  for i in range(nTeta*(2*nRaio-2),nTeta*(2*nRaio-1)):
    f.write("%9d %6d\n"%(i+1,1))
  f.write("end constraintemp\n")
  f.write("nodalsources\n")
  for i in range(0,nTeta):
    f.write("%9d %16.6f\n"%(i+1,600.0))
  for i in range(nTeta*(2*nRaio-2),nTeta*(2*nRaio-1)):
    f.write("%9d %16.6f\n"%(i+1,200.0))
  f.write("end nodalsources\n")
#**********************************************************************

#condicao de contorno termico
  f.write("initialtemp\n")
  for i in range(0,nTeta*(2*nRaio-1)):
    f.write("%9d %16.6f\n"%(i+1,25.0))
  f.write("end initialtemp\n")
#**********************************************************************

#condicao de contorno mecanico
  f.write("constraindisp\n")
  for i in range(0,2*nRaio-1):
    no1 = 1+i*nTeta
    f.write("%9d %2d %2d\n"%(no1,0,1))
  for i in range(1,2*nRaio):
    no1 = i*nTeta 
    f.write("%9d %2d %2d\n"%(no1,1,0))
  f.write("end constraindisp\n")
#**********************************************************************
  print 'nCoor',nCoor-1,'nElm',nElm-1
  f.write("return\n")
  f.close()

#**********************************************************************
if __name__ == '__main__':
  sys.exit(main(sys.argv))
