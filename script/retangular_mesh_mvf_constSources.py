#!/usr/bin/python
import sys
import math as m

#**********************************************************************
def main(argv):
  progArg = (['largura altura nl nh fileOut'])
#checando os agumentos  
  nArgs = len(argv)
  if nArgs < 4:
    sys.stderr.write("Usage: %s "%argv[0])
    for arg in progArg:
      print arg 
    return 1

# ... atribuindo variaveis
  l       = float(argv[1])
  h       = float(argv[2])
  nl      = int(argv[3])
  nh      = int(argv[4])
  fileOut = argv[5]
# ... 
  dl = l/float(nl-1)
  dh = h/float(nh-1)
  x0 = 0.0
  y0 = 0.0
  print dl,dh
  print nl,nh
  print 'nCoor: ',nl*nh,'nEl: ',(nl-1)*(nh-1)
# ... coordenadas
  fileName =fileOut +'.dat'
#  fileName = fileOut + '_coor.dat'
  f = open(fileName,'w')
  f.write("coordinates\n")
  nCoor = 1
  y = y0
  for i in range(0,nh):
    for j in range(0,nl):
      x = x0 + float(j*dl)
      f.write("%9d %16.6f %16.6f\n"%(nCoor,x,y))
      nCoor = nCoor + 1
    y = y + float(dh) 
  f.write("endCoordinates\n")
#  f.close()

# ... elementos
#  fileName = fileOut + '_elem.dat'
#  f = open(fileName,'w')
  f.write("cells\n")
  nElm = 1
  for i in range(0,nh-1):
    for j in range(0,nl-1):
      no1 = 1+j+i*nl
      no2 = no1 + 1
      no4 = (i+1)*nl + j + 1 
      no3 = no4 + 1 
      f.write("%9d %2d %2d %2d %2d %6d %6d %6d %6d\n"%(nElm,1,3,4,4,no1,no2,no3,no4))
      nElm  = nElm + 1
  f.write("endCells\n")
#  f.close()
# ... condicoe de contorno
#  fileName = fileOut + '_res.dat'
#  f = open(fileName,'w')
  f.write("faceRd1\n")
#parede inferior

#parede superior
  nElm = 1
  for i in range(0,nh-1):
    for j in range(0,nl-1):
      if nElm == i*(nl-1)+1 :
        f.write("%9d %2d %2d %2d %2d %2d %2d\n"%(nElm,5,0,0,0,1,2))
      else:
        f.write("%9d %2d %2d %2d %2d %2d %2d\n"%(nElm,5,0,0,0,0,2))
      nElm  = nElm + 1
  f.write("endFaceRd1\n")
  
  f.write("faceSd1\n")
  nElm = 1
  for i in range(0,nh-1):
    for j in range(0,nl-1):
      f.write("%9d %2d %2d %2d %2d %2d %2d\n"%(nElm,5,0,0,0,0,1))
      nElm  = nElm + 1

  f.write("endFaceSd1\n")
#
# f.write("elveloc\n")
# for i in range(1,(nl-1)*(nh-1)):
#   f.write("%9d %.2f %.2f\n"%(i,2.0,0))
# f.write("end elveloc\n")
  f.write("endMesh")
  f.close()
      
#**********************************************************************

if __name__ == '__main__':
  sys.exit(main(sys.argv))
