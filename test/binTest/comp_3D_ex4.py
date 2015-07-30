#!/bin/python
import sys 
import math as m

def main(argv):
  
  progArg = (['fileIn csv'])

  nArgs = len(argv)
  if nArgs < 2:
    sys.stderr.write("Usage: %s "%argv[0])
    for arg in progArg:
      print arg 
    return 1
  
  T1 = 100.0
  Ta = 10.0
  h  = 25.0
  L  = 1.0
  k  = 2.0

  ab =(h*(Ta-T1)/(k+h*L))
  

  fileMvf = open(sys.argv[1],"r")

  mvf     =[]
  gmvf    =[[]]
  x       =[[]]
  er      =[]
  gxer     =[]
  gyer     =[]
  gzer     =[]
  exato   =[]
  gxexato =[]
  gyexato =[]
  gzexato =[]
 
  fileMvf.readline()
  nCell = 0

  for f in fileMvf:
    s = f.split(',')
#...
    x.append([])
    x[nCell].append(float(s[1]))
    x[nCell].append(float(s[2]))
    x[nCell].append(float(s[3]))
#...
    mvf.append(float(s[4]))
#...
    gmvf.append([])
    gmvf[nCell].append(float(s[5]))
    gmvf[nCell].append(float(s[6]))
    gmvf[nCell].append(float(s[7]))
#...
    nCell+=1


  for i in range(0,nCell):
    x1 = x[i][0] 
    x2 = x[i][1]
    x3 = x[i][2]
# ... solucao
    t = ab*(x3+0.5e0) + T1
    exato.append(t)

# ... derivada
    gx =0.0
    gy =0.0
    gz =ab 

    gxexato.append(gx)
    gyexato.append(gy)
    gzexato.append(gz)
    



  erC   = 0.e0      
  gxerC = 0.e0      
  gyerC = 0.e0      
  gzerC = 0.e0      
  for i in range(0,nCell):
# ... solucao
    er.append((exato[i] - mvf[i])**2)
    erC += er[i]

# ... derivada
    gxer.append((gxexato[i] - gmvf[i][0])**2)
    gyer.append((gyexato[i] - gmvf[i][1])**2)
    gzer.append((gzexato[i] - gmvf[i][2])**2)
    gxerC += gxer[i]
    gyerC += gyer[i]
    gzerC += gzer[i]

  exatoMax   = max(exato)     
  gxexatoMax = max(gxexato)     
  gyexatoMax = max(gyexato)     
  gzexatoMax = max(gzexato)     

  erMax   = er[0]
  gxerMax = gxer[0]
  gyerMax = gyer[0]
  gzerMax = gzer[0]
  n=0
  nx=0
  ny=0
  nz=0
  for i in range(1,nCell):
    if erMax <= er[i]:
      erMax = er[i]
      n = i 
    if gxerMax <= gxer[i]:
      xgerMax = gxer[i]
      nx = i 
    if gyerMax <= gyer[i]:
      gyerMax = gyer[i]
      ny = i 
    if gzerMax <= gzer[i]:
      gzerMax = gzer[i]
      nz = i 

  erC   = erC/exatoMax**2 
  gxerC = gxerC/(gxexatoMax**2 + gyexatoMax**2 + gzexatoMax**2)
  gyerC = gyerC/(gxexatoMax**2 + gyexatoMax**2 + gzexatoMax**2)
  gzerC = gzerC/(gxexatoMax**2 + gyexatoMax**2 + gzexatoMax**2)
  print "Solucao:"
  print "erro acumulado normalizado %e" %(erC)  
  print "erro Maximo     er[%d]=%e " %(n+1,erMax) 
  
  print "Derivada:" 
  print "x: erro acumulado normalizado %e" %(gxerC)  
  print "y: erro acumulado normalizado %e" %(gyerC)  
  print "z: erro acumulado normalizado %e" %(gzerC)  
  print "x: erro Maximo     er[%d]=%e " %(nx+1,gxerMax) 
  print "y: erro Maximo     er[%d]=%e " %(ny+1,gyerMax) 
  print "z: erro Maximo     er[%d]=%e " %(nz+1,gzerMax) 
  
  if erC  > 5.0e-1 :
    return 1 
  
# if gxerC  > 1.e-6 :
#   return 1 
 
# if gyerC  > 5.e-6 :
#   return 1 
  
# if gzerC  > 5.e-6 :
#   return 1 

if __name__ == '__main__':
  sys.exit(main(sys.argv))
