#!/bin/python3.6

import Combustao as comb
import readFile  as rf      
import sys
import matplotlib.pyplot as plt 
import math as m

class Mistura(object):

  def __init__(self, reac, yFuel, yAir,fEst=True,density=1.0):
    self.reac  = reac
# ... concentracao estequimetroca
    if fEst:
      yFuel, yAir = self.arFuelEstiometrico()
      
    yO2, yN2   = self.fracOxigenAndNitrogen(yAir)

    self.indexSp = {'Fuel':0,'O2':1,'N2':2,'CO2':3,'H2O':4}
    
    self.yFrac = [yFuel, yO2, yN2, 0.0, 0.0]

    self.density = density

  def __str__(self):
    y = self.yFrac
    nF = self.reac.nameFuel()
    str1  ='{0:<6} : {1:.14f}\n'.format(nF   , y[0])
    str1 +='{0:<6} : {1:.14f}\n'.format('O2' , y[1])
    str1 +='{0:<6} : {1:.14f}\n'.format('N2' , y[2])
    str1 +='{0:<6} : {1:.14f}\n'.format('CO2', y[3])
    str1 +='{0:<6} : {1:.14f}\n'.format('H2O', y[4])

    return str1

  def arFuelEstiometrico(self):

    sO2,sN2,_,_,_ = self.reac.massRel()
    
    pO2  = self.reac.arComp[0]
    pN2  = self.reac.arComp[1]   
    mWO2 = self.reac.species['O2'].molarMass
    mWN2 = self.reac.species['N2'].molarMass

    tmp = 1.0 + sO2 + sN2

    yFuel = 1.0/tmp
    return yFuel,1-yFuel

  def airMolarMass(self):
    pO2  = self.reac.arComp[0]
    pN2  = self.reac.arComp[1]
    spO2 = self.reac.species['O2']
    spN2 = self.reac.species['N2']

    return pO2*spO2.molarMass + pN2*spN2.molarMass

  def fracOxigenAndNitrogen(self, yAir):

    pO2  = self.reac.arComp[0]
    pN2  = self.reac.arComp[1]

    mWAir = self.airMolarMass()
    mWO2  = self.reac.species['O2'].molarMass
    mWN2  = self.reac.species['N2'].molarMass
            

    return yAir*pO2*mWO2/mWAir, yAir*pN2*mWN2/mWAir

  def molarMassMix(self):

    mFuel = self.reac.species['Fuel'].molarMass
    mO2   = self.reac.species['O2'].molarMass
    mN2   = self.reac.species['N2'].molarMass
    mCO2  = self.reac.species['CO2'].molarMass
    mH2O  = self.reac.species['H2O'].molarMass

    yFrac     = self.yFrac
    molarMass = [mFuel, mO2, mN2, mCO2, mH2O]

    tmp = 0.0
    for y,mM in zip(yFrac, molarMass):
      tmp += y*mM                  

    return tmp

  def getyFrac(self):
    return self.yFrac

  def comb(self,istep=10,dt=1.0,tMix=1.0,eps = 1.0e-3,maxIt = 500,fTime='Rk4'):

# ...
    def updateRk4(x,h,k):

      if isinstance(x,(float,int)):
        xNew = x + h*k
        
      else:
        n = len(x)
        xNew = n*[0.0]

        for i in range(0,n):
          xNew[i] = x[i] + h*k[i]  

      return xNew
# ......................................................................    

#...
    def yUpdateRk4(x,k1,k2,k3,k4,h):
      
      if isinstance(x,(float,int)):
        xNew = x + (h/6.0)*(k1+2*k2+2*k3+k4)  
      else:
        n = len(x)
        xNew = n*[0.0]
        for i in range(0,n):
          xNew[i] = x[i] + (h/6.0)*(k1[i]+2*k2[i]+2*k3[i]+k4[i])

      return xNew
# ......................................................................

# ...
    def func(w,density):

      if isinstance(w,(float,int)):
        k = w/density
      else:
        n = len(w)
        k = n*[0.0]
        for i in range(0,len(w)):
          k[i] = w[i]/density

      return k
# ......................................................................    

# ...   
    def euler(q,y,dt):

      if isinstance(y,(float,int)):
        y1 = y + dt*q
      else:
        n = len(y)
        y1 = n*[0.0]
        for j in range(0,n):
          y1[j] = y[j] + dt*q[j]

      return y1
# ......................................................................

        
# ... massa relativa massa do combustivel
    sO2,_,sCO2,sH2O,_ = self.reac.massRel()
# ......................................................................

# ... abrindo arquivos
    f1 = open('data_frac.out','w')
    f2 = open('data_w.out','w')
    f3 = open('data_hs.out','w')
    f4 = open('data_q.out','w')
    f5 = open('data_t.out','w')
    f6 = open('data_conc.out','w')
# .......................................................................

# ...
    densityMix = self.density
    yFrac0 = self.yFrac
    yFraci = self.yFrac
    w1     = len(self.yFrac)*[0.0]
    w2     = len(self.yFrac)*[0.0]
    w3     = len(self.yFrac)*[0.0]
    w4     = len(self.yFrac)*[0.0]
    tempRef= self.reac.T0
    temp   = self.reac.T1
    hs     = self.entalpia(temp)
    t      = 0.0  
# .......................................................................

# ... escrevendo os resultados
    Mistura.fileRes(f1,0,t,self.yFrac)
    Mistura.fileRes(f3,0,t,[hs])
    Mistura.fileRes(f5,0,t,[temp])
    Mistura.fileRes(f6,0,t,self.concMixture())
# ......................................................................

# ...
    gap = 0
    for i in range(1,istep+1):
# ...
      hs0    = hs
      yFrac0 = self.yFrac
      temp0  = temp
# ...
      hsi    = hs0
      yFraci = yFrac0
      tempi  = temp0
# ... 
      for j in range(0,maxIt):
        
# ...
        if fTime == 'Rk4':
# ... equacao das fracoes massicas         
          Mistura.taxaConsumo(w1,temp0,densityMix,self.reac,yFraci\
                             ,sO2,sCO2,sH2O,tMix,fReac='Arrhenius')
          k1 = func(w1,densityMix)

          yFracRk = updateRk4(yFrac0,0.5*dt,k1)

          Mistura.taxaConsumo(w2,temp0,densityMix,self.reac,yFracRk\
                           ,sO2,sCO2,sH2O,tMix,fReac='Arrhenius')

          k2 = func(w2,densityMix)
          yFracRk = updateRk4(yFrac0,0.5*dt,k2)

          Mistura.taxaConsumo(w3,temp0,densityMix,self.reac,yFracRk\
                           ,sO2,sCO2,sH2O,tMix,fReac='Arrhenius')
          k3 = func(w3,densityMix)
          yFracRk = updateRk4(yFrac0,dt,k3)

          Mistura.taxaConsumo(w4,temp0,densityMix,self.reac,yFracRk\
                             ,sO2,sCO2,sH2O,tMix,fReac='Arrhenius')
          k4 = func(w4,densityMix)

          yFracii = yUpdateRk4(yFrac0,k1,k2,k3,k4,dt)

# ... equacao da entalpia sensivel
          q1 = Mistura.taxaLibEnergia(w1,self.reac,tempRef,tempi)
          k1 = func(q1,densityMix)
          hsRk = updateRk4(hs0,0.5*dt,k1)

          q2 = Mistura.taxaLibEnergia(w2,self.reac,tempRef,tempi)
          k2 = func(q2,densityMix)
          hsRk = updateRk4(hs0,0.5*dt,k2)

          q3 = Mistura.taxaLibEnergia(w3,self.reac,tempRef,tempi)
          k3 = func(q3,densityMix)
          hsRk = updateRk4(hs0,dt,k3)

          q4 = Mistura.taxaLibEnergia(w4,self.reac,tempRef,tempi)
          k4 = func(q4,densityMix)

          hsii = yUpdateRk4(hs0,k1,k2,k3,k4,dt)
# ......................................................................

          
        elif fTime == 'Euler':
# ... equacao das fracoes massicas             
          Mistura.taxaConsumo(w1,tempi,densityMix,self.reac,yFraci\
                           ,sO2,sCO2,sH2O,tMix,fReac='Arrhenius')
          k1  = func(w1,densityMix)
          yFracii = euler(k1,yFrac0,dt)
# ......................................................................

# ... equacao da entalpia sensivel
          q1 = Mistura.taxaLibEnergia(w1,self.reac,tempRef,tempi)
          k1 = func(q1,densityMix)

          hsii = euler(k1,hs0,dt)
# ......................................................................


# ... transformando a entalpia sensivel em temperatura
        tempii = self.entalpyForTemp(hsii,tempi)        
# ......................................................................

# ...
        dyFrac = 0.0
        
        for x1,x0 in zip(yFracii,yFraci):
          dyFrac += (x1-x0)**2

        dyFrac = m.sqrt(dyFrac)

        dhs   = abs(hsii - hsi)
        dtemp = abs(tempii - tempi)
        if dhs < eps and dyFrac < eps:
#          print ('step = {0:9} it = {1:3d} dhs = {2:e} dyFrac = {3:e} dtemp = {4:e}'\
#               .format(i,j,dhs,dyFrac,dtemp))
          break
# ......................................................................
      
# ...
        hsi    = hsii
        yFraci = yFracii
        tempi  = tempii
# .....................................................................
        
# ...
      self.yFrac = yFracii
      hs         = hsii
      temp       = tempii
# .....................................................................

# ...
      t += dt
# ......................................................................      
      
#     print('Time: {0:.3f} {1:.3f}'.format(t,temp))

# ...
      gap+=1
      if gap == 100:
        print ('step = {0:9} t = {1:.8f} temp = {2:.6f}'.format(i,t,temp))
        gap = 0  
# ......................................................................

# ... escrevendo os resultados
      Mistura.fileRes(f1,i,t,self.yFrac)
      Mistura.fileRes(f2,i,t,w1)
      Mistura.fileRes(f3,i,t,[hs])
      Mistura.fileRes(f4,i,t,[q1])
      Mistura.fileRes(f5,i,t,[temp])
      Mistura.fileRes(f6,i,t,self.concMixture())
# ......................................................................
      
# ... fechando os arquivos
    f1.close()
    f2.close()
    f3.close()
    f4.close()
    f5.close()
# ......................................................................

    return t,temp,hs

  def concMixture(self):
    x0 = self.concSpecies('Fuel')
    x1 = self.concSpecies('O2')
    x2 = self.concSpecies('N2')
    x3 = self.concSpecies('CO2')
    x4 = self.concSpecies('H2O')

    return [x0,x1,x2,x3,x4]
    

  def entalpia(self,t,t0=298.15):

    yFrac = self.yFrac

    sp0 = self.reac.species['Fuel']
    sp1 = self.reac.species['O2']
    sp2 = self.reac.species['N2']
    sp3 = self.reac.species['CO2']
    sp4 = self.reac.species['H2O']

    cp0 = sp0.cp.integral(t0,t)
    cp1 = sp1.cp.integral(t0,t)
    cp2 = sp2.cp.integral(t0,t)
    cp3 = sp3.cp.integral(t0,t)
    cp4 = sp4.cp.integral(t0,t)

    tmp = yFrac[0]*cp0+yFrac[1]*cp1+yFrac[2]*cp2+yFrac[3]*cp3+yFrac[4]*cp4


    return tmp

  def entalpiaK(sp,t,t0=298.15):

    dhs = sp.cp.integral(t0,t)
    hf  = sp.entalpyFormation()/sp.molarMass
    
    return dhs + hf


  def cpMix(self,t):

    yFrac = self.yFrac

    sp0 = self.reac.species['Fuel']
    sp1 = self.reac.species['O2']
    sp2 = self.reac.species['N2']
    sp3 = self.reac.species['CO2']
    sp4 = self.reac.species['H2O']

    cp0 = sp0.cp.value(t)
    cp1 = sp1.cp.value(t)
    cp2 = sp2.cp.value(t)
    cp3 = sp3.cp.value(t)
    cp4 = sp4.cp.value(t)

    v = yFrac[0]*cp0+yFrac[1]*cp1+yFrac[2]*cp2+yFrac[3]*cp3+yFrac[4]*cp4 
 
    return v

  def entalpyForTemp(self,hs,t):

    ti = t
    for i in range(0,100):
      f  = hs - self.entalpia(ti)
      df = -self.cpMix(ti)
      ti-= f/df
      if abs(f/df) < 1.e-8:        
        break
    
    return ti

  def stoichiometryAirFuel(self):

    pO2  = self.reac.arComp[0]
    nO2, _, _, _, _ = self.reac.stoichCoeff

# ... numero de mols de ar par um mol de O2
    nAir = 1.0/pO2

    mWair = self.airMolarMass()

    mAir  = nO2*nAir*mWair

    mFuel = self.reac.species['Fuel'].molarMass

    
    afSt = mAir/mFuel

    return afSt

  def ratioEquivalente(self):
    afSt = self.stoichiometryAirFuel()

    yFuel = self.yFrac[0]
    yO2   = self.yFrac[1]
    yN2   = self.yFrac[2]
    
    af = (yO2 + yN2)/yFuel

    return afSt/af
  
# ......................................................................

# ...
  def concSpecies(self,name):

   density = self.density 
   mW = self.reac.species[name].molarMass
   y  = self.yFrac[self.indexSp[name]]

   return density*y/mW
# ......................................................................

# ...
  @staticmethod
  def taxaConsumo(w,t,densityMix,reac,yFrac,sO2,sCO2,sH2O,tMix,fReac='edu'):
# ...
    if fReac == 'edu':
# ...
      def edu(densityMix,yFuel,yO2,s=1.0,tMix=1.0):
        return densityMix*min(yFuel,yO2/sO2)/tMix
# ......................................................................

# ...
      omega = edu(densityMix,yFrac[0],yFrac[1],sO2,tMix)

      w[0] = -omega
      w[1] = -sO2*omega
      w[2] =  0.0
      w[3] =  sCO2*omega
      w[4] =  sH2O*omega
# ......................................................................

# ...
    elif fReac == 'Arrhenius':
# ...
      def oneStep(density,yFuel,yO2,mWfuel,mWO2,t):        
        k1f = 1.30e+09*m.exp(-24358.0/t)
        cCh4 = density*yFuel/mWfuel
        cO2  = density*yO2/mWO2
        try:
          f = k1f*m.pow(cCh4,-0.3e0)*m.pow(cO2,1.3e0)
        except:
          print('[ch4] = {0}\n[O2] = {0}'.format(cCh4,cO2))
          f = 0.0
        return f
# ......................................................................

      nO2p, nN2r, nCO2p, nH2Op, nN2p = reac.stoichCoeff

# ...
#      print (yFrac[0],yFrac[1])
      mWfuel = reac.species['Fuel'].molarMass
      mWO2   = reac.species['O2'].molarMass
      mWCO2  = reac.species['CO2'].molarMass
      mWH2O  = reac.species['H2O'].molarMass
      omega = oneStep(densityMix,yFrac[0],yFrac[1],mWfuel,mWO2,t)
#      print ('omega',t,omega)


      w[0] = -mWfuel*omega
      w[1] = -mWO2*nO2p*omega
      w[2] =  (nN2p-nN2r)*omega
      w[3] =  mWCO2*nCO2p*omega
      w[4] =  mWH2O*nH2Op*omega
# ......................................................................
        
# ......................................................................

# ...   
  @staticmethod
  def taxaLibEnergia(w,reac,t0,t1,fEntalpy = 'hForm'):

#...
    if fEntalpy == 'hComb':
      hC = reac.hCombustionC(t0,t1)/reac.species['Fuel'].molarMass
      q  = hC*w[0]
      print('c',q,hC)
# ......................................................................

# ...
    elif fEntalpy == 'hForm':
      sp0 = reac.species['Fuel']
      sp1 = reac.species['O2']
      sp2 = reac.species['N2']
      sp3 = reac.species['CO2']
      sp4 = reac.species['H2O']
      h0  = Mistura.entalpiaK(sp0,t1)
      h1  = Mistura.entalpiaK(sp1,t1)
      h2  = Mistura.entalpiaK(sp2,t1)
      h3  = Mistura.entalpiaK(sp3,t1)
      h4  = Mistura.entalpiaK(sp4,t1)
      q = -(h0*w[0] + h1*w[1] + h2*w[2] + h3*w[3] + h4*w[4])
# ......................................................................

# ...
    return q
# ......................................................................  

# ...   
  @staticmethod
  def fileRes(f,istep,t,x):

    f.write(' {0:10d} {1:.14f} '.format(istep,t))
    for xi in x:
      f.write(' {0:.14e} '.format(xi))
    f.write('\n')
# ......................................................................  
   

def plot(name,yLable):

  with open(name,'r') as f:
    data=f.read()

  lines = data.splitlines()
  
  t    = []
  xF   = []
  xO2  = []
  xN2  = []
  xCO2 = []
  xH2O = []

  for line in lines:
    _,ti,y1,y2,y3,y4,y5 = line.split()

    t.append(float(ti))
    xF.append(float(y1))
    xO2.append(float(y2))
    xN2.append(float(y3))
    xCO2.append(float(y4))
    xH2O.append(float(y5))

  fig, ax = plt.subplots()
  ax.plot(t,xF,label='Fuel')
  ax.plot(t,xO2,label='O2')
  ax.plot(t,xN2,label='N2')
  ax.plot(t,xCO2,label='CO2')
  ax.plot(t,xH2O,label='H2O')

  ax.set_xlim(0.036,0.037)
  ax.set(xlabel='t(s)',ylabel=yLable)

  ax.legend()
  ax.grid()
  plt.savefig('fig/'+name+'.png')

def plot2(name,yLable):

  with open(name,'r') as f:
    data=f.read()

  lines = data.splitlines()
  
  t = []
  x = []

  for line in lines:
    _,ti,x1 = line.split()

    t.append(float(ti))
    x.append(float(x1))


  fig, ax = plt.subplots()
  ax.set_xlim(0.036,0.037)
  ax.plot(t,x)
  ax.grid()

  ax.set(xlabel='t(s)',ylabel=yLable)
  plt.savefig('fig/'+name+'.png')


def plotTemp(name,temps,yLable):

  with open(name,'r') as f:
    data=f.read()

  lines = data.splitlines()
  
  t  = []
  x  = []

  for line in lines:
    _,ti,x1 = line.split()

    t.append(float(ti))
    x.append(float(x1))


  fig, ax = plt.subplots()
  ax.set_xlim(0.036,0.037)
  ax.plot(t,x,label='Temp')

  x1 = (0.036,0.0037)
  y1 = (temps[0],temps[0])
  ax.plot(x1,y1,label='Max Temp 1')
  y1 = (temps[1],temps[1])
  ax.plot(x1,y1,label='Max Temp 2')
  y1 = (temps[2],temps[2])
  ax.plot(x1,y1,label='Max Temp 3')
  
  ax.set(xlabel='t(s)',ylabel=yLable)
  ax.legend()
  ax.grid()
  plt.savefig('fig/'+name+'.png')
  
def main(argv):

    fileName = argv[1]

    reac,dt,istep = rf.readFile(fileName)
    m    = Mistura(reac,0.26530612,0.7346939,fEst=True,density=7.84)
    print(m.concSpecies('Fuel'),m.concSpecies('O2'))
    print('Mistura inicial\n')
    print('sAir = {0:16.8f}\nphi  = {1:16.8f}\n'.format(m.stoichiometryAirFuel(),m.ratioEquivalente()))
    return
#    t, temp , hs = 0.0, 0.0, 0.0
    t, temp, hs = m.comb(istep=istep,dt=dt,tMix=0.125,eps = 1.0e-4,maxIt = 100,fTime='Euler')

    print('Mistura Final\n')

    tK1, tC1 = m.reac.tempAbiaboticHeatCombustion()
    tK2, tC2 = m.reac.tempAbiaboticHeatFormation()
    tK3, tC3 = m.reac.tempAbiaboticHeatFormationNewtonRapshon()

    dTe1 = abs(tK1 - temp)
    dTe2 = abs(tK2 - temp)
    dTe3 = abs(tK3 - temp)
    
    
    print('{0:<10} = {1}'.format('Time(s)'  , t))
    print('{0:<10} = {1}'.format('hs(KJ/KG)', hs))
    print('{0:<10} = {1}'.format('Temp(K)'  , temp))    

    print('1 {0:.2f} K {1:.2f} C'.format(tK1,tC1))
    print('2 {0:.2f} K {1:.2f} C'.format(tK2,tC2))
    print('3 {0:.2f} K {1:.2f} C'.format(tK3,tC3))
    
    print('dTe1 = {0:.2f}\ndTe2 = {1:.2f}\ndTe3 = {2:.2f}'.format(dTe1,dTe2,dTe3))
    
    plot('data_frac.out','yFrac')
    plot('data_conc.out','kmol/m3')
    plot('data_w.out'   ,'kg/s')
    plot2('data_q.out'  ,'J/s')                      
    plot2('data_hs.out' ,'J/ms')
    plotTemp('data_t.out'  ,(tK1,tK2,tK3),'K')


    yFrac = m.getyFrac()
    
    print(sum(yFrac))

if __name__ == "__main__":    
#    main(sys.argv)
        
    main(('','ch4.dat'))
