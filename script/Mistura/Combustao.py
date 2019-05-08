#!/bin/python3.6

import matplotlib.pyplot as plt
import numpy as np
import Polinomio as pol


def atomMolarMass():

    return {'H':1.0,'O':16,'C':12,'N':14}
      
class specie(object):

    @staticmethod
    def eConv(hf,units,molarMass):
       
        if(units == 'KJ/KMOL'):
          return hf
        else:
          print('conv hf para KJ/KMOL')
          return hf*molarMass

    @staticmethod
    def cpConPol(polCp,units,molarMass):
       
        if(units == 'KJ/KMOLK'):
          print('conv cp para KJ/KGK')
          return [x/molarMass for x in polCp]
        else:     
          return polCp
        

    def __init__(self,cpUnit,hfUnit,chon,polCp,polHf,T0,name):
        self.chon        = chon
        self.molarMass   = self.cMolarMass()
        self.cp          = pol.polinomio(specie.cpConPol(polCp,cpUnit,self.molarMass))
        self.hf          = pol.polinomio(polHf)

        self.name        = name
#        self.hFormation  = specie.eConv(self.hf.value(T0),hfUnit,self.molarMass)

    def entalpyFormation(self):        
        
        return self.hFormation

# ... plot
    def plotCp(self,t0=200.0,tf=3000.0):
        name = self.name
        temp = np.linspace(t0,tf,10)
        y = []
        for t in temp:
            y.append(self.cp.value(t))

        plt.plot(temp,y,label=name)
        plt.legend()
       
# ...
    def isName(self,name):
        return name == self.name
   
    def cMolarMass(self):

        c,h,o,n = self.chon

        W = atomMolarMass()
        m = c*W['C'] + h*W['H'] + o*W['O'] + n*W['N']
        return m
    
# ...
    def __str__(self):        
        str1  = "name: {0}\n".format(self.name)
        str1 += "hFormation: {0:10.2f} kJ/kmol\n".format(self.hFormation)
        str1 += "MolaMass  : {0:10.2f} kg/kmol\n".format(self.molarMass)
        str1 += "Formula   :  C{0}H{1}O{2}N{3}\n".format(*self.chon)
        str1 += self.cp.__str__()
        return str1

  

class reation(object):   

    def __init__(self,sp,pO2,pN2,T0=0.0,T1=0.0):
        self.species     = sp
        self.nSpecies    = len(sp)
        self.arComp      = pO2, pN2
        self.T0          = T0
        self.T1          = T1
        self.stoichCoeff = self.stoichiometricCoeff()
#        self.hCombustion = self.hCombustionC(T0,T1)
        
    def __str__(self):

      str1 = ''
      for sp in self.species:
          str1 += ' '+sp.__str__()

      return str1

    def plotCp(self,t0=298.0,t1=3000.0):
        for _,sp in self.species.items():
            sp.plotCp(t0,t1) 


    def tempAbiaboticHeatCombustion(self):

# determinando as especies
        spFuel = self.species['Fuel']
        spO2   = self.species['O2']
        spCO2  = self.species['CO2']
        spH2O  = self.species['H2O']
        spN2   = self.species['N2']

# 
        nO2, nN2, nCo2r, nH2Or, nN2r = self.stoichCoeff

# massa das especies KG
        mFuel = spFuel.molarMass       
        mCO2  = nCo2r*spCO2.molarMass
        mH2O  = nH2Or*spH2O.molarMass
        mN2   = nN2r*spN2.molarMass

        str1  = '\nTotal mass (kg):\n'
        str1 += 'Fuel : {0}\nCO2  : {1}\nH2O  : {2}\nN2   : {3}'.format(mFuel,mCO2,mH2O,mN2)
        print(str1)

# massa das especies KG reletiva do combustivel        
        mCO2r = mCO2/mFuel
        mH2Or = mH2O/mFuel
        mN2r  = mN2/mFuel

        str1  = '\nMass for mass of fuel (kg/kg):\n'
        str1 += 'CO2 : {0}\nH2O : {1}\nN2  : {2}'.format(mCO2r,mH2Or,mN2r)
        print(str1)

        T0   =  self.T0
        T1   =  self.T1
        T2   =  6000.0

        hC = self.hCombustion/spFuel.molarMass

        print ('\nEntalpia de formacao padrao:\nDHf = {0:10.2f} KJ/KMOL\nDHf = {1:10.2f} KJ/KG'.format(self.hCombustion,hC))

        print ('\nCalculo iterativo')
        
        for i in range(1,50):
# rhs
            Ti   = (T2 + T1)*0.5
            soma = mCO2r*spCO2.cp.value(Ti)+ mH2Or*spH2O.cp.value(Ti) + mN2r*spN2.cp.value(Ti)
            rhs  = -hC + soma*T1
#
            dT = abs(T2 - rhs/soma)
            if (dT < 1.e-8):
                break
            T2 = rhs/soma

#           print ('it = {0:4} T2 = {1:10.4f}'.format(i,T2))

        return T2,reation.convCelsus(T2)


    def tempAbiaboticHeatFormation(self):

# determinando as especies
        spFuel = self.species['Fuel']
        spO2   = self.species['O2']
        spCO2  = self.species['CO2']
        spH2O  = self.species['H2O']
        spN2   = self.species['N2']
# 
        nO2, nN2, nCO2p, nH2Op, nN2p = self.stoichCoeff

# massa das especies KG
        mFuel = spFuel.molarMass       
        mCO2p = nCO2p*spCO2.molarMass
        mH2Op = nH2Op*spH2O.molarMass
        mN2p  = nN2p*spN2.molarMass
        mO2   = nO2*spO2.molarMass
        mN2   = nN2*spN2.molarMass

        str1  = '\nTotal mass (kg):\n'
        str1 += 'Fuel : {0}\nO2 : {0}\nN2 : {1}\nCO2 : {2}\nH2O : {3}\nN2  : {4}'.format(mFuel,mO2,mN2,mCO2p,mH2Op,mN2p)
        print(str1)

# massa das especies KG reletiva do combustivel        
        mCO2pr = mCO2p/mFuel
        mH2Opr = mH2Op/mFuel
        mN2pr  = mN2p/mFuel
        mO2r   = mO2/mFuel
        mN2r   = mN2/mFuel

        str1  = '\nMass for mass of fuel (kg/kg):\n'
        str1 += 'O2 : {0}\nN2 : {1}\nCO2 : {2}\nH2O : {3}\nN2  : {4}'.format(mO2r,mN2r,mCO2pr,mH2Opr,mN2pr)
        print(str1)


        T0   =  self.T0
        T1   =  self.T1
        T2   =  6000.0

        hCO2  = spCO2.entalpyFormation()
        hH2O  = spH2O.entalpyFormation()
        hFuel = spFuel.entalpyFormation() 

        Qm = (nCO2p*hCO2 + nH2Op*hH2O - hFuel)

        print ('\nQm:\nDHf = {0:10.2f} KJ/KMOL\nDHf = {1:10.2f} KJ/KG'.format(Qm,Qm/spFuel.molarMass))

# ... integral do calor especifico de T0 a T1
        intCH4 = spFuel.cp.integral(T0,T1)*spFuel.molarMass
        intO2  = spO2.cp.integral(T0,T1)*spO2.molarMass
        intN2  = spN2.cp.integral(T0,T1)*spN2.molarMass
        Qs      = intCH4 + nN2*intN2 + nO2*intO2
        print ('\nQs:\nQs = {0:10.2f} KJ/KMOL\n'.format(Qs))
        Q = Qs - Qm
        print ('Qs:\nQ = {0:10.2f} KJ/KMOL\n'.format(Q))
        print ('\nCalculo iterativo')
        for i in range(1,50):
# rhs
            Ti   = (T2 + T0)*0.5
            soma = nCO2p*spCO2.cp.value(Ti)*spCO2.molarMass\
                 + nH2Op*spH2O.cp.value(Ti)*spH2O.molarMass\
                 + nN2p*spN2.cp.value(Ti)*spN2.molarMass
            rhs  = Q + soma*T0
#
            dT = abs(T2 - rhs/soma)
            if (dT < 1.e-08):
                break
            T2 = rhs/soma

            print ('it = {0:4} T2 = {1:10.4f}'.format(i,T2))

        return T2,reation.convCelsus(T2)  

    def stoichiometricCoeff(self):

        m = self.species['Fuel'].chon[0]   
        n = self.species['Fuel'].chon[1]
        l = self.species['Fuel'].chon[2]
        p = self.species['Fuel'].chon[2]

        pO2 = self.arComp[0]
        pN2 = self.arComp[1]
        
        nO2  = m + n/4.0 - l/2.0

        nN2  = nO2*pN2/pO2 
        nAir = nO2*(1.e0/pO2)

        nN2r  = nO2*(pN2/pO2) + p
        nCO2r = m
        nH2Or = n/2.0

        nTotal = nCO2r + nH2Or + nN2r;
        pCO2r = nCO2r/nTotal
        pH2Or = nH2Or/nTotal
        pN2r  = nN2r/nTotal
#       nProd = (nN2r/pN2r+nCO2r/pCO2r+nH2Or/pH2Or)/3.0
        nProd = nH2Or/pH2Or

        print('Reacao quimica')

        print('Fuel + {0:.2f} O2 + {1:.2f} N2 -> {2:.2f} CO2 + {3:.2f} H2O + {4:.2f} N2'.format(nO2,nN2,nCO2r,nH2Or,nN2r))

        print('Fuel + {0:.2f} ( O2 +  {1:.2f} N2) -> {2:.2f} CO2 + {3:.2f} H2O + {4:.2f} N2'.format(nAir*pO2,pN2/pO2,nCO2r,nH2Or,nN2r))

        print('Fuel + {0:.2f} ( {1:.2f} O2 +  {2:.2f} N2) -> {6:.2f} ({3:.2f} CO2 + {4:.2f} H2O + {5:.2f} N2)'.format(nAir,pO2,pN2,pCO2r,pH2Or,pN2r,nProd))

        return nO2,nN2,nCO2r,nH2Or,nN2r
    
    def tempAbiaboticHeatFormationNewtonRapshon(self):

# determinando as especies
        spFuel = self.species['Fuel']
        spO2   = self.species['O2']
        spCO2  = self.species['CO2']
        spH2O  = self.species['H2O']
        spN2   = self.species['N2']
# 
        nO2, nN2, nCO2p, nH2Op, nN2p = self.stoichCoeff

# massa das especies KG
        mFuel = spFuel.molarMass       
        mCO2p = nCO2p*spCO2.molarMass
        mH2Op = nH2Op*spH2O.molarMass
        mN2p  = nN2p*spN2.molarMass
        mO2   = nO2*spO2.molarMass
        mN2   = nN2*spN2.molarMass

        str1  = '\nTotal mass (kg):\n'
        str1 += 'Fuel : {0}\nO2 : {0}\nN2 : {1}\nCO2 : {2}\nH2O : {3}\nN2  : {4}'.format(mFuel,mO2,mN2,mCO2p,mH2Op,mN2p)
        print(str1)

# massa das especies KG reletiva do combustivel        
        mCO2pr = mCO2p/mFuel
        mH2Opr = mH2Op/mFuel
        mN2pr  = mN2p/mFuel
        mO2r   = mO2/mFuel
        mN2r   = mN2/mFuel

        str1  = '\nMass for mass of fuel (kg/kg):\n'
        str1 += 'O2 : {0}\nN2 : {1}\nCO2 : {2}\nH2O : {3}\nN2  : {4}'.format(mO2r,mN2r,mCO2pr,mH2Opr,mN2pr)
        print(str1)


        T0   =  self.T0
        T1   =  self.T1
        T2   =  2000.0

        hCO2  = spCO2.entalpyFormation()
        hH2O  = spH2O.entalpyFormation()
        hFuel = spFuel.entalpyFormation() 

        Q0 = hFuel - nCO2p*hCO2 - nH2Op*hH2O

        print ('\nQ0:\nQ0 = {0:10.2f} KJ/KMOL\nQ0 = {1:10.2f} KJ/KG'.format(Q0,Q0/spFuel.molarMass))

# ... integral do calor especifico de T0 a T1
        intCH4 = spFuel.cp.integral(T0,T1)*spFuel.molarMass
        intO2  = spO2.cp.integral(T0,T1)*spO2.molarMass
        intN2  = spN2.cp.integral(T0,T1)*spN2.molarMass
        Qs     = intCH4 + nN2*intN2 + nO2*intO2
        print ('\nQs:\nQs = {0:10.2f} KJ/KMOL\n'.format(Qs))
        Q = Qs + Q0
        print ('Qs:\nQ = {0:10.2f} KJ/KMOL\n'.format(Q))
        print ('\nCalculo iterativo')
        Ti = T0
        for i in range(1,50):
            intCO2 = spCO2.cp.integral(T0,Ti)*spCO2.molarMass
            intH2O = spH2O.cp.integral(T0,Ti)*spH2O.molarMass
            intN2  = spN2.cp.integral(T0,Ti)*spN2.molarMass
#
            soma = nCO2p*intCO2 + nH2Op*intH2O + nN2p*intN2 
#
            f = Q - soma 
#
            soma = nCO2p*spCO2.cp.value(Ti)*spCO2.molarMass\
                 + nH2Op*spH2O.cp.value(Ti)*spH2O.molarMass\
                 + nN2p*spN2.cp.value(Ti)*spN2.molarMass
            df = -soma 
# rhs
            Ti -= f/df
# ...
            if Ti > 6000.0:
                Ti = 6000.0
# .........................................................................................

#
            dT = abs(f/df)
            if (dT < 1.e-08):
                break

            print ('it = {0:4} Ti = {1:10.4f} f ={2:+10.4e} dT ={3:+10.4e}'.format(i,Ti,f,dT))

        T2 = Ti

        return T2,reation.convCelsus(T2)  

    def massRel(self):
      
# determinando as especies
        spFuel = self.species['Fuel']
        spO2   = self.species['O2']
        spCO2  = self.species['CO2']
        spH2O  = self.species['H2O']
        spN2   = self.species['N2']
# 
        nO2, nN2, nCO2p, nH2Op, nN2p = self.stoichCoeff

# massa das especies KG
        mFuel = spFuel.molarMass       
        mCO2p = nCO2p*spCO2.molarMass
        mH2Op = nH2Op*spH2O.molarMass
        mN2p  = nN2p*spN2.molarMass
        mO2   = nO2*spO2.molarMass
        mN2   = nN2*spN2.molarMass
        nF    = self.nameFuel()


        str1  = '\nTotal mass (kg):\n'
        str1 += nF
        str1 += ' : {0}\nO2  : {1}\nN2  : {2}\nCO2 : {3}\nH2O : {4}\nN2  : {5}'.format(mFuel,mO2,mN2,mCO2p,mH2Op,mN2p)
        print(str1)

# massa das especies KG reletiva do combustivel        
        mCO2pr = mCO2p/mFuel
        mH2Opr = mH2Op/mFuel
        mN2pr  = mN2p/mFuel
        mO2r   = mO2/mFuel
        mN2r   = mN2/mFuel

        print('Reacao quimica Kg/Kg de Fuel')

        print(nF + ' + {0:.6f} O2 + {1:.6f} N2 -> {2:.6f} CO2 + {3:.6f} H2O + {4:.6f} N2'.format(mO2r,mN2r,mCO2pr,mH2Opr,mN2pr))   

        return mO2r,mN2r,mCO2pr,mH2Opr,mN2pr

    def nameFuel(self):
        spFuel       = self.species['Fuel']
        c, h, o, n   = spFuel.chon

        name = ''
        if c != 0:
            name += 'C'+str(c)
        if h != 0:
            name += 'H'+str(h)
        if o != 0:
            name += 'H'+str(o)
        if n != 0:
            name += 'N'+str(n)

        return name.strip()

 

    def hCombustionC(self,t0,t1):

        nO2, nN2, nCO2r, nH2Or, nN2r = self.stoichCoeff

# especies
        sCO2  = self.species['CO2']
        sH2O  = self.species['H2O']
        sN2   = self.species['N2']
        sO2   = self.species['O2']
        sFuel = self.species['Fuel']  

# entalpia de formacao padrao (KJ/KG)      
        hCO2  = sCO2.entalpyFormation()
        hH2O  = sH2O.entalpyFormation()
        hFuel = sFuel.entalpyFormation()    

# integral de cp
        intFuel = sFuel.cp.integral(t0,t1)*sFuel.molarMass
        intCO2  = sCO2.cp.integral(t0,t1)*sCO2.molarMass
        intH2O  = sH2O.cp.integral(t0,t1)*sH2O.molarMass
        intN2   = sN2.cp.integral(t0,t1)*sN2.molarMass
        intO2   = sO2.cp.integral(t0,t1)*sO2.molarMass
      
        dh0 = (nCO2r*hCO2 + nH2Or*hH2O) - hFuel

        dCp = intCO2 + nH2Or*intH2O + nN2r*intN2\
             -(intFuel + nO2*intO2 + nN2*intN2)
       
        return dh0 + dCp
   
    @staticmethod
    def convCelsus(tK):
        return tK - 273.15e0

    @staticmethod
    def convKelvin(tC):
        return tC + 273.15e0
        
def plotGraph():
    plt.xlabel('T(Kelvin)')
    plt.ylabel('cp(kJ/kg K)')
    plt.title('Specific heat')
    plt.legend()
    plt.show()

