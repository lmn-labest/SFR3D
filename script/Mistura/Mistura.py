import readFile  as rf  

def atomMolarMass():

    return {'H': 1.00794e0,'O': 15.9990e0,'C':12.01070e0,'N':14.00670e0 }

 
def cMolarMass(chon):
    """
    Data de criacao  : 07/05/2019
    Data de modificao: 00/00/0000
    ---------------------------------------------------------------
    cMolarMass : Calcula os massa molar
    ---------------------------------------------------------------
    Parametro de entra:
    ---------------------------------------------------------------
    chon   -> tupla/lista com a composicao C H O N do combustivel
    ---------------------------------------------------------------
    Parametro de saida:
    ---------------------------------------------------------------
    retorna a massa molar
    ---------------------------------------------------------------
    """
    c,h,o,n = chon

    W = atomMolarMass()
    m = c*W['C'] + h*W['H'] + o*W['O'] + n*W['N']
    return m

def stoichiometricCoeff(chon,arComp):
    """
    Data de criacao  : 07/05/2019
    Data de modificao: 00/00/0000
    ---------------------------------------------------------------
    stoichiometricCoeff : Calcula os coeficientes estequiometricos
    ---------------------------------------------------------------
    Parametro de entra:
    ---------------------------------------------------------------
    chon   -> tupla/lista com a composicao C H O N do combustivel
    arComp -> tupla/lista com a porcentagem O e N no ar
    ---------------------------------------------------------------
    Parametro de saida:
    ---------------------------------------------------------------
    tupla  -> ceoficientes estequiometricos
              nO2,nN2,nCO2r,nH2Or,nN2r
    ---------------------------------------------------------------
    """
    m,n, l, p = chon 

    pO2,pN2 = arComp
        
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
    if (pN2 == 0.0e0):
        nProd = (nCO2r/pCO2r+nH2Or/pH2Or)/2.0
    else:
        nProd = (nN2r/pN2r+nCO2r/pCO2r+nH2Or/pH2Or)/3.0
#        nProd = nH2Or/pH2Or

    print('Reacao quimica')

    print('Fuel + {0:.2f} O2 + {1:.2f} N2 -> {2:.2f} CO2 + {3:.2f} H2O + {4:.2f} N2'.format(nO2,nN2,nCO2r,nH2Or,nN2r))

    print('Fuel + {0:.2f} ( O2 +  {1:.2f} N2) -> {2:.2f} CO2 + {3:.2f} H2O + {4:.2f} N2'.format(nAir*pO2,pN2/pO2,nCO2r,nH2Or,nN2r))

    print('Fuel + {0:.2f} ( {1:.2f} O2 +  {2:.2f} N2) -> {6:.2f} ({3:.2f} CO2 + {4:.2f} H2O + {5:.2f} N2)'.format(nAir,pO2,pN2,pCO2r,pH2Or,pN2r,nProd))

    return nO2,nN2,nCO2r,nH2Or,nN2r

def arFuelEstiometrico(chon,arComp):
    """
    Data de criacao  : 07/05/2019
    Data de modificao: 00/00/0000
    ---------------------------------------------------------------
    arFuelEstiometrico : Calcula a relacao ar combustivel
    estequiometrico
    ---------------------------------------------------------------
    Parametro de entra:
    ---------------------------------------------------------------
    chon   -> tupla/lista com a composicao C H O N do combustivel
    arComp -> tupla/lista com a porcentagem O e N no ar
    ---------------------------------------------------------------
    Parametro de saida:
    ---------------------------------------------------------------
    tupla  -> fracoes massicas
              yFuel,yAir,yO2,yN2
    ---------------------------------------------------------------
    """

    nO2,nN2, _, _, _ = stoichiometricCoeff(chon,arComp)

    MwO2   = cMolarMass((0,0,2,0))
    MwN2   = cMolarMass((0,0,0,2))
    MwFuel = cMolarMass(chon)
 
    mO2 = nO2*MwO2
    mN2 = nN2*MwN2
    mFuel  = MwFuel

# ... Massa de combustivel / Massa Total 
    yFuel = mFuel/(mFuel + mO2 + mN2)

# ... Massa de O2 / Massa Total 
    yO2   = mO2/(mFuel + mO2 + mN2)
# ... Massa de N2 / Massa Total 
    yN2   = mN2/(mFuel + mO2 + mN2)
    
    return yFuel,1-yFuel,yO2 ,1 - yO2 - yFuel

def fracOxigenAndNitrogen(yAir,arC):

    pO2,pN2  = arC

    mWO2   = cMolarMass((0,0,2,0))
    mWN2   = cMolarMass((0,0,0,2))

    mWAir = pO2*mWO2 + pN2*mWN2

    c1 = pO2*mWO2/mWAir
    c2 = pN2*mWN2/mWAir

    return yAir*c1, yAir*c2

def main(argv):

    fileName = argv[1]

    chon,arC = rf.readFile(fileName)

    yFuel,yAir,yO2,yN2 = arFuelEstiometrico(chon,arC)
#
    print('Mistura estequiometria Ar computivel:')
    print('yFuel = {0:16.8f}\nyAir  = {1:16.8f}\nyO2   = {2:16.8f}\nyN2   = {3:16.8f}\n'.format(yFuel,yAir,yO2,yN2))
#
    yFuel = 0.2
    yO2,yN2 = fracOxigenAndNitrogen(1-yFuel,arC)
    print('yFuel = {0:16.8f}\nO2    = {1:16.8f}\nyN2   = {2:16.8f}\n'.format(yFuel,yO2,yN2))
if __name__ == "__main__":    
#    main(sys.argv)        
    main(('','ch4.dat'))
