#!/usr/bin/python

def main():

    g    = float(input("Gravidade: "))
    cp   = float(input("Calor especifico: "))
    ro   = float(input("Massa especifica: "))
    mi   = float(input("Viscosidade dinamica: "))
    k    = float(input("Condutividade termica: "))
    tq   = float(input("Temperatura da face aquecida: "))
    tf   = float(input("Temperatura da face fria: "))
    tRef = float(input("Temperatura de referencia: "))
    l    = float(input("Distancia entre as faces: "))

    mat = "IdealGas"

    if mat == "water":
        beta = 1.0/ro
    elif mat == "IdealGas":
        beta = 1.0/tRef


    Gr = (ro**2)*beta*g*abs(tq-tf)*(l**3)/(mi**2)
    Pr = mi*cp/k

    print("Prandtl=%e"%Pr) 
    print("Grashof=%e"%Gr)
    print("Rayleigh=%e"%(Gr*Pr))
    
if __name__ == "__main__":
    main()
