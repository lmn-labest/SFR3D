import math as m

def main():
  
  Q   = 30e+03  #w
  rho = 1.2     # kg/m3
  cp  = 1000    # J/kg K
  T   = 293     # K
  g   = 9.81    # m/s2 

  D = (Q/(rho*cp*T*m.sqrt(g)))**(2.0/5.0)

  print('Grosseira {0} cm\nFina      {1} cm'.format(100.0*D/4.0,100.0*D/16.0))

main()



