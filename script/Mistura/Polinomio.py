class polinomio(object):

    def __init__(self,a):
        self.a      = a
        self.degree = len(a)-1

# ... calcula o calor especifico
    def value(self,t):
        a = self.a 
        n = self.degree

        tmp = a[0]

        for i in range(1,n+1):
            tmp += a[i]*t**i

        return tmp

# ... 
    def __str__(self):
        
        a    = self.a
        n    = self.degree

        str1  = "pol(T) = a0 + a1*T + a2*T² + ... \n"
        str1 += "nPol = {0}\n".format(n)
        for i in range(0,n+1):
            str1 += "a{0} = {1}\n".format(i,a[i]) 
        return str1 

# ... integral definida do polinomio
    def integral(self,t0,t1):
        a = self.a 
        n = self.degree
        tmp = 0.0
        for i in range(0,n+1):
            d = float(i+1)
            A = a[i]/d
            tmp += A*(pow(t1,d)-pow(t0,d))

        return tmp

# ... integral definida do polinomio
    def integralNum(self,t0,t1,m=10):
            
        h = (t1-t0)/m
        x = t0
        soma = self.value(x)
        for i in range(1,m):
            x += h
            soma += 2.0*self.value(x)
        soma += self.value(t1)
                   
        return soma*h*0.5


    
