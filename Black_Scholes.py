"""Pricing a Black Scholes Option"""
   
#opt_type - C, P      
#S - Stock Price
#K - Option Strike
#v - Volatility, annualized
#r - Risk free rate, annualized
#d - div yield
#t - days to expiration
import math    
from scipy.stats import norm  

# S = 100
# K = 100
# v = .23
# r = 0.05
# d = 0.01
# t = 30
        
class BSOption(object): # creates a new class like int, str, list, etc
    def __init__(self, opt_type,S=100,K=100,v=.2,r=.05,d=.01,t=30):  #creates a constructor method needs to be here
        self.opt_type = opt_type
        self.S = S
        self.K = K
        self.v = v
        self.r = r
        self.d = d
        self.t = t
        print('Pricing a Black Scholes {}'.format(opt_type))
       
    def d1(self):
        d1 = (math.log(self.S / self.K) + (self.r - self.d + 0.5 * self.v ** 2)* (self.t/252)) / (self.v * math.sqrt(self.t/252))
        return d1

    def d2(self):
        d2 = self.d1() - self.v*math.sqrt(self.t/252)
        return d2
       
    def Nd1c(self):
        Nd1c = norm.cdf(self.d1())
        return Nd1c
    
    def Nd2c(self):
        Nd2c = norm.cdf(self.d2())
        return Nd2c
    
    def Nd1p(self):
        Nd1p = norm.cdf(-1*self.d1())
        return Nd1p
    
    def Nd2p(self):
        Nd2p = norm.cdf(-1*self.d2())
        return Nd2p      
         
    def price(self):
                
        if self.opt_type == 'C':
            price = self.S*math.exp(-self.d*self.t/252)*self.Nd1c() - self.K*math.exp(-self.r*self.t/252)*self.Nd2c()
        elif self.opt_type == 'P':
            price = self.K*math.exp(-self.r*self.t/252)*self.Nd2p()-self.S*math.exp(-self.d*self.t/252)*self.Nd1p()
        
        return price
    
    def delta(self):

        if self.opt_type == 'C':
            delta = math.exp(-self.d*self.t/252)*self.Nd1c()
        elif self.opt_type == 'P':
            delta = math.exp(-self.d*self.t/252)*(self.Nd1c()-1)
        
        return delta
    
    def gamma(self):      
        gamma = (math.exp(-self.d*self.t/252)/(self.S*self.v*math.sqrt(self.t/252)))*norm.pdf(self.d1())    
        return gamma
    
    def vega(self):      
        vega = 1/100 * self.S*math.exp(-self.d*self.t/252)*math.sqrt(self.t/252)*norm.pdf(self.d1()) 
        return vega
    
    def theta(self):
     
        if self.opt_type == 'C':       
     
            c1 = -(self.S*self.v*math.exp(-self.d*self.t/252)*self.Nd1c())/(2*math.sqrt(self.t/252))
            c2 = -self.r*self.K*math.exp(-self.r*self.t/252)*self.Nd2c()
            c3 = self.d*self.S*math.exp(-self.d*self.t/252)*self.Nd1c()
            theta = (1/252)*(c1+c2+c3)
                           
        elif self.opt_type == 'P':
                    
            p1 = -(self.S*self.v*math.exp(-self.d*self.t/252)*self.Nd1c())/(2*math.sqrt(self.t/252))
            p2 = self.r*self.K*math.exp(-self.r*self.t/252)*self.Nd2p()
            p3 = -self.d*self.S*math.exp(-self.d*self.t/252)*self.Nd1p()
            theta = (1/252)*(p1+p2+p3)
        
        return theta
    
    def rho(self):

        if self.opt_type == 'C':
            rho = (1/100)*self.K*(self.t/252)*math.exp(-self.r*self.t/252)*self.Nd2c()
        elif self.opt_type == 'P':
            rho = -(1/100)*self.K*(self.t/252)*math.exp(-self.r*self.t/252)*self.Nd2p()
        return rho
    
    