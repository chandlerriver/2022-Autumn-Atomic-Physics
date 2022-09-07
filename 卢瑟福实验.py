import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp

#Physics costant#
Epsilion = 8.854187817 * 10**(-12)
G = 6.754 * 10**(-11)
RAM = 1.9927 * 10**(-26) #relatice atmoic mass
c = 299792458
pi = 3.141592653
dt = 0.0000000000000001
#Physics constant end #


#
realratio = 10**(-5)
#

class gold():
    def __init__(self,mass=0,radius=0,charge=0,position=0):
        self.mass = 197 * RAM
        self.radius = 134 * 10**(-12) * realratio
        self.charge = 79
        self.position = [50,0]

gold_atomic = gold()

class He4():
    def __init__(self,mass=0,radius=0,charge=0,position=0,verb=0,acceleration=0):
        self.mass = 4.002602 * RAM
        self.radius = 1.67824 * 10**(-15)
        self.charge = 2
        self.position = [-1000,0+10*self.radius]
        self.verb = [0.1*c,0]
        self.acceleration = [0,0]

    def EMF(self,gold_atomic):#电磁作用力
        Ax = self.position[0] - gold_atomic.position[0] 
        Ay = self.position[1] - gold_atomic.position[1]
        A = [Ax,Ay]
        l = (Ax**2 + Ay**2)**0.5
        k = self.charge * gold_atomic.charge /4 /pi /Epsilion /self.mass /(l**3)
        A = [i*k for i in A]
        #print("EMA",end="  ")
        #print(A)
        return A

    def UG(self,gold_atomic):#万有引力
        Ax = self.position[0] - gold_atomic.position[0] 
        Ay = self.position[1] - gold_atomic.position[1]
        A = [Ax,Ay]
        l = (Ax**2 + Ay**2)**0.5
        k = G * gold_atomic.mass /(l**3)
        A = [i*k for i in A]
        #print("UG",end="    ")
        #print(A)
        return A
    
    def cal_acce(self):
        self.acceleration = [0,0]
        E = self.EMF(gold_atomic)
        G = self.UG(gold_atomic)
        
        self.acceleration[0] += E[0]
        self.acceleration[1] += E[1]
        self.acceleration[0] += G[0]
        self.acceleration[1] += G[1]

    def run(self):
        P = [self.verb[i]*dt+self.acceleration[i]*(dt**2) for i in range(2)]
        self.position = [self.position[i]+P[i] for i in range(2)]

He4 = He4()
positionx = []
positiony = []
for i in range(40000000):
    #print("   ")
    He4.cal_acce()
    #print("加速度",end = "  ")
    #print(He4.acceleration)
    He4.run()
    #print("位置",end = "  ")
    #print(He4.position)
    positionx.append(He4.position[0])
    positiony.append(He4.position[1])
plt.plot(positionx,positiony)
plt.scatter(positionx[0],positiony[0],c="red")
plt.scatter(gold_atomic.position[0],gold_atomic.position[1])
plt.show()
