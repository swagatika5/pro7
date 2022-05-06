import numpy as np
import matplotlib.pyplot as plt
import sys

#CONSTANTS

el = 1.6e-19 #charge of electron
k = 1.38e-23
h = 6.626e34
S = 0.5 #spin for fermion
u = 1 # chemical potential for fermions
v = 1
m = 9.1e-31
c = 3e8
pi = 4*np.arctan(1.0)

#FERMIONS

def f_non_rel():
    
    #NON RELATIVISTIC FERMIONS
    def g(E): # Density of state
        return (2*S + 1)*((2*pi*v)*((2*m)**(3/2)))*((E)**(1/2))/(h**3)

    def n(E,T): #F-D distribution
        return 1/(np.exp((E-u)*el/(k*T))+1)

    def f(E,T): #Density of Particle
        return g(E)*n(E,T)
   
    def plotgraph(x,y,xlabelstr,ylabelstr,titlestr):
        plt.plot(x, y, lw=3)
        plt.xlabel(xlabelstr, size=12)
        plt.ylabel(ylabelstr, size=12)
        plt.title(titlestr)
    
    p=2000
    e = np.linspace(0,2,p)
    tl = 100
    th = 1000

    g1 = np.zeros(p)
    for i in range(p):
        g1[i] = g(e[i]) #Density of state
    plotgraph(e,g1,"Energy(ev)","g(e)","Non - Relativistic fermions") 
    plt.show()

    n1 = np.zeros(p)
    n2 = np.zeros(p)
    for i in range(p):
        n1[i] = n(e[i],tl) #F-D distribution #low temp
        n2[i] = n(e[i],th) #high temp
    plotgraph(e,n1,"Energy","n(e)"," ") 
    plotgraph(e,n2,"Energy(ev)","n(e)","Non - Relativistic fermions")
    
    plt.plot(e,n1)
    plt.plot(e,n2)
    plt.legend(["100k","1000k"])
    plt.show()

    f1 = np.zeros(p)
    f2 = np.zeros(p)
    for i in range(p):
        f1[i] = f(e[i],tl) #F-D distribution #low temp
        f2[i] = f(e[i],th) #high temp
    plotgraph(e,f1,"Energy","f(e)"," ") 
    plotgraph(e,f2,"Energy(ev)","f(e)","dN/dE v E for Non - Relativistic Fermions(s = 1/2  u = 1ev)")

    plt.plot(e,f1)
    plt.plot(e,f2)
    plt.legend(["100k","1000k"])
    plt.show()
    
    return None #end non rel fermi

def f_rel():
    #Relativistic Fermions
    #chemical potential as well as energy will be in Mev
    def gr(E):
        return (2*S)*(4*pi*v*(E**2))/((h**3)*(c**3))
    
    def nr(E,T):
        return 1/(np.exp((E-u)*el*(10**6)/(k*T))+1)
    
    def fr(E,T):
        return gr(E)*nr(E,T)
     
    def plotgraph(x,y,xlabelstr,ylabelstr,titlestr):
        plt.plot(x, y, lw=3)
        plt.xlabel(xlabelstr, size=12)
        plt.ylabel(ylabelstr, size=12)
        plt.title(titlestr)
    
    p=2000
    e = np.linspace(0,2,p)
    tlr = 1e8
    thr = 1e9
    
    
    gr1 = np.zeros(p)
    for i in range(p):
        gr1[i] = gr(e[i]) #Density of state
    plotgraph(e,gr1,"Energy(Mev)","g(e)","Relativistic fermions") 
    plt.show()
    
    nr1 = np.zeros(p)
    nr2 = np.zeros(p)
    for i in range(p):
        nr1[i] = nr(e[i],tlr) #F-D distribution #low temp
        nr2[i] = nr(e[i],thr) #high temp
    plotgraph(e,nr1,"Energy","n(e)"," ") 
    plotgraph(e,nr2,"Energy(Mev)","n(e)","Relativistic fermions")
    
    plt.plot(e,nr1)
    plt.plot(e,nr2)
    plt.legend(["10^8k","10^9k"])
    plt.show()
    
    fr1 = np.zeros(p)
    fr2 = np.zeros(p)
    for i in range(p):
        fr1[i] = fr(e[i],tlr) #F-D distribution #low temp
        fr2[i] = fr(e[i],thr) #high temp
    plotgraph(e,fr1,"Energy","f(e)"," ") 
    plotgraph(e,fr2,"Energy(Mev)","f(e)","dN/dE v E for Relativistic Fermions s = 1/2  u = 1Mev")
    
    plt.plot(e,fr1)
    plt.plot(e,fr2)
    plt.legend(["10^8k","10^9k"])
    plt.show()
    
    
    
    return None #end rel fermi

#BOSONS

Sb = 1 #spin for bosons
ub = -1 #chemical potential for bosons
mb = 4*1.66e-27
def b_non_rel():
    #NON RELATIVISTIC BOSONS
    def g(E): # Density of state
        return (2*Sb + 1)*((2*pi*v)*((2*mb)**(3/2)))*((E)**(1/2))/(h**3)

    def n(E,T): #F-D distribution
        return 1/(np.exp((E-ub)*el/(k*T))-1)

    def f(E,T): #Density of Particle
        return g(E)*n(E,T)
   
    def plotgraph(x,y,xlabelstr,ylabelstr,titlestr):
        plt.plot(x, y, lw=3)
        plt.xlabel(xlabelstr, size=12)
        plt.ylabel(ylabelstr, size=12)
        plt.title(titlestr)
    
    p=500
    e = np.linspace(0,0.5,p)
    tl = 100
    th = 1000

    g1 = np.zeros(p)
    for i in range(p):
        g1[i] = g(e[i]) #Density of state
    plotgraph(e,g1,"Energy(ev)","g(e)","Non - Relativistic bosons") 
    plt.show()

    n1 = np.zeros(p)
    n2 = np.zeros(p)
    for i in range(p):
        n1[i] = n(e[i],tl) #F-D distribution #low temp
        n2[i] = n(e[i],th) #high temp
    plotgraph(e,n1,"Energy","n(e)"," ") 
    plt.plot(e,n1)
    plt.legend(["100k"])
    plt.show()
    plotgraph(e,n2,"Energy(ev)","n(e)","Non - Relativistic bosons")
    plt.plot(e,n2)
    plt.legend(["1000k"])
    plt.show()
    
    f1 = np.zeros(p)
    f2 = np.zeros(p)
    for i in range(p):
        f1[i] = f(e[i],tl) #F-D distribution #low temp
        f2[i] = f(e[i],th) #high temp
    plotgraph(e,f1,"Energy","f(e)"," ") 
    plt.plot(e,f1)
    plt.legend(["100k"])
    plt.show()
    plotgraph(e,f2,"Energy(ev)","f(e)","dN/dE v E for Non - Relativistic bosons(s = 1  u = -1ev)")
    plt.plot(e,f2)
    plt.legend(["1000k"])
    plt.show()
    
    
    return None #end non rel bosons

def b_rel():
    
    #RELATIVISTIC BOSONS
    #chemical p[otential and energy in Mev]
    def gr(E): # Density of state
        return (2*Sb)*((4*pi*v)*((E)**(2))/((h**3)*(c**3)))

    def nr(E,T): #F-D distribution
        return 1/(np.exp((E-ub)*el*(10**6)/(k*T))-1)

    def fr(E,T): #Density of Particle
        return gr(E)*nr(E,T)
   
    def plotgraph(x,y,xlabelstr,ylabelstr,titlestr):
        plt.plot(x, y, lw=3)
        plt.xlabel(xlabelstr, size=12)
        plt.ylabel(ylabelstr, size=12)
        plt.title(titlestr)
    
    p=6001
    e = np.linspace(0,6,p)
    tlr = 1e9
    thr = 1e11

    gr1 = np.zeros(p)
    for i in range(p):
        gr1[i] = gr(e[i]) #Density of state
    plotgraph(e,gr1,"Energy(Mev)","g(e)","Relativistic bosons") 
    plt.show()

    nr1 = np.zeros(p)
    nr2 = np.zeros(p)
    for i in range(p):
        nr1[i] = nr(e[i],tlr) #F-D distribution #low temp
        nr2[i] = nr(e[i],thr) #high temp
    plotgraph(e,nr1,"Energy(Mev)","n(e)","Relativistic bosons ") 
    plt.plot(e,nr1)
    plt.legend(["10^9k"])
    plt.show()
    plotgraph(e,nr2,"Energy(Mev)","n(e)","Relativistic bosons")
    plt.plot(e,nr2)
    plt.legend(["10^11k"])
    plt.show()
    
    fr1 = np.zeros(p)
    fr2 = np.zeros(p)
    for i in range(p):
        fr1[i] = fr(e[i],tlr) #b-e distribution #low temp
        fr2[i] = fr(e[i],thr) #high temp
    plotgraph(e,fr1,"Energy(Mev)","f(e)","Relativistic bosons ") 
    plt.plot(e,fr1)
    plt.legend(["10^9k"])
    plt.show()
    plotgraph(e,fr2,"Energy(Mev)","f(e)","dN/dE v E for  Relativistic bosons(s = 1  u = -1Mev)")
    plt.plot(e,fr2)
    plt.legend(["10^11k"])
    plt.show()    
    
    return None #end rel bosons

def default():
    return "Invalid Input"

switcher ={    #switch
    1: f_non_rel,
    2: f_rel,
    3: b_non_rel,
    4: b_rel
    }

def switch(w):
    return switcher.get(w,default)()
z=2
while z>1:
    print("For Non Relativistic Fermion Press......1")
    print("For Relativistic Fermion Press..........2")
    print("For Non Relativistic Boson Press........3")
    print("For Relativistic Boson Press............4")
    print("To exit Press...........................0")
    x =int(input("X = "))
    if x == 1:
        print(switch(1))
    elif x == 2:
        print(switch(2))
    elif x == 3:
        print(switch(3))
    elif x == 4:
        print(switch(4))
    elif x == 0:
        sys.exit()
    else:
        print("Invalid Input......Try again")