##Difference finie explicit
##Exemple : intrusion d'un dike
import numpy as np
import matplotlib.pyplot as plt
import time
from pylab import *

ion()#mode d'interaction 

##Paramètre physique 
L=100 #longueur du modèle en m
Tmagma=1200 #temperature du magma °C
Trock=300 #température de l'encaissant °C
kappa=1e-6 #constante de diffusivité thermal de la roche m²/s
W=50 #epaisseur en m

xmin=-100. #mettre des virgule pour manipuler des float de maniere automatique
xmax=100.
dx=2.
nstep=100
dt=10*3600*24. #soit 10 jours

x=np.linspace(xmin,xmax,nstep)

T=np.ones(len(x))*Trock #initial temperature vector of rocks

#bordure du dike
xleft=0-W/2
xright=0+W/2

for i in range(0,len(x)):
    if xleft<=x[i] and xright>=x[i]:
        T[i]=Tmagma
##        print (T)

t=0 #initialisation du temps

for n in range (0,len(x)):    
    Tnew=np.zeros(len(x))
    Tnew[0]=T[0]
    Tnew[len(x)-1]=T[len(x)-1]
    
for i in range(1,len(x)-1):
    cfl=kappa*dt/(dx**2) #le carré est symbolisé par un double asterix
    Tnew[i]=T[i]+cfl*(T[i+1]-2*T[i]+T[i-1])
    #figure
plt.plot(x,Tnew,'r-')
plt.show() #activation du plot interactif
time.sleep(0.1)
T=Tnew
t=t+dt

####Difference finie explicit
####Exemple : intrusion d'un dike
##import numpy as np
##import matplotlib.pyplot as plt
##import time
##from pylab import *
##
##ion()#mode d'interaction 
##
####Paramètre physique 
##L=100 #longueur du modèle en m
##Tmagma=1200 #temperature du magma °C
##Trock=300 #température de l'encaissant °C
##kappa=1e-6 #constante de diffusivité thermal de la roche m²/s
##W=50 #epaisseur en m
##
##xmin=-100. #mettre des virgule pour manipuler des float de maniere automatique
##xmax=100.
##dx=2.
##nstep=100
##dt=10*3600*24. #soit 10 jours
##
##x=np.linspace(xmin,xmax,nstep)
##
##T=np.ones(len(x))*Trock #initial temperature vector of rocks
##
###bordure du dike
##xleft=0-W/2
##xright=0+W/2
##
##for i in range(0,len(x)):
##    if xleft<=x[i] and xright>=x[i]:
##        T[i]=Tmagma
####        print (T)
##
##t=0 #initialisation du temps
##
##for n in range (0,len(x)):    
##    Tnew=np.zeros(len(x))
##    Tnew[0]=T[0]
##    Tnew[len(x)-1]=T[len(x)-1]
##    
##for i in range(1,len(x)-1):
##    cfl=kappa*dt/(dx**2) #le carré est symbolisé par un double asterix
##    Tnew[i]=T[i]+cfl*(T[i+1]-2*T[i]+T[i-1])
##    #figure
##    draw()#force le dessin de la figure
##    pause(0.05)
##T=Tnew
##t=t+dt #recalcul pour les nouvelles valeurs
##
##ioff() #mode d'interaction off
##show()





    
    
    
    

  
    
    
    



    
    
    
    

  
    
    
    
