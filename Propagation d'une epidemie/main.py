from tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np

import os
os.chdir("C:/Users/Anthony/Desktop/COV_A5/")
import methodes_resolution as res
import modelisation as mod
import affichage as aff

modeles = mod.modeles
noms_modeles = mod.noms_modeles

methodes = res.methodes
noms_methodes = res.noms_methodes

## Paramètres

# Parametres d'affichage
duree = 80 # Duree (jour)
nb_pts = 10000
pas = duree/nb_pts
temps = np.linspace(0,duree,nb_pts)

# Parametres de modelistaion
N = 100000 # Population totale
beta_0 = 1.4 # Contacts par pers. infectee par jour
gamma_0 = 1/12 # Inverse de la duree de guerison
Lambda = 0.012/365*N # Natalite (nb/jour)
mu = 0.009/365 # Taux de mortalite (nb/pers/jour)
alpha_0 = 1/5 # Inverse de la duree d'incubation
# R_0 represente le nombre de contact par pers. infectee, varie selon le choix de f.
D_0 = 1
dimension_espace = 2
diffusion = 0 # 0 si non 1 si oui

parametres = [N,beta_0,gamma_0,Lambda,mu,alpha_0,D_0]

methode = res.Euler_explicite # Methode de resolution

nb_pts_x = 100

A = np.zeros((nb_pts_x,nb_pts_x))
for i in range(nb_pts_x-1):
    A[i,i]=-2
    A[i,i+1]=1
    A[i+1,i]=1
A[0,0]=-1
A[-1,-1]=-1
A/=pas**2

u0=np.array([[0.9*N,0,0.1*N,0] for i in range(nb_pts_x)])
U = res.Euler_explicite(u0,modeles[3],parametres,temps,pas)
print(U)

#aff.open(methode,N,duree,pas,temps,Lambda,mu,dimension_espace)