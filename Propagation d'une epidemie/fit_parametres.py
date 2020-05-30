#Ce script a pour objectif de fitter les paramètres (beta,gamma) d'un modèle SIR

import numpy as np
import matplotlib.pyplot as plt

import os
os.chdir("C:/Users/belfe/OneDrive/Documents/Cours/Ponts_1A_2019-2020/PIR/PIR_Calcul_Scientifique-master/Propagation d'une epidemie/")
import modelisation as mod
import methodes_resolution as met

def f_SIR(u,parametres,A):
    """Calcule f(u) selon le modele SIR sans dynamique demographique."""
    #R_0 = beta/gamma
    S = u[0]
    I = u[2]
    R = u[3]
    N = parametres[0]
    beta = parametres[1]
    gamma = parametres[2]
    dS_dt = -beta*I*S/N
    dE_dt = 0
    dI_dt = beta*I*S/N - gamma*I
    dR_dt = gamma*I
    return np.array([dS_dt,dE_dt,dI_dt,dR_dt])

beta=np.linspace(0.1,0.5,41)
#Beta est souvent pris vers 0.3
gamma=np.linspace(1/20,1/4,41)
#le temps de guérison est entre 4 et 20 jours pour le coronavirus
N_pop=65000000
duree = 120 # Duree (jour)
nb_pts_t = 120
pas_t = duree/nb_pts_t
temps = np.linspace(0,duree,nb_pts_t)

def Euler_explicite(u0,f,parametres,temps,pas_t):
    """Resout du/dt=f(u) avec la methode d'Euler explicite (ou Runge-Kutta d'ordre 1)."""
    U = [u0]
    u = u0

    # Calcul de la matrice A, necessaire uniquement pour la diffusion
    nb_pts_x = len(u0)
    pas_x = 1 / nb_pts_x # car la distance est unitaire
    A = np.zeros((nb_pts_x,nb_pts_x))
    for i in range(nb_pts_x-1):
        A[i,i]=-2
        A[i,i+1]=1
        A[i+1,i]=1
    A[0,0]=-1
    A[-1,-1]=-1
    A/=pas_x**2
    # Résolution
    for t in temps[:-1]:
        u = u + pas_t*f(u,parametres,A)
        U.append(u)
    return np.array(U)

def erreur(t1,t2):
    #Calcule l'erreur quadratique entre les points de t1 et ceux de t2
    return (np.sum((t1-t2)**2))**(1/2)

def calcul_erreur(data):
    #Calcul l'erreur de chaque couple beta,gamma et la stocke dans un tableau
    u0=data[0]
    f=f_SIR
    longueur=data.shape[0]
    A=np.zeros((len(beta),len(gamma)))
    for i in range (len(beta)):
        b=beta[i]
        for j in range(len(gamma)):
            g=gamma[j]
            parametres=[N_pop,b,g]
            U=Euler_explicite(u0,f,parametres,temps,pas_t)
            S,I,R=U[:,0],U[:,2],U[:,3]
            A[i,j]=erreur(S[:longueur],data[:,0])+erreur(I[:longueur],data[:,2])+erreur(R[:longueur],data[:,3])
    return A

def erreurs_inferieures(t, x):
    #Renvoie tous les couples d'indices (i,j) tel que t[i,j]<=x
    L=[]
    (long,larg)=t.shape
    for i in range (long):
        for j in range (larg):
            if t[i,j]<=x:
                L.append((i,j))
    return L

def cas_extremes(indices_admissibles):
    #Parmi les cas admissibles, on cherche celui de plus grand R0=beta/gamma
    #Et celui de plus petit pour avoir le meilleur et le pire cas
    k_min=indices_admissibles[0]
    k_max=k_min
    min=beta[k_min[0]]/gamma[k_min[1]]
    max=beta[k_max[0]]/gamma[k_max[1]]
    for k in indices_admissibles:
        i,j=k
        if beta[i]/gamma[j]<min:
            k_min=k
            min=beta[k_min[0]]/gamma[k_min[1]]
        if beta[i]/gamma[j]>max:
            k_max=k
            max=beta[k_max[0]]/gamma[k_max[1]]
    return [k_min,k_max]


def dessine_probables(data):
    #Dessine les courbes admissibles (dont l'erreur est < à 2*err_opt)
    u0=data[0]
    f=f_SIR
    erreurs=calcul_erreur(data)
    k=np.argmin(erreurs)
    nb_col=erreurs.shape[0]
    i,j=k//nb_col,k%nb_col
    parametres=[N_pop,beta[i],gamma[j]]
    longueur=data.shape[0]

    U=Euler_explicite(u0,f,parametres,temps,pas_t)
    S,I,R=U[:,0],U[:,2],U[:,3]
    err_opt=erreurs[i,j]
    indices_admissibles=erreurs_inferieures(erreurs,2*err_opt)

    #Affichage des données déjà récoltées
    p1=plt.subplot(221)
    p4=plt.subplot(224)
    p1.plot(temps[:longueur],data[:,0],color='blue',label='S')
    p1.plot(temps[:longueur],data[:,2],color='red',label='I')
    p1.plot(temps[:longueur],data[:,3],color='green',label='R')
    p1.legend()

    extremes=cas_extremes(indices_admissibles)
    kmin,kmax=extremes[0],extremes[1]
    #Affichage du cas minimal
    b,g=beta[kmin[0]],gamma[kmin[1]]
    parametres=[N_pop,b,g]
    U=Euler_explicite(u0,f,parametres,temps,pas_t)
    S,I,R=U[:,0],U[:,2],U[:,3]

    p2=plt.subplot(222)
    p2.plot(temps,S,color='blue',label='S')
    p2.plot(temps,I,color='red',label='I')
    p2.plot(temps,R,color='green',label='R')
    p2.legend()

    p4.text(0, 0.8, "beta_min="+str(b))
    p4.text(0,0.7,"gamma_min="+str(g))
    p4.text(0, 0.6, "nombre infectés total :"+str(R[-1]))
    #Affichage du cas maximal
    b,g=beta[kmax[0]],gamma[kmax[1]]
    parametres=[N_pop,b,g]
    U=Euler_explicite(u0,f,parametres,temps,pas_t)
    S,I,R=U[:,0],U[:,2],U[:,3]
    p3=plt.subplot(223)
    p3.plot(temps,S,color='blue',label='S')
    p3.plot(temps,I,color='red',label='I')
    p3.plot(temps,R,color='green',label='R')
    p3.legend()

    p4.text(0, 0.4,"beta_max="+str(b))
    p4.text(0,0.3,"gamma_max="+str(g))
    p4.text(0, 0.2, "nombre infectés total :"+str(R[-1]))

    plt.show()
    return indices_admissibles

data=np.zeros((25,4))
#Données officielles jour par jour à compter du 22/03
#Pour éviter le bruit du début
data[0]=[N_pop-16018,0,16018-(2200+674),2200+674]
data[1]=[N_pop-19856,0,19856-(2200+860),2200+860]
data[2]=[N_pop-22304,0,22304-(3243+1100),3243+1100]
data[3]=[N_pop-25233,0,25233-(3900+1331),3900+1331]
data[4]=[N_pop-29155,0,29155-(4948+1696),4948+1696]
data[5]=[N_pop-32964,0,32964-(5700+1995),5700+1995]
data[6]=[N_pop-37575,0,37575-(5700+2314),5700+2314]
data[7]=[N_pop-40174,0,40174-(7202+2606),7202+2606]
data[8]=[N_pop-44550,0,44550-(7927+3024),7927+3024]
data[9]=[N_pop-52128,0,52128-(9444+3523),9444+3523]
data[10]=[N_pop-56989,0,56989-(10934+4403),10934+4403]
data[11]=[N_pop-59105,0,59105-(12428+5387),12428+5387]
data[12]=[N_pop-64338,0,64338-(14008+6507),14008+6507]
data[13]=[N_pop-68605,0,68605-(15438+7560),15438+7560]
data[14]=[N_pop-70478,0,70478-(16183+8078),16183+8078]
data[15]=[N_pop-74930,0,74930-(17250+8911),17250+8911]
data[16]=[N_pop-78167,0,78167-(19337+10328),19337+10328]
data[17]=[N_pop-82048,0,82048-(21254+10869),21254+10869]
data[18]=[N_pop-86334,0,86334-(23206+12210),23206+12210]
data[19]=[N_pop-90676,0,90676-(24932+13197),24932+13197]
data[20]=[N_pop-93790,0,93790-(26391+13832),26391+13832]
data[21]=[N_pop-120633,0,120633-(27186+14393),27186+14393]
data[22]=[N_pop-124298,0,124298-(27718+14967),27718+14967]
data[23]=[N_pop-129257,0,129257-(28512+15712),28512+15712]
data[24]=[N_pop-132473,0,132473-(30440+17148),30440+17148]
