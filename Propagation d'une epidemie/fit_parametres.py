
## Fitter les paramètres (beta,gamma) d'un modèle SIR

import numpy as np
import matplotlib.pyplot as plt

import os
import modelisation as mod
import methodes_resolution as res

beta = np.linspace(0.1,0.5,41) # beta est souvent pris vers 0.3
gamma = np.linspace(1/20,1/4,41) # le temps de guérison est entre 4 et 20 jours pour le coronavirus
N_pop = 65000000
duree = 120 # Duree (jour)
nb_pts_t = 120
pas_t = duree/nb_pts_t
temps = np.linspace(0, duree, nb_pts_t)

def erreur(t1,t2):
    """Calcule l'erreur quadratique entre les points de t1 et ceux de t2."""
    return np.sum((t1-t2)**2)**(1/2)

def calcul_erreur(data):
    """Calcul l'erreur de chaque couple (beta,gamma) et la stocke dans un tableau."""

    u0 = data[0]
    longueur = len(data)
    A = np.zeros((len(beta),len(gamma)))
    for i,b in enumerate(beta):
        for j,g in enumerate(gamma):
            parametres = [N_pop,0,b,g]
            U = res.Euler_explicite(u0,mod.f_SIR,parametres,temps,pas_t)
            S,I,R = U[:,0],U[:,2],U[:,3]
            A[i,j] = sum(erreur(U[:longueur,k],data[:,k]) for k in [0,2,3])
    return A

def erreurs_inferieures(t,x):
    """Renvoie tous les couples d'indices (i,j) tel que t[i,j]<=x."""
    return [(i,j) for i in range(t.shape[0]) for j in range(t.shape[1]) if t[i,j]<=x]

def cas_extremes(indices_admissibles):
    """Parmi les cas admissibles, on cherche celui de plus grand et de
    plus petit R_0=beta/gamma pour avoir un encadrement de R_0."""

    k_min = indices_admissibles[0]
    k_max = k_min
    R_0_min = beta[k_min[0]]/gamma[k_min[1]]
    R_0_max = beta[k_max[0]]/gamma[k_max[1]]
    for i,j in indices_admissibles:
        R_0 = beta[i]/gamma[j]
        if R_0<R_0_min:
            k_min = i,j
            R_0_min = beta[i]/gamma[j]
        elif R_0>R_0_max:
            k_max = i,j
            R_0_max = beta[i]/gamma[j]
    return [k_min,k_max]


def dessine_probables(data):
    """Dessine les courbes admissibles (dont l'erreur est < à 2*err_opt)."""

    u0 = data[0]
    f = mod.f_SIR
    erreurs = calcul_erreur(data)
    k = np.argmin(erreurs)
    nb_col = erreurs.shape[0]
    i,j = k//nb_col,k%nb_col
    parametres = [N_pop,0,beta[i],gamma[j]]
    longueur = data.shape[0]

    U = res.Euler_explicite(u0,f,parametres,temps,pas_t)
    S,I,R = U[:,0],U[:,2],U[:,3]
    err_opt = erreurs[i,j]
    indices_admissibles = erreurs_inferieures(erreurs,2*err_opt)

    # Affichage des données déjà récoltées
    p1 = plt.subplot(221)
    p4 = plt.subplot(224)
    p1.plot(temps[:longueur], data[:,0], color='b', label='S')
    p1.plot(temps[:longueur], data[:,2], color='r', label='I')
    p1.plot(temps[:longueur], data[:,3], color='g', label='R')
    p1.legend()

    [kmin,kmax] = cas_extremes(indices_admissibles)

    # Affichage du cas minimal
    b,g = beta[kmin[0]],gamma[kmin[1]]
    parametres = [N_pop,0,b,g]
    U = res.Euler_explicite(u0,f,parametres,temps,pas_t)
    S,I,R = U[:,0],U[:,2],U[:,3]

    p2 = plt.subplot(222)
    p2.plot(temps, S, color='b', label='S')
    p2.plot(temps, I, color='r', label='I')
    p2.plot(temps, R, color='g', label='R')
    p2.legend()

    p4.text(0, 0.8, "beta_min = " + "%.3f"%b)
    p4.text(0, 0.7, "gamma_min = " + "%.3f"%g)
    p4.text(0, 0.6, "nombre infectés total : " + "%.3f"%R[-1])

    # Affichage du cas maximal
    b,g = beta[kmax[0]],gamma[kmax[1]]
    parametres = [N_pop,0,b,g]
    U = res.Euler_explicite(u0,f,parametres,temps,pas_t)
    S,I,R = U[:,0],U[:,2],U[:,3]
    p3 = plt.subplot(223)
    p2.plot(temps, S, color='b', label='S')
    p2.plot(temps, I, color='r', label='I')
    p2.plot(temps, R, color='g', label='R')
    p3.legend()

    p4.text(0, 0.4, "beta_max = " + "%.3f"%b)
    p4.text(0, 0.3, "gamma_max = " + "%.3f"%g)
    p4.text(0, 0.2, "nombre infectés total : " + "%.3f"%R[-1])

    plt.show()
    return indices_admissibles


# Données officielles jour par jour à compter du 22/03
# Pour éviter le bruit du début
data = np.zeros((25,5))
data[0]=[N_pop-16018,0,16018-(2200+674),2200+674,0]
data[1]=[N_pop-19856,0,19856-(2200+860),2200+860,0]
data[2]=[N_pop-22304,0,22304-(3243+1100),3243+1100,0]
data[3]=[N_pop-25233,0,25233-(3900+1331),3900+1331,0]
data[4]=[N_pop-29155,0,29155-(4948+1696),4948+1696,0]
data[5]=[N_pop-32964,0,32964-(5700+1995),5700+1995,0]
data[6]=[N_pop-37575,0,37575-(5700+2314),5700+2314,0]
data[7]=[N_pop-40174,0,40174-(7202+2606),7202+2606,0]
data[8]=[N_pop-44550,0,44550-(7927+3024),7927+3024,0]
data[9]=[N_pop-52128,0,52128-(9444+3523),9444+3523,0]
data[10]=[N_pop-56989,0,56989-(10934+4403),10934+4403,0]
data[11]=[N_pop-59105,0,59105-(12428+5387),12428+5387,0]
data[12]=[N_pop-64338,0,64338-(14008+6507),14008+6507,0]
data[13]=[N_pop-68605,0,68605-(15438+7560),15438+7560,0]
data[14]=[N_pop-70478,0,70478-(16183+8078),16183+8078,0]
data[15]=[N_pop-74930,0,74930-(17250+8911),17250+8911,0]
data[16]=[N_pop-78167,0,78167-(19337+10328),19337+10328,0]
data[17]=[N_pop-82048,0,82048-(21254+10869),21254+10869,0]
data[18]=[N_pop-86334,0,86334-(23206+12210),23206+12210,0]
data[19]=[N_pop-90676,0,90676-(24932+13197),24932+13197,0]
data[20]=[N_pop-93790,0,93790-(26391+13832),26391+13832,0]
data[21]=[N_pop-120633,0,120633-(27186+14393),27186+14393,0]
data[22]=[N_pop-124298,0,124298-(27718+14967),27718+14967,0]
data[23]=[N_pop-129257,0,129257-(28512+15712),28512+15712,0]
data[24]=[N_pop-132473,0,132473-(30440+17148),30440+17148,0]

dessine_probables(data)




