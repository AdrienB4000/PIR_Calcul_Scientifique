
## Fitter les parametres (beta,gamma) d'un modele SIR

import numpy as np
import matplotlib.pyplot as plt

import os
import modelisation as mod
import methodes_resolution as res

beta = np.linspace(0.01,0.5,50) # beta est souvent pris vers 0.3
gamma = np.linspace(0.01,0.5,50) # le temps de guerison est entre 4 et 20 jours pour le coronavirus
N_pop = 65000000
duree = 100 # Duree (jour)
nb_pts_t = 100
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
    """Dessine les courbes admissibles (dont l'erreur est < a 2*err_opt)."""
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
    indices_admissibles = erreurs_inferieures(erreurs,1.5*err_opt)

    # Affichage des donnees deja recoltees
    p1 = plt.subplot(221)
    p4 = plt.subplot(224)
    # Le plot de S est inutile puisque S est quasi constant
    #p1.plot(temps[:longueur], data[:,0], color='b', label='S')
    p1.plot(temps[:longueur], data[:,2], color='r', label='I')
    p1.plot(temps[:longueur], data[:,3], color='g', label='R')

    # On trace la meilleure approximation trouvee
    b,g = beta[i],gamma[j]
    print("%.3f"%b)
    print("%.3f"%g)
    parametres = [N_pop,0,b,g]
    U = res.Euler_explicite(u0,f,parametres,temps,pas_t)
    S,I,R=U[:,0],U[:,2],U[:,3]
    p1.plot(temps[:longueur],I[:longueur],linestyle=':',color='red',label="Ie")
    p1.plot(temps[:longueur],R[:longueur],linestyle=':',color='green',label="Ie")

    p1.legend()
    p1.set_title("Les donnees")

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
    p2.set_title("Cas optimiste")

    p4.text(0, 0.8, "beta_min = " + "%.3f"%b)
    p4.text(0, 0.7, "gamma_min = " + "%.3f"%g)
    p4.text(0, 0.6, "nombre infectes total : " + "%.3f"%(N_pop-min(N_pop/(b/g),N_pop)))

    # Affichage du cas maximal
    b,g = beta[kmax[0]],gamma[kmax[1]]
    parametres = [N_pop,0,b,g]
    U = res.Euler_explicite(u0,f,parametres,temps,pas_t)
    S,I,R = U[:,0],U[:,2],U[:,3]
    p3 = plt.subplot(223)
    p3.plot(temps, S, color='b', label='S')
    p3.plot(temps, I, color='r', label='I')
    p3.plot(temps, R, color='g', label='R')
    p3.legend()
    p3.set_title("Cas pessimiste")

    p4.text(0, 0.4, "beta_max = " + "%.3f"%b)
    p4.text(0, 0.3, "gamma_max = " + "%.3f"%g)
    p4.text(0, 0.2, "nombre infectes total : " + "%.3f"%(N_pop-round(min(N_pop/(b/g),N_pop))))

    plt.tight_layout()
    plt.show()
    return indices_admissibles

data=np.zeros((51,5))
#Donnees officielles jour par jour a compter du 22/03
#Pour eviter le bruit du debut
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
data[25]=[N_pop-144944,0,144944-(32297+17901),32297+17901,0]
data[26]=[N_pop-146923,0,146923-(33834+18661),33834+18661,0]
data[27]=[N_pop-146906,0,146906-(35379+19303),35379+19303,0]
data[28]=[N_pop-151808,0,151808-(35973+19694),35973+19694,0]
data[29]=[N_pop-154188,0,154188-(36782+20240),36782+20240,0]
data[30]=[N_pop-156921,0,156921-(38543+20765),38543+20765,0]
data[31]=[N_pop-154715,0,154715-(39988+21309),39988+21309,0]
data[32]=[N_pop-157026,0,157026-(41414+21825),41414+21825,0]
data[33]=[N_pop-158636,0,158636-(42715+22214),42715+22214,0]
data[34]=[N_pop-160292,0,160292-(43816+22583),43816+22583,0]
data[35]=[N_pop-160847,0,160847-(44125+22825),44125+22825,0]
data[36]=[N_pop-164589,0,164589-(44733+23262),44733+23262,0]
data[37]=[N_pop-167605,0,167605-(45997+23629),45997+23629,0]
data[38]=[N_pop-165093,0,165093-(47338+24056),47338+24056,0]
data[39]=[N_pop-165764,0,165764-(48572+24345),48572+24345,0]
data[40]=[N_pop-165764,0,165764-(49300+24563),49300+24563,0]
data[41]=[N_pop-166976,0,166976-(49751+24729),49751+24729,0]
data[42]=[N_pop-167272,0,167272-(49973+24864),49973+24864,0]
data[43]=[N_pop-167886,0,167886-(50438+25168),50438+25168,0]
data[44]=[N_pop-168935,0,168935-(51803+25498),51803+25498,0]
data[45]=[N_pop-172465,0,172465-(53022+25772),53022+25772,0]
data[46]=[N_pop-173040,0,173040-(54076+25949),54076+25949,0]
data[47]=[N_pop-174318,0,174318-(54770+26192),54770+26192,0]
data[48]=[N_pop-174758,0,174758-(54886+26271),54886+26271,0]
data[49]=[N_pop-175027,0,175027-(55062+26341),55062+26341,0]
data[50]=[N_pop-175479,0,175479-(55569+26604),55569+26604,0]
