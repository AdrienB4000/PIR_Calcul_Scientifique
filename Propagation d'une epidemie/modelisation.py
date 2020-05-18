import numpy as np

"""Modèle épidémologique SIR, SIS, MSIR, SEIR, SEIS, MSEIR, MSEIRS :
S : (susceptible) population susceptible d'etre contaminee
I : (infected) population infectee
R : (recovered) population guerie
E : (exposed) population exposee

Avec u = [S,E,I,R], le modèle donne l'equation : du/dt = f(u) + D*delta_U
avec f qui varie selon la complexite du modele. """

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

def f_SEIR(u,parametres,A):
    """Calcule f(u) selon le modele SEIR avec dynamique demographique."""
    #R_0 = alpha/(mu+alpha)*beta/(mu+gamma)
    S = u[0]
    E = u[1]
    I = u[2]
    R = u[3]
    N = parametres[0]
    beta = parametres[1]
    gamma = parametres[2]
    alpha = parametres[3]
    dS_dt = -beta*I*S/N
    dE_dt = beta*I*S/N - alpha*E
    dI_dt = alpha*E - gamma*I
    dR_dt = gamma*I
    return np.array([dS_dt,dE_dt,dI_dt,dR_dt])



def f_SIR_D(u,parametres,A):
    """Calcule f(u) selon le modele SIR sans dynamique demographique."""
    #R_0 = beta/gamma
    S = u[:,0]
    I = u[:,2]
    R = u[:,3]
    N = parametres[0]
    beta = parametres[1]
    gamma = parametres[2]
    D = parametres[4]
    dS_dt = -beta*I*S/N
    dE_dt = 0*I
    dI_dt = beta*I*S/N - gamma*I
    dR_dt = gamma*I
    # La diffusion se calcule uniquement avec la matrice A, avec condition de gradient nul au bord
    return np.array([[dS_dt[i],dE_dt[i],dI_dt[i],dR_dt[i]] for i in range(len(I))]) + D*np.dot(A,u)

def f_SEIR_D(u,parametres,A):
    """Calcule f(u) selon le modele SEIR avec dynamique demographique."""
    #R_0 = alpha/(mu+alpha)*beta/(mu+gamma)
    S = u[:,0]
    E = u[:,1]
    I = u[:,2]
    R = u[:,3]
    N = parametres[0]
    beta = parametres[1]
    gamma = parametres[2]
    alpha = parametres[3]
    D = parametres[4]
    dS_dt = -beta*I*S/N
    dE_dt = beta*I*S/N - alpha*E
    dI_dt = alpha*E - gamma*I
    dR_dt = gamma*I
    # La diffusion se calcule uniquement avec la matrice A, avec condition de gradient nul au bord
    return np.array([[dS_dt[i],dE_dt[i],dI_dt[i],dR_dt[i]] for i in range(len(I))]) + D*np.dot(A,u)

modeles = [f_SIR,f_SEIR,f_SIR_D,f_SEIR_D]
noms_modeles = ["SIR","SEIR","SIR","SEIR"]

## Jacobienne des fonctions de chaque modèle pour la méthode d'Euler implicite
def d_SIR(u,parametres):
    S = u[0]
    I = u[2]
    N = parametres[0]
    beta = parametres[1]
    gamma = parametres[2]
    return np.array([[-beta*I/N,0,-beta*S/N,0],[0,0,0,0],[beta*I/N,0,beta*S/N-gamma,0],[0,0,gamma,0]])

def d_SEIR(u,parametres):
    S = u[0]
    I = u[2]
    N = parametres[0]
    beta = parametres[1]
    gamma = parametres[2]
    alpha = parametres[3]
    return np.array([[-beta*I/N,0,-beta*S/N,0],[beta*I/N,-alpha,beta*S/N,0],[0,alpha,-gamma,0],[0,0,gamma,0]])

def d_SIR_D(u,parametres):
    S = u[0]
    I = u[2]
    N = parametres[0]
    beta = parametres[1]
    gamma = parametres[2]
    D = parametres[4]
    return np.array([[-beta*I/N,0,-beta*S/N,0],[0,0,0,0],[beta*I/N,0,beta*S/N-gamma,0],[0,0,gamma,0]])

def d_SEIR_D(u,parametres):
    S = u[0]
    I = u[2]
    N = parametres[0]
    beta = parametres[1]
    gamma = parametres[2]
    alpha = parametres[3]
    D = parametres[4]
    return np.array([[-beta*I/N,0,-beta*S/N,0],[beta*I/N,-alpha,beta*S/N,0],[0,aplha,-gamma,0],[0,0,gamma,0]])

jacobiennes = [d_SIR,d_SEIR,d_SIR_D,d_SEIR_D]