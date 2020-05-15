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

def f_SIR_dyn(u,parametres,A):
    """Calcule f(u) selon le modele SIR avec dynamique demographique."""
    #R_0 = beta*Lambda/mu/(mu+gamma)
    S = u[0]
    I = u[2]
    R = u[3]
    N = parametres[0]
    beta = parametres[1]
    gamma = parametres[2]
    Lambda = parametres[3]
    mu = parametres[4]
    dS_dt = Lambda - mu*S - beta*I*S/N
    dE_dt = 0
    dI_dt = beta*I*S/N - gamma*I - mu*I
    dR_dt = gamma*I - mu*R
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
    Lambda = parametres[3]
    mu = parametres[4]
    alpha = parametres[5]
    dS_dt = Lambda - mu*S - beta*I*S/N
    dE_dt = beta*I*S/N-(mu+alpha)*E
    dI_dt = alpha*E - (gamma+mu)*I
    dR_dt = gamma*I - mu*R
    return np.array([dS_dt,dE_dt,dI_dt,dR_dt])



def f_SIR_diff(u,parametres,A):
    """Calcule f(u) selon le modele SIR sans dynamique demographique."""
    #R_0 = beta/gamma
    S = np.array([u[i][0] for i in range(len(A))])
    I = np.array([u[i][2] for i in range(len(A))])
    R = np.array([u[i][3] for i in range(len(A))])
    N = parametres[0]
    beta = parametres[1]
    gamma = parametres[2]
    D = parametres[6]
    dS_dt = -beta*I*S/N
    dE_dt = 0*I
    dI_dt = beta*I*S/N - gamma*I
    dR_dt = gamma*I
    return np.array([[dS_dt[i],dE_dt[i],dI_dt[i],dR_dt[i]] for i in range(len(I))]) + np.dot(A,u)

def f_SEIR_diff(u,parametres,A):
    """Calcule f(u) selon le modele SEIR avec dynamique demographique."""
    #R_0 = alpha/(mu+alpha)*beta/(mu+gamma)
    S = u[0]
    E = u[1]
    I = u[2]
    R = u[3]
    N = parametres[0]
    beta = parametres[1]
    gamma = parametres[2]
    Lambda = parametres[3]
    mu = parametres[4]
    alpha = parametres[5]
    D = parametres[6]
    dS_dt = Lambda - mu*S - beta*I*S/N
    dE_dt = beta*I*S/N-(mu+alpha)*E
    dI_dt = alpha*E - (gamma+mu)*I
    dR_dt = gamma*I - mu*R
    return np.array([[dS_dt[i],dE_dt[i],dI_dt[i],dR_dt[i]] for i in range(len(I))]) + np.dot(A,u)

modeles = [f_SIR,f_SIR_dyn,f_SEIR,f_SIR_diff,f_SEIR_diff]
noms_modeles = ["SIR","SIR","SEIR","SIR","SEIR"]


