import numpy as np

"""Modèle épidémologique SIR, SIS, MSIR, SEIR, SEIS, MSEIR, MSEIRS :
S : (susceptible) population susceptible d'etre contaminee
E : (exposed) population exposee
I : (infected) population infectee
R : (recovered) population guerie
T : (treated) population traitée

Avec u = [S,E,I,R,T], le modèle donne l'equation : du/dt = f(u) + D*laplacien_U
avec f qui varie selon la complexite du modele. """

def f_SIR(u,parametres,A):
    """Calcule f(u) selon le modele SIR."""
    S = u[0]
    I = u[2]
    R = u[3]
    N = parametres[0]
    demo = parametres[1]
    beta = parametres[2]
    gamma = parametres[3]
    dS_dt = demo*N - beta*I*S/N - demo*S
    dI_dt = beta*I*S/N - gamma*I - demo*I
    dR_dt = gamma*I - demo*R
    return np.array([dS_dt,0,dI_dt,dR_dt,0])

def f_SEIR(u,parametres,A):
    """Calcule f(u) selon le modele SEIR."""
    S = u[0]
    E = u[1]
    I = u[2]
    R = u[3]
    N = parametres[0]
    demo = parametres[1]
    beta = parametres[2]
    gamma = parametres[3]
    alpha = parametres[4]
    dS_dt = demo*N - beta*I*S/N - demo*S
    dE_dt = beta*I*S/N - alpha*E - demo*E
    dI_dt = alpha*E - gamma*I - demo*I
    dR_dt = gamma*I - demo*R
    return np.array([dS_dt,dE_dt,dI_dt,dR_dt,0])

def f_SIRT(u,parametres,A):
    """Calcule f(u) selon le modele SIRT."""
    S = u[0]
    I = u[2]
    R = u[3]
    T = u[4]
    N = parametres[0]
    demo = parametres[1]
    beta = parametres[2]
    gamma = parametres[3]
    delta = parametres[5]
    eta = parametres[6]
    tau = parametres[7]
    dS_dt = demo*N - beta*(I+delta*T)*S/N - demo*S
    dI_dt = beta*(I+delta*T)*S/N - (tau+gamma)*I - demo*I
    dR_dt = gamma*I + eta*T - demo*R
    dT_dt = tau*I - eta*T - demo*T

    return np.array([dS_dt,0,dI_dt,dR_dt,dT_dt])



def f_SIR_D(u,parametres,A):
    """Calcule f(u) selon le modele SIR."""
    #R_0 = beta/gamma
    S = u[:,0]
    I = u[:,2]
    R = u[:,3]
    N = parametres[0]
    demo = parametres[1]
    beta = parametres[2]
    gamma = parametres[3]
    D = parametres[8]
    dS_dt = demo*N - beta*I*S/N - demo*S
    dI_dt = beta*I*S/N - gamma*I - demo*I
    dR_dt = gamma*I - demo*R
    # La diffusion se calcule uniquement avec la matrice A, avec condition de gradient nul au bord
    return np.array([[dS_dt[i],0,dI_dt[i],dR_dt[i],0] for i in range(len(I))]) + D*np.dot(A,u)

def f_SEIR_D(u,parametres,A):
    """Calcule f(u) selon le modele SEIR."""
    #R_0 = alpha/(mu+alpha)*beta/(mu+gamma)
    S = u[:,0]
    E = u[:,1]
    I = u[:,2]
    R = u[:,3]
    N = parametres[0]
    demo = parametres[1]
    beta = parametres[2]
    gamma = parametres[3]
    alpha = parametres[4]
    D = parametres[8]
    dS_dt = demo*N - beta*I*S/N - demo*S
    dE_dt = beta*I*S/N - alpha*E - demo*E
    dI_dt = alpha*E - gamma*I - demo*I
    dR_dt = gamma*I - demo*R
    # La diffusion se calcule uniquement avec la matrice A, avec condition de gradient nul au bord
    return np.array([[dS_dt[i],dE_dt[i],dI_dt[i],dR_dt[i],0] for i in range(len(I))]) + D*np.dot(A,u)

def f_SIRT_D(u,parametres,A):
    """Calcule f(u) selon le modele SIRT."""
    S = u[:,0]
    I = u[:,2]
    R = u[:,3]
    T = u[:,4]
    N = parametres[0]
    demo = parametres[1]
    beta = parametres[2]
    gamma = parametres[3]
    delta = parametres[5]
    eta = parametres[6]
    tau = parametres[7]
    D = parametres[8]
    dS_dt = demo*N - beta*(I+delta*T)*S/N - demo*S
    dI_dt = beta*(I+delta*T)*S/N - (tau+gamma)*I - demo*I
    dR_dt = gamma*I + eta*T - demo*R
    dT_dt = tau*I - eta*T - demo*T
    # La diffusion se calcule uniquement avec la matrice A, avec condition de gradient nul au bord
    return np.array([[dS_dt[i],0,dI_dt[i],dR_dt[i],dT_dt[i]] for i in range(len(I))]) + D*np.dot(A,u)

modeles = [f_SIR,f_SEIR,f_SIRT,f_SIR_D,f_SEIR_D,f_SIRT_D]
noms_modeles = ["SIR","SEIR","SIRT","SIR","SEIR","SIRT"]

## Jacobienne des fonctions de chaque modèle pour la méthode d'Euler implicite
def d_SIR(u,parametres):
    """Calcule la jacobienne de f pour le modèle SIR sans diffusion en u."""
    S = u[0]
    I = u[2]
    N = parametres[0]
    demo = parametres[1]
    beta = parametres[2]
    gamma = parametres[3]
    return np.array([[-beta*I/N,0,-beta*S/N,0,0],[0,0,0,0,0],[beta*I/N,0,beta*S/N-gamma,0,0],[0,0,gamma,0,0],[0,0,0,0,0]])

def d_SEIR(u,parametres):
    """Calcule la jacobienne de f pour le modèle SEIR sans diffusion en u."""
    S = u[0]
    I = u[2]
    N = parametres[0]
    demo = parametres[1]
    beta = parametres[2]
    gamma = parametres[3]
    alpha = parametres[4]
    return np.array([[-beta*I/N,0,-beta*S/N,0,0],[beta*I/N,-alpha,beta*S/N,0,0],[0,alpha,-gamma,0,0],[0,0,gamma,0,0],[0,0,0,0,0]])

def d_SIRT(u,parametres):
    """Calcule la jacobienne de f pour le modèle SIRT sans diffusion en u."""
    S = u[0]
    I = u[2]
    T = u[4]
    N = parametres[0]
    demo = parametres[1]
    beta = parametres[2]
    gamma = parametres[3]
    delta = parametres[5]
    eta = parametres[6]
    tau = parametres[7]
    return np.array([[-beta*(I+delta*T)/N,0,-beta*S/N,0,-beta*delta*S/N],[0,0,0,0,0],[beta*(I+delta*T)/N,0,beta*S/N-(tau+gamma),0,beta*delta*S/N],[0,0,gamma,0,eta],[0,0,tau,0,-eta]])

def d_SIR_D(u,parametres):
    """Calcule la jacobienne de f pour le modèle SIR avec diffusion en u."""
    # La jacobienne est diagonale par blocs
    taille = len(u)
    jac = np.zeros((taille,taille))
    for i in range(int(taille/5)):
        jac[5*i:5*i+5,5*i:5*i+5] = d_SIR(u[5*i:5*i+5],parametres)
    return jac

def d_SEIR_D(u,parametres):
    """Calcule la jacobienne de f pour le modèle SEIR avec diffusion en u."""
    # La jacobienne est diagonale par blocs
    taille = len(u)
    jac = np.zeros((taille,taille))
    for i in range(int(taille/5)):
        jac[5*i:5*i+5,5*i:5*i+5] = d_SEIR(u[5*i:5*i+5],parametres)
    return jac

def d_SIRT_D(u,parametres):
    """Calcule la jacobienne de f pour le modèle SIRT avec diffusion en u."""
    # La jacobienne est diagonale par blocs
    taille = len(u)
    jac = np.zeros((taille,taille))
    for i in range(int(taille/5)):
        jac[5*i:5*i+5,5*i:5*i+5] = d_SIRT(u[5*i:5*i+5],parametres)
    return jac

jacobiennes = [d_SIR,d_SEIR,d_SIRT,d_SIR_D,d_SEIR_D,d_SIRT_D]
