import numpy as np
import os
os.chdir("C:/Users/Anthony/Desktop/")
import modelisation as mod



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

def Newton_lent(x0,h):
    """Resout l'equation h(x)=0 en partant de x0 avec la methode de Newton."""
    x = x0
    epsilon = 1e-5
    nb_max = 100
    delta = 1e-7
    for i in range(nb_max):
        dh = np.array([(h([x[j]+delta*(j==i) for j in range(len(x))])-h(x))/delta for i in range(len(x))])
        norm_dh2 = np.linalg.norm(dh)**2
        if norm_dh2<epsilon:
            print("Probleme : impossible de resoudre le modele avec la methode d'Euler implicite : division par 0 lors de la methode de Newton !\n")
            break
        else:
            dx = -h(x)*dh/norm_dh2
            x = x + dx
            if np.linalg.norm(dx) < epsilon:
                break
    return x

def Newton(x0,h,jacobienne,pas_t,parametres):
    """Resout l'equation h(x)=0 en partant de x0 avec la methode de Newton."""
    x = x0
    epsilon = 1e-4
    nb_max = 100
    for i in range(nb_max):
        dh = np.identity(4) - pas_t*jacobienne(x,parametres)
        if np.linalg.matrix_rank(dh)!=4:
            print("Probleme : impossible de resoudre le modele avec la methode d'Euler implicite\n")
            dh = dh + 0.1*np.identity(4)
            break
        else:
            dx = -np.linalg.solve(dh,h(x))
            x = x + dx
            if np.linalg.norm(dx) < epsilon:
                break
    return x

def Euler_implicite(u0,f,parametres,temps,pas_t):
    """Resout du/dt=f(u) avec la methode d'Euler implicite."""
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
        #u = Newton_lent(u,lambda x : np.linalg.norm(u-x+pas_t*f(u,parametres,A)))
        u = Newton(u,lambda x : u-x+pas_t*f(u,parametres,A),mod.jacobiennes[mod.modeles.index(f)],pas_t,parametres)
        U.append(u)
    return np.array(U)

def Heun(u0,f,parametres,temps,pas_t):
    """Resout du/dt=f(u) avec la methode de Heun."""
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
        u = u + pas_t/2*(f(u,parametres,A)+f(u+pas_t*f(u,parametres,A),parametres,A))
        U.append(u)
    return np.array(U)

def Runge_Kutta_2(u0,f,parametres,temps,pas_t):
    """Resout du/dt=f(u) avec la methode de Runge-Kutta d'ordre 2."""
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
        u = u + pas_t*f(u+pas_t/2*f(u,parametres,A),parametres,A)
        U.append(u)
    return np.array(U)

def Runge_Kutta_4(u0,f,parametres,temps,pas_t):
    """Resout du/dt=f(u) avec la methode de Runge-Kutta d'ordre 4."""
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
        k1 = f(u,parametres,A)
        k2 = f(u+pas_t/2*k1,parametres,A)
        k3 = f(u+pas_t/2*k2,parametres,A)
        k4 = f(u+pas_t*k3,parametres,A)
        u = u + pas_t/6*(k1+2*k2+2*k3+k4)
        U.append(u)
    return np.array(U)

methodes = [Euler_explicite,Euler_implicite,Heun,Runge_Kutta_2,Runge_Kutta_4]
noms_methodes = ["Euler_explicite","Euler_implicite","Heun","Runge_Kutta_2","Runge_Kutta_4"]



