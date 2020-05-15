import numpy as np

def Euler_explicite(u0,f,parametres,temps,pas):
    """Resout du/dt=f(u) avec la methode d'Euler explicite (ou Runge-Kutta d'ordre 1)."""
    U = [u0]
    u = u0
    nb_pts_x = len(u0)
    A = np.zeros((nb_pts_x,nb_pts_x))
    for i in range(nb_pts_x-1):
        A[i,i]=-2
        A[i,i+1]=1
        A[i+1,i]=1
    A[0,0]=-1
    A[-1,-1]=-1
    A/=pas**2
    for t in temps[:-1]:
        u = u + pas*f(u,parametres,A)
        U.append(u)
    return np.array(U)

def Newton(x0,h):
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

def Euler_implicite(u0,f,parametres,temps,pas):
    """Resout du/dt=f(u) avec la methode d'Euler implicite."""
    U = [u0]
    u = u0
    nb_pts_x = len(u0)
    A = np.zeros((nb_pts_x,nb_pts_x))
    for i in range(nb_pts_x-1):
        A[i,i]=-2
        A[i,i+1]=1
        A[i+1,i]=1
    A[0,0]=-1
    A[-1,-1]=-1
    A/=pas**2
    for t in temps[:-1]:
        u = Newton(u,lambda x : np.linalg.norm(u-x+pas*f(u,parametres,A)))
        U.append(u)
    return np.array(U)

def Heun(u0,f,parametres,temps,pas):
    """Resout du/dt=f(u) avec la methode de Heun."""
    U = [u0]
    u = u0
    nb_pts_x = len(u0)
    A = np.zeros((nb_pts_x,nb_pts_x))
    for i in range(nb_pts_x-1):
        A[i,i]=-2
        A[i,i+1]=1
        A[i+1,i]=1
    A[0,0]=-1
    A[-1,-1]=-1
    A/=pas**2
    for t in temps[:-1]:
        u = u + pas/2*(f(u,parametres,A)+f(u+pas*f(u,parametres,A),parametres,A))
        U.append(u)
    return np.array(U)

def Runge_Kutta_2(u0,f,parametres,temps,pas):
    """Resout du/dt=f(u) avec la methode de Runge-Kutta d'ordre 2."""
    U = [u0]
    u = u0
    nb_pts_x = len(u0)
    A = np.zeros((nb_pts_x,nb_pts_x))
    for i in range(nb_pts_x-1):
        A[i,i]=-2
        A[i,i+1]=1
        A[i+1,i]=1
    A[0,0]=-1
    A[-1,-1]=-1
    A/=pas**2
    for t in temps[:-1]:
        u = u + pas*f(u+pas/2*f(u,parametres,A),parametres,A)
        U.append(u)
    return np.array(U)

def Runge_Kutta_4(u0,f,parametres,temps,pas):
    """Resout du/dt=f(u) avec la methode de Runge-Kutta d'ordre 4."""
    U = [u0]
    u = u0
    nb_pts_x = len(u0)
    A = np.zeros((nb_pts_x,nb_pts_x))
    for i in range(nb_pts_x-1):
        A[i,i]=-2
        A[i,i+1]=1
        A[i+1,i]=1
    A[0,0]=-1
    A[-1,-1]=-1
    A/=pas**2
    for t in temps[:-1]:
        k1 = f(u,parametres,A)
        k2 = f(u+pas/2*k1,parametres,A)
        k3 = f(u+pas/2*k2,parametres,A)
        k4 = f(u+pas*k3,parametres,A)
        u = u + pas/6*(k1+2*k2+2*k3+k4)
        U.append(u)
    return np.array(U)

methodes = [Euler_explicite,Euler_implicite,Heun,Runge_Kutta_2,Runge_Kutta_4]
noms_methodes = ["Euler_explicite","Euler_implicite","Heun","Runge_Kutta_2","Runge_Kutta_4"]



