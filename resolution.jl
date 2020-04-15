"""Modèle épidémologique SIR, SIS, MSIR, SEIR, SEIS, MSEIR, MSEIRS :
S : (susceptible) population susceptible d'etre contaminee
I : (infected) population infectee
R : (recovered) population guerie
M : (maternal) nouveaux-nes immunise maternellement pour une duree provisoire
E : (exposed) population exposee
V : (vaccinated) population vaccinee

Avec u = [M,S,E,I,R,V], le modèle donne l'equation : du/dt = f(u)
avec f qui varie selon la complexite du modele. """

using PyPlot

## Fonctions f du modele
function f_SIR(u,parametres)
    """Calcule f(u) selon le modele SIR sans dynamique demographique."""
    #R_0 = beta/gamma
    M = u[1]
    S = u[2]
    E = u[3]
    I = u[4]
    R = u[5]
    V = u[6]
    N = parametres[1]
    beta = parametres[2]
    gamma = parametres[3]
    Lambda = parametres[4]
    mu = parametres[5]
    delta = parametres[6]
    epsilon = parametres[7]
    sigma = parametres[8]
    rho = parametres[9]
    dM_dt = 0
    dS_dt = -beta*I*S/N
    dE_dt = 0
    dI_dt = beta*I*S/N - gamma*I
    dR_dt = gamma*I
    dV_dt = 0
    return [dM_dt,dS_dt,dE_dt,dI_dt,dR_dt,dV_dt]
end

function f_SIR_dyn(u,parametres)
    """Calcule f(u) selon le modele SIR avec dynamique demographique."""
    #R_0 = beta*Lambda/mu/(mu+gamma)
    M = u[1]
    S = u[2]
    E = u[3]
    I = u[4]
    R = u[5]
    V = u[6]
    N = parametres[1]
    beta = parametres[2]
    gamma = parametres[3]
    Lambda = parametres[4]
    mu = parametres[5]
    delta = parametres[6]
    epsilon = parametres[7]
    sigma = parametres[8]
    rho = parametres[9]
    dM_dt = 0
    dS_dt = Lambda - mu*S - beta*I*S/N
    dE_dt = 0
    dI_dt = beta*I*S/N - gamma*I - mu*I
    dR_dt = gamma*I - mu*R
    dV_dt = 0
    return [dM_dt,dS_dt,dE_dt,dI_dt,dR_dt,dV_dt]
end

function f_SIS(u,parametres)
    """Calcule f(u) selon le modele SIS sans dynamique demographique."""
    #R_0 = beta/gamma
    M = u[1]
    S = u[2]
    E = u[3]
    I = u[4]
    R = u[5]
    V = u[6]
    N = parametres[1]
    beta = parametres[2]
    gamma = parametres[3]
    Lambda = parametres[4]
    mu = parametres[5]
    delta = parametres[6]
    epsilon = parametres[7]
    sigma = parametres[8]
    rho = parametres[9]
    dM_dt = 0
    dS_dt = -beta*I*S/N + gamma*I
    dE_dt = 0
    dI_dt = beta*I*S/N - gamma*I
    dR_dt = 0
    dV_dt = 0
    return [dM_dt,dS_dt,dE_dt,dI_dt,dR_dt,dV_dt]
end

function f_MSIR(u,parametres)
    """Calcule f(u) selon le modele MSIR avec dynamique demographique."""
    M = u[1]
    S = u[2]
    E = u[3]
    I = u[4]
    R = u[5]
    V = u[6]
    N = parametres[1]
    beta = parametres[2]
    gamma = parametres[3]
    Lambda = parametres[4]
    mu = parametres[5]
    delta = parametres[6]
    epsilon = parametres[7]
    sigma = parametres[8]
    rho = parametres[9]
    dM_dt = Lambda - delta*M - mu*M
    dS_dt = delta*M - beta*I*S/N - mu*S
    dE_dt = 0
    dI_dt = beta*I*S/N - (gamma+mu)*I
    dR_dt = gamma*I - mu*R
    dV_dt = 0
    return [dM_dt,dS_dt,dE_dt,dI_dt,dR_dt,dV_dt]
end

function f_SEIR(u,parametres)
    """Calcule f(u) selon le modele SEIR avec dynamique demographique."""
    #R_0 = a/(mu+a)*beta/(mu+gamma)
    M = u[1]
    S = u[2]
    E = u[3]
    I = u[4]
    R = u[5]
    V = u[6]
    N = parametres[1]
    beta = parametres[2]
    gamma = parametres[3]
    Lambda = parametres[4]
    mu = parametres[5]
    delta = parametres[6]
    epsilon = parametres[7]
    sigma = parametres[8]
    rho = parametres[9]
    dM_dt = 0
    dS_dt = Lambda - mu*S - beta*I*S/N
    dE_dt = beta*I*S/N-(mu+epsilon)*E
    dI_dt = epsilon*E - (gamma+mu)*I
    dR_dt = gamma*I - mu*R
    dV_dt = 0
    return [dM_dt,dS_dt,dE_dt,dI_dt,dR_dt,dV_dt]
end

function f_SEIS(u,parametres)
    """Calcule f(u) selon le modele SEIS avec dynamique demographique."""
    M = u[1]
    S = u[2]
    E = u[3]
    I = u[4]
    R = u[5]
    V = u[6]
    N = parametres[1]
    beta = parametres[2]
    gamma = parametres[3]
    Lambda = parametres[4]
    mu = parametres[5]
    delta = parametres[6]
    epsilon = parametres[7]
    sigma = parametres[8]
    rho = parametres[9]
    dM_dt = 0
    dS_dt = Lambda - mu*S - beta*I*S/N +gamma*I
    dE_dt = beta*I*S/N-(mu+epsilon)*E
    dI_dt = epsilon*E - (gamma+mu)*I
    dR_dt = 0
    dV_dt = 0
    return [dM_dt,dS_dt,dE_dt,dI_dt,dR_dt,dV_dt]
end

function f_MSEIR(u,parametres)
    """Calcule f(u) selon le modele MSEIR avec dynamique demographique."""
    M = u[1]
    S = u[2]
    E = u[3]
    I = u[4]
    R = u[5]
    V = u[6]
    N = parametres[1]
    beta = parametres[2]
    gamma = parametres[3]
    Lambda = parametres[4]
    mu = parametres[5]
    delta = parametres[6]
    epsilon = parametres[7]
    sigma = parametres[8]
    rho = parametres[9]
    dM_dt = Lambda - (delta+mu)*M
    dS_dt = delta*M - mu*S - beta*I*S/N
    dE_dt = beta*I*S/N-(mu+epsilon)*E
    dI_dt = epsilon*E - (gamma+mu)*I
    dR_dt = gamma*I - mu*R
    dV_dt = 0
    return [dM_dt,dS_dt,dE_dt,dI_dt,dR_dt,dV_dt]
end

function f_MSEIRS(u,parametres)
    """Calcule f(u) selon le modele MSEIRS avec dynamique demographique."""
    M = u[1]
    S = u[2]
    E = u[3]
    I = u[4]
    R = u[5]
    V = u[6]
    N = parametres[1]
    beta = parametres[2]
    gamma = parametres[3]
    Lambda = parametres[4]
    mu = parametres[5]
    delta = parametres[6]
    epsilon = parametres[7]
    sigma = parametres[8]
    rho = parametres[9]
    dM_dt = Lambda - (delta+mu)*M
    dS_dt = delta*M - mu*S - beta*I*S/N + sigma*R
    dE_dt = beta*I*S/N-(mu+epsilon)*E
    dI_dt = epsilon*E - (gamma+mu)*I
    dR_dt = gamma*I - (mu+sigma)*R
    dV_dt = 0
    return [dM_dt,dS_dt,dE_dt,dI_dt,dR_dt,dV_dt]
end

modeles = [f_SIR,f_SIR_dyn,f_SIS,f_MSIR,f_SEIR,f_SEIS,f_MSEIR,f_MSEIRS]
noms_modeles = ["SIR","SIR","SIS","MSIR","SEIR","SEIS","MSEIR","MSEIRS"]
modeles_avec_dyn = [f_SIR_dyn,f_MSIR,f_SEIR,f_SEIS,f_MSEIR,f_MSEIRS]
modeles_sans_dyn = [f_SIR,f_SIS]


## Methodes de resolution
function Euler_explicite(u0,f,parametres)
    """Resout le modele SIR avec la methode d'Euler explicite (ou Runge-Kutta d'ordre 1)."""
    U=zeros(Float64, nb_pts+1, length(u0))
    U[1,:]=u0
    u = u0
    i=1
    for t in temps[1:length(temps)-1]
        i+=1
        u = u + pas*f(u,parametres)
        U[i,:]=u
    end
    return U
end

function Newton(x0,h)
    """Resout l'equation h(x)=0 en partant de x0 avec la methode de Newton."""
    x = x0
    eps = 1e-5
    nb_max = 100
    e = 1e-7
    for i in range(nb_max)
        dh = [(h([x[j]+e*(j==i) for j in 1:len(x)])-h(x))/e for i in 1:len(x)]
        norm_dh2 = norm(dh)^2
        if norm_dh2 < eps
            print("Probleme : impossible de resoudre le modele avec la methode d'Euler implicite : division par 0 lors de la methode de Newton !\n")
            break
        else
            dx = -h(x)*dh/norm_dh2
            x = x + dx
            if norm(dx) < eps
                break
            end
        end
    end
    return x
end

function Euler_implicite(u0,f,parametres)
    """Resout le modele SIR avec la methode d'Euler implicite."""
    U=zeros(Float64, nb_pts+1, length(u0))
    U[1,:]=u0
    u = u0
    i=1
    for t in temps[1:length(temps)-1]
        i+=1
        u = Newton(u, x -> norm(u-x+pas*f(u,parametres)))
        U[i,:]=u
    end
    return U
end

function Heun(u0,f,parametres)
    """Resout le modele SIR avec la methode de Heun."""
    U=zeros(Float64, nb_pts+1, length(u0))
    U[1,:]=u0
    u = u0
    i=1
    for t in temps[1:length(temps)-1]
        i+=1
        u = u + pas/2*(f(u,parametres)+f(u+pas*f(u,parametres),parametres))
        U[i,:]=u
    end
    return U
end

function Runge_Kutta_2(u0,f,parametres)
    """Resout le modele SIR avec la methode de Runge-Kutta d'ordre 2."""
    U=zeros(Float64, nb_pts+1, length(u0))
    U[1,:]=u0
    u = u0
    i=1
    for t in temps[1:length(temps)-1]
        i+=1
        u = u + pas*f(u+pas/2*f(u,parametres),parametres)
        U[i,:]=u
    end
    return U
end

function Runge_Kutta_4(u0,f,parametres)
    """Resout le modele SIR avec la methode de Runge-Kutta d'ordre 4."""
    U=zeros(Float64, nb_pts+1, length(u0))
    U[1,:]=u0
    u = u0
    i=1
    for t in temps[1:length(temps)-1]
        i+=1
        k1 = f(u,parametres)
        k2 = f(u+pas/2*k1,parametres)
        k3 = f(u+pas/2*k2,parametres)
        k4 = f(u+pas*k3,parametres)
        u = u + pas/6*(k1+2*k2+2*k3+k4)
        U[i,:]=u
    end
    return U
end

methodes = [Euler_explicite,Euler_implicite,Heun,Runge_Kutta_2,Runge_Kutta_4]
noms_methodes = ["Euler_explicite","Euler_implicite","Heun","Runge_Kutta_2","Runge_Kutta_4"]

## Trace
function trace(U,modele=nothing,methode=nothing,cree=true)
    """Trace les solutions obtenues u en ouvrant une nouvelle fenetre si cree=True."""
    M = U[:,1]
    S = U[:,2]
    E = U[:,3]
    I = U[:,4]
    R = U[:,5]
    V = U[:,6]
    titre = "Evolution de l'epidemie au cours du temps"
    if cree
        plot(temps,M,color="yellow",label="M",linestyle="--")
        plot(temps,S,color="blue",label="Susceptibles",linestyle=":")
        plot(temps,E,color="orange",label="Exposés",linestyle=":")
        plot(temps,I,color="red",label="Infectés",linestyle=":")
        plot(temps,R,color="green",label="Guéris",linestyle=":")
        plot(temps,V,color="magenta",label="Vaccinés",linestyle=":")
        if modele != nothing
            titre *= " modele " * (noms_modeles[findall(i->i==modele,modeles)])[1]
        end
        if methode != nothing
            titre *= " (methode " * (noms_methodes[findall(i->i==methode,methode)])[1] * ")"
        end
        PyPlot.legend(loc="upper right")
    else
        plot(temps,M,color="yellow")
        plot(temps,S,color="blue")
        plot(temps,E,color="orange")
        plot(temps,I,color="red")
        plot(temps,R,color="green")
        plot(temps,V,color="magenta")
    end
    PyPlot.title(titre)
    grid("on")
    xlabel("Temps en jours")
    ylabel("Population")
end


## Resolution

function solution(modele,u0,cree)
    """Trace la solution exacte du modele selon la fonction f du modele choisie en ouvrant une nouvelle fenetre si cree=True."""
    U=[u0]
    for t in temps[1:length(temps)-1]
        #A completer
        U.append(u0)
    end
    trace(U,modele,cree)
    return U
end

function resolution(methode,modele,parametres,u0,cree)
    """Resout numeriquement le modele a partir de u0 selon la fonction f du modele et avec la methode de
    resolution choisies et trace la solution en ouvrant une nouvelle fenetre si cree=True."""
    U = methode(u0,modele,parametres)
    trace(U,modele,nothing,cree)
    return U
end


## Parametres

# Parametres d'affichage
duree = 80 # Duree (jour)
nb_pts = 10000
pas = duree/nb_pts
temps = collect(0:pas:duree)

# Parametres de modelistaion
N = 100000 # Population totale
beta = 1.2 # Contact par pers. infectee par jour
gamma = 1/8 # Inverse de la duree de guerison
Lambda = 0.012/365*N # Natalite (nb/jour)
mu = 0.009/365 # Taux de mortalite (nb/pers/jour)
delta = 0.012/365 # Taux de desimmunisation (nb/pers/jour)
epsilon = 1/5 # Inverse de la duree d'incubation
sigma = 1/30 # Inverse de la duree d'immunite (modele MSEIRS)
rho = 0.01 # Taux de letalite de l'epidemie
# R_0 represente le nombre de contact par pers. infectee, varie selon le choix de f.

parametres = [N,beta,gamma,Lambda,mu,delta,epsilon,sigma,rho]

# Parametres de resolution
I0 = 10
proportion_M0_dans_N = 0.1
M0 = 0*proportion_M0_dans_N*N
u0 = [M0,N-I0-M0,0,I0,0,0] # [M,S,E,I,R,V]
methode = Runge_Kutta_4 # Methode de resolution
modele = f_SEIS # Modele d'epidemie
cree = true # Affiche le resultat sur une nouvelle fenetre


## Affichage

resolution(methode,modele,parametres,u0,cree)
