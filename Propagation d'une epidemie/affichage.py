from tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np
import matplotlib.pyplot as plt

import os
os.chdir("C:/Users/Anthony/Desktop/")
import methodes_resolution as res
import modelisation as mod

# Paramètres par défaut
alpha_def = 0.2
beta_def = 1.4
gamma_def = 0.1
delta_def = 0.5
eta_def = 0.8
D_def = 0.05
# On suppose qu'à l'instant initial il n'y a que S ou I.
S_def = 0.9
E_def = 0.
I_def = 0.1
R_def = 0.
T_def = 0.
u_def = np.array([S_def,E_def,I_def,R_def,T_def])

modeles = mod.modeles # Les modèles sont classés d'abord sans puis avec diffusion

def open(methode, N_pop, duree, pas_t, diffusion, u0):
    class Boutons:
        def __init__(self, master, methode, N_pop, duree, pas_t, diffusion, u0):
            temps = np.linspace(0,duree,int(duree/pas_t))

            # Titre fenetre
            champ_label = Label(master, anchor='n', text="Modélisation de l'évolution d'une population de " + str(N_pop) + " personnes sur " + str(duree) + " jours " + diffusion*"avec" + (1-diffusion)*"sans" + " diffusion")
            champ_label.pack()

            # Choix du modèle
            self.liste = Listbox(master, height=3, width=25)
            self.liste.pack(side=RIGHT)
            self.liste.insert(0, "Modèle SIR")
            self.liste.insert(1, "Modèle SEIR")
            self.liste.insert(2, "Modèle SIRT")
            self.liste.select_set(0)

            # Choix des paramètres d'évolution
            F = Frame(master)
            F1 = LabelFrame(F, text="Paramètres d'évolution")

            self.beta=Scale(F1, orient="horizontal", from_=0, to=4, resolution=0.1, label="Contacts par pers. infectee par jour", tickinterval=1., length=250)
            self.beta.set(beta_def)
            self.beta.pack(side=TOP)

            self.gamma=Scale(F1, orient="horizontal", from_=0, to=1, resolution=0.05, label="Inverse de la duree de guerison", tickinterval=0.25, length=250)
            self.gamma.set(gamma_def)
            self.gamma.pack(side=BOTTOM)

            self.alpha=Scale(F1, orient="horizontal", from_=0, to=1, resolution=0.05, label="Inverse de la duree d'incubation", tickinterval=0.25, length=250)
            self.alpha.set(alpha_def)
            self.alpha.pack(side=BOTTOM)

            self.delta=Scale(F1, orient="horizontal", from_=0, to=1, resolution=0.05, label="Taux d'infection des personnes traitees", tickinterval=0.25, length=250)
            self.delta.set(delta_def)
            self.delta.pack(side=BOTTOM)

            self.eta=Scale(F1, orient="horizontal", from_=0, to=4, resolution=0.1, label="Taux de guerison de traitement", tickinterval=1., length=250)
            self.eta.set(eta_def)
            self.eta.pack(side=BOTTOM)

            # Choix du paramètre de diffusion
            if diffusion == 1:
                self.D=Scale(F1, orient="horizontal", from_=0, to=2, resolution=0.05, label="Parametre de diffusion", tickinterval=2., length=250)
                self.D.set(D_def)
                self.D.pack(side=BOTTOM)

            F1.pack(side=TOP)

            # Choix des paramètres initiaux (sans diffusion uniquement, sinon c'est trop compliqué à choisir)
            if diffusion == 0:
                F2 = LabelFrame(F, text="Paramètre initial")

                self.I=Scale(F2, orient='horizontal', from_=0, to=1, resolution=0.01, label="Taux d'infectés initial", tickinterval=0.25, length=250)
                self.I.set(I_def)
                self.I.pack(side=BOTTOM)

                F2.pack(side=BOTTOM)

            F.pack(side=LEFT)

            # Plot initial
            modele = modeles[3*diffusion] # correspond à SIR car les modèles sont classés d'abord sans puis avec diffusion

            # On affiche une première version par défaut
            fig = Figure()

            if diffusion == 0:
                U = methode(N_pop*u_def,modele,[N_pop,beta_def,gamma_def,alpha_def,delta_def,eta_def],temps,pas_t) / N_pop
                S = U[:,0]
                E = U[:,1]
                I = U[:,2]
                R = U[:,3]
                T = U[:,4]
                self.a = fig.add_subplot(111)
                self.p1,=self.a.plot(temps,S,':',color="blue",label="S")
                self.p2,=self.a.plot(temps,E,':',color="orange",label="E")
                self.p3,=self.a.plot(temps,I,':',color="red",label="I")
                self.p4,=self.a.plot(temps,R,':',color="green",label="R")
                self.p5,=self.a.plot(temps,T,':',color="black",label="T")
                self.a.axis([0-duree/20,duree+duree/20,0-0.05,1+0.05])
                self.a.set_xlabel("Nombre de jours")
                self.a.legend()
                self.a.set_ylabel("Population")
                self.a.set_title("Methode " + res.noms_methodes[res.methodes.index(methode)])
            else:
                U = methode(u0,modele,[N_pop,beta_def,gamma_def,alpha_def,delta_def,eta_def,D_def],temps,pas_t) / N_pop

                temps_a_afficher = [[0,5],[15,30]]
                nb_temps_ligne = len(temps_a_afficher)
                nb_temps_colonne = len(temps_a_afficher[0])
                for i,temps_ligne in enumerate(temps_a_afficher):
                    for j,t in enumerate(temps_ligne):
                        S = U[t][:,0]
                        E = U[t][:,1]
                        I = U[t][:,2]
                        R = U[t][:,3]
                        T = U[t][:,4]
                        self.a = fig.add_subplot(nb_temps_ligne,nb_temps_colonne,i*nb_temps_colonne+j+1)
                        nb_pts_x = len(u0)
                        distance = np.array(range(nb_pts_x))/(nb_pts_x-1)
                        self.p1,=self.a.plot(distance,S,':',color="blue",label="S")
                        self.p2,=self.a.plot(distance,E,':',color="orange",label="E")
                        self.p3,=self.a.plot(distance,I,':',color="red",label="I")
                        self.p4,=self.a.plot(distance,R,':',color="green",label="R")
                        self.p5,=self.a.plot(distance,T,':',color="black",label="T")
                        self.a.axis([0-0.05,1+0.05,0-0.05,1+0.05])
                        self.a.set_xlabel("Distance")
                        self.a.legend()
                        self.a.set_ylabel("Population")
                        self.a.set_title("t = "+"%.3f"%temps[t])
            self.canvas = FigureCanvasTkAgg(fig,master=master)
            self.canvas.get_tk_widget().pack()

            # Bouton pour quitter
            quitter = Button(master, text="Quitter ce super programme", command=master.destroy)
            quitter.pack(side=BOTTOM)

            # Bouton pour actualiser
            actualiser = Button(master, text="Actualiser", command=self.affiche)
            actualiser.pack(side=BOTTOM)

            # Ajout d'un design de qualité
            #master.iconbitmap("enpc_favicon.ico")

        def affiche(self):
            temps = np.linspace(0,duree,int(duree/pas_t))
            # Selection du modele
            modele = modeles[self.liste.curselection()[0] + 3*diffusion] # car dans modeles il y a d'abord 3 modèles sans diffusion puis 3 avec.

            if diffusion == 0:
                # On récupère les choix de l'utilisateur
                parametres = [N_pop,self.beta.get(),self.gamma.get(),self.alpha.get(),self.delta.get(),self.eta.get()]
                u_0 = N_pop*np.array([1-self.I.get(),0,self.I.get(),0,0])

                U = methode(u_0,modele,parametres,temps,pas_t) / N_pop
                S = U[:,0]
                E = U[:,1]
                I = U[:,2]
                R = U[:,3]
                T = U[:,4]

                self.p1.set_ydata(S)
                self.p2.set_ydata(E)
                self.p3.set_ydata(I)
                self.p4.set_ydata(R)
                self.p5.set_ydata(T)
            else:
                # On récupère les choix de l'utilisateur
                parametres = [N_pop,self.beta.get(),self.gamma.get(),self.alpha.get(),self.delta.get(),self.eta.get(),self.D.get()]

                U = methode(u0,modele,parametres,temps,pas_t) / N_pop
                # pour obtenir S pour chaque distance à l'instant t on utilise U[t][:,0]
                temps_a_afficher = [[0,5],[15,30]]
                nb_temps_ligne = len(temps_a_afficher)
                nb_temps_colonne = len(temps_a_afficher[0])
                fig = Figure()
                for i,temps_ligne in enumerate(range(nb_temps_ligne)):
                    for j,t in enumerate(range(nb_temps_colonne)):
                        self.a = fig.add_subplot(nb_temps_ligne,nb_temps_colonne,i*nb_temps_colonne+j+1)
                        S = U[t][:,0]
                        E = U[t][:,1]
                        I = U[t][:,2]
                        R = U[t][:,3]
                        T = U[t][:,4]

                        self.p1.set_ydata(S)
                        self.p2.set_ydata(E)
                        self.p3.set_ydata(I)
                        self.p4.set_ydata(R)
                        self.p5.set_ydata(T)
            self.canvas.draw()

    #Ouverture de la fenetre
    fenetre = Tk()
    fenetre.title("Modélisation de la propagation d'une épidémie")
    fenetre.geometry("1300x600+30+30")

    boutons = Boutons(fenetre, methode, N_pop, duree, pas_t, diffusion, u0)
    fenetre.mainloop()
