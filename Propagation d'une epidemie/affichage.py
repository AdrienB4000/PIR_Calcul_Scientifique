from tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np

import os
import methodes_resolution as res
import modelisation as mod

# Paramètres par défaut
beta_def = 1.4
gamma_def = 0.1
alpha_def = 0.2
delta_def = 0.5
eta_def = 0.8
tau_def = 0.2
D_def = 0.05
# On suppose qu'à l'instant initial il n'y a que S ou I.
S_def = 0.9
E_def = 0.
I_def = 0.1
R_def = 0.
T_def = 0.
u_def = np.array([S_def,E_def,I_def,R_def,T_def])

labels = ["S","E","I","R","T"]
colors = ["blue","orange","red","green","black"]

def open(N_pop, duree, pas_t, u0):
    class Boutons:
        def __init__(self, master, N_pop, duree, pas_t, u0):
            temps = np.linspace(0,duree,int(duree/pas_t))
            self.diffusion = 0 # par défaut sans diffusion

            # Titre fenetre
            titre = "Modélisation de l'évolution d'une population de " + str(N_pop) + " personnes sur " + str(duree) + " jours " + self.diffusion*"avec" + (1-self.diffusion)*"sans" + " diffusion"
            champ_label = Label(master, anchor='n', text=titre)
            champ_label.pack()

            # Choix du modèle
            self.liste = Listbox(master, height=3, width=25, exportselection=0)
            self.liste.pack(side=RIGHT)
            self.liste.insert(0, "Modèle SIR")
            self.liste.insert(1, "Modèle SEIR")
            self.liste.insert(2, "Modèle SIRT")
            self.liste.select_set(0)

            # Choix de la méthode de résolution
            self.liste_methode = Listbox(master, height=len(res.noms_methodes), width=20, exportselection=0)
            self.liste_methode.pack(side=RIGHT)
            for i in range(len(res.noms_methodes)):
                self.liste_methode.insert(i, res.noms_methodes[i])
            self.liste_methode.select_set(0)

            # Choix des paramètres d'évolution
            F = Frame(master)
            F1 = LabelFrame(F, text="Paramètres d'évolution")

            # Choix de la diffusion
            self.var1 = IntVar()
            self.checkbox_diffusion = Checkbutton(F1, text="diffusion", onvalue=1, offvalue=0, variable=self.var1)
            self.checkbox_diffusion.pack(side=TOP)

            self.beta=Scale(F1, orient="horizontal", from_=0, to=4, resolution=0.1, label="Contacts par pers. infectee par jour", tickinterval=1., length=320)
            self.beta.set(beta_def)
            self.beta.pack(side=TOP)

            self.gamma=Scale(F1, orient="horizontal", from_=0, to=1, resolution=0.05, label="Taux de guerison", tickinterval=0.25, length=320)
            self.gamma.set(gamma_def)
            self.gamma.pack(side=TOP)

            self.alpha=Scale(F1, orient="horizontal", from_=0, to=1, resolution=0.05, label="Taux d'incubation (SEIR)", tickinterval=0.25, length=320)
            self.alpha.set(alpha_def)
            self.alpha.pack(side=TOP)

            self.delta=Scale(F1, orient="horizontal", from_=0, to=1, resolution=0.05, label="Taux d'infection des personnes traitees (SIRT)", tickinterval=0.25, length=320)
            self.delta.set(delta_def)
            self.delta.pack(side=TOP)

            self.eta=Scale(F1, orient="horizontal", from_=0, to=2, resolution=0.1, label="Taux de guerison du traitement (SIRT)", tickinterval=1., length=320)
            self.eta.set(eta_def)
            self.eta.pack(side=TOP)

            self.tau=Scale(F1, orient="horizontal", from_=0, to=2, resolution=0.1, label="Taux de traitement (SIRT)", tickinterval=1., length=320)
            self.tau.set(tau_def)
            self.tau.pack(side=TOP)

            # Choix du paramètre de diffusion
            self.D=Scale(F1, orient="horizontal", from_=0, to=0.1, resolution=0.005, label="Parametre de diffusion", tickinterval=2., length=320)
            self.D.set(D_def)
            self.D.pack(side=BOTTOM)

            F1.pack(side=TOP)

            # Choix des paramètres initiaux (sans diffusion uniquement, sinon c'est trop compliqué à choisir)
            F2 = LabelFrame(F, text="Paramètre initial")

            self.I=Scale(F2, orient='horizontal', from_=0, to=1, resolution=0.01, label="Taux d'infectés initial", tickinterval=0.25, length=320)
            self.I.set(I_def)
            self.I.pack(side=BOTTOM)

            F2.pack(side=BOTTOM)

            F.pack(side=LEFT)

            # Affichage par défaut

            # Méthode de résolution par defaut
            self.methode = res.Euler_explicite

            # Modèle par défaut
            modele = mod.modeles[3*self.diffusion] # correspond à SIR car les modèles sont classés d'abord sans puis avec diffusion

            # Affichage
            self.fig = Figure()

            U = self.methode(N_pop*u_def,modele,[N_pop,beta_def,gamma_def,alpha_def,delta_def,eta_def,tau_def],temps,pas_t) / N_pop
            self.a = self.fig.add_subplot(111)
            self.p1,=self.a.plot(temps, U[:,0], ':', color="blue", label="S")
            self.p3,=self.a.plot(temps, U[:,2], ':', color="red", label="I")
            self.p4,=self.a.plot(temps, U[:,3], ':', color="green", label="R")
            self.a.axis([0-duree/20,duree+duree/20,0-0.05,1+0.05])
            self.a.set_xlabel("Nombre de jours")
            self.a.set_ylabel("Population")
            self.a.legend()

            self.canvas = FigureCanvasTkAgg(self.fig,master=master)
            self.canvas.get_tk_widget().pack()

            # Bouton pour quitter
            quitter = Button(master, text="Quitter ce super programme", command=master.destroy)
            quitter.pack(side=BOTTOM)

            # Bouton pour actualiser
            actualiser = Button(master, text="Actualiser", command=self.affiche)
            actualiser.pack(side=BOTTOM)

        def update_methode(self, i):
            self.methode = res.methodes[i]
            return None

        def affiche(self):
            self.update_methode(int(self.liste_methode.curselection()[0]))

            self.diffusion = self.var1.get()
            temps = np.linspace(0, duree, int(duree/pas_t))

            # Selection du modele
            modele = mod.modeles[self.liste.curselection()[0] + 3*self.diffusion] # car dans modeles il y a d'abord 3 modèles sans diffusion puis 3 avec.

            if self.diffusion == 0:
                self.fig.clf()
                # On récupère les choix de l'utilisateur
                parametres = [N_pop,self.beta.get(),self.gamma.get(),self.alpha.get(),self.delta.get(),self.eta.get(),self.tau.get()]
                u_0 = N_pop*np.array([1-self.I.get(),0,self.I.get(),0,0])
                U = self.methode(u_0,modele,parametres,temps,pas_t) / N_pop
                # On affiche uniquement les courbes utiles
                self.a = self.fig.add_subplot(111)
                self.p1,=self.a.plot(temps, U[:,0], ':', color="blue", label="S")
                if self.liste.curselection()[0]%3==1:
                    self.p2,=self.a.plot(temps, U[:,1], ':', color="orange", label="E")
                self.p3,=self.a.plot(temps, U[:,2], ':', color="red", label="I")
                self.p4,=self.a.plot(temps, U[:,3], ':', color="green", label="R")
                if self.liste.curselection()[0]%3==2:
                    self.p5,=self.a.plot(temps, U[:,4], ':', color="black", label="T")
                self.a.axis([0-duree/20,duree+duree/20,0-0.05,1+0.05])
                self.a.set_xlabel("Nombre de jours")
                self.a.set_ylabel("Population")
                self.a.legend()
            else:
                # On récupère les choix de l'utilisateur
                parametres = [N_pop,self.beta.get(),self.gamma.get(),self.alpha.get(),self.delta.get(),self.eta.get(),self.tau.get(),self.D.get()]
                U = self.methode(u0,modele,parametres,temps,pas_t) / N_pop
                # Selection des temps à afficher
                nb_lignes = 3
                nb_colonnes = 3
                temps_a_afficher = np.array([0,int(len(U)/50),int(len(U)/38),int(len(U)/20),int(len(U)/12),int(len(U)/8),int(len(U)/5),int(len(U)/4),int(len(U)/2)])
                self.fig.clf()
                nb_pts_x = len(u0)
                distance = np.array(range(nb_pts_x))/(nb_pts_x-1)
                for i in range(len(temps_a_afficher)):
                    t = int(temps_a_afficher[i])
                    # On affiche uniquement les courbes utiles
                    self.a = self.fig.add_subplot(nb_lignes,nb_colonnes,i+1,)
                    self.fig.tight_layout()
                    self.a.axis([0-0.05,1+0.05,0-0.05,1+0.05])
                    self.a.set_xlabel("Distance")
                    self.a.set_ylabel("Population")
                    self.a.set_title("t = "+"%.3f"%temps[t])
                    self.p1,=self.a.plot(distance, U[t][:,0], ':', color="blue", label="S")
                    if self.liste.curselection()[0]%3==1:
                        self.p2,=self.a.plot(distance, U[t][:,1], ':', color="orange", label="E")
                    self.p3,=self.a.plot(distance, U[t][:,2], ':', color="red", label="I")
                    self.p4,=self.a.plot(distance, U[t][:,3], ':', color="green", label="R")
                    if self.liste.curselection()[0]%3==2:
                        self.p5,=self.a.plot(distance, U[t][:,4], ':', color="black", label="T")

            self.canvas.draw()

    #Ouverture de la fenetre
    fenetre = Tk()
    fenetre.wm_iconbitmap('enpc_favicon.ico')
    fenetre.title("Modélisation de la propagation d'une épidémie")
    fenetre.geometry("1400x900+30+30")

    boutons = Boutons(fenetre, N_pop, duree, pas_t, u0)
    fenetre.mainloop()
