from tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np
import matplotlib.pyplot as plt

import os
os.chdir("C:/Users/Anthony/Desktop/COV_A5/")
import methodes_resolution as res
import modelisation as mod

# Paramètres par défaut
alpha_def = 0.2
beta_def = 1.4
gamma_def = 0.1
D_def = 0.01
S_def = 0.9
E_def = 0.
I_def = 0.1
R_def = 0.
u_def = np.array([S_def,E_def,I_def,R_def])

modeles = mod.modeles # Les modèles sont classés d'abord sans diffusion puis ensuite avec.

def open(methode, N_pop, duree, pas_t, temps, diffusion, u0):
    class Boutons:
        def __init__(self, master, methode, N_pop, duree, pas_t, temps, diffusion, u0):
            # Titre fenetre
            champ_label = Label(master, anchor='n', text="Modélisation de l'évolution d'une population de " + str(N_pop) + " personnes sur " + str(duree) + " jours " + diffusion*"avec" + (1-diffusion)*"sans" + " diffusion")
            champ_label.pack()

            # Choix du modèle
            self.liste = Listbox(master, height=2, width=25)
            self.liste.pack(side=RIGHT)
            self.liste.insert(0, "Modèle SIR")
            self.liste.insert(1, "Modèle SEIR")
            self.liste.select_set(0)

            # Choix des paramètres d'évolution
            F = Frame(master)
            F1 = LabelFrame(F, text="Paramètres d'évolution")

            self.alpha=Scale(F1, orient="horizontal", from_=0, to=1, resolution=0.05, label="Inverse de la duree d'incubation", tickinterval=0.25, length=250)
            self.alpha.set(alpha_def)
            self.alpha.pack(side=TOP)

            self.beta=Scale(F1, orient="horizontal", from_=0, to=4, resolution=0.1, label="Contacts par pers. infectee par jour", tickinterval=1., length=250)
            self.beta.set(beta_def)
            self.beta.pack()

            self.gamma=Scale(F1, orient="horizontal", from_=0, to=1, resolution=0.05, label="Inverse de la duree de guerison", tickinterval=0.25, length=250)
            self.gamma.set(gamma_def)
            self.gamma.pack(side=BOTTOM)

            # Choix du paramètre de diffusion
            if diffusion == 1:
                self.D=Scale(F1, orient="horizontal", from_=0, to=10, resolution=0.1, label="Parametre de diffusion", tickinterval=2., length=250)
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
            modele = modeles[2*diffusion] # car dans modeles il y a d'abord 3 modèles sans diffusion puis 2 avec.

            # On affiche une première version par défaut
            if diffusion == 0:
                U = methode(N_pop*u_def,modele,[N_pop,beta_def,gamma_def,alpha_def],temps,pas_t) / N_pop
                S = U[:,0]
                E = U[:,1]
                I = U[:,2]
                R = U[:,3]
            else:
                U = methode(u0,modele,[N_pop,beta_def,gamma_def,alpha_def,D_def],temps,pas_t) / N_pop
                S = U[0][:,0]
                E = U[0][:,1]
                I = U[0][:,2]
                R = U[0][:,3]

            fig = Figure()
            self.a = fig.add_subplot(111)
            if diffusion == 0:
                self.p1,=self.a.plot(temps,S,':',color="blue",label="S")
                self.p2,=self.a.plot(temps,E,':',color="orange",label="E")
                self.p3,=self.a.plot(temps,I,':',color="red",label="I")
                self.p4,=self.a.plot(temps,R,':',color="green",label="R")
                self.a.axis([0-duree/20,duree+duree/20,0-0.05,1+0.05])
                self.a.set_xlabel("Nombre de jours")
            else:
                nb_pts_x = len(u0)
                distance = np.array(range(nb_pts_x))/(nb_pts_x-1)
                self.p1,=self.a.plot(distance,S,':',color="blue",label="S")
                self.p2,=self.a.plot(distance,E,':',color="orange",label="E")
                self.p3,=self.a.plot(distance,I,':',color="red",label="I")
                self.p4,=self.a.plot(distance,R,':',color="green",label="R")
                self.a.axis([0-0.05,1+0.05,0-0.05,1+0.05])
                self.a.set_xlabel("Nombre de jours")
            self.a.legend()
            self.a.set_ylabel("Population")
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
            # Selection du modele
            modele = modeles[self.liste.curselection()[0] + 2*diffusion] # car dans modeles il y a d'abord 3 modèles sans diffusion puis 2 avec.

            if diffusion == 0:
                # On récupère les choix de l'utilisateur
                parametres = [N_pop,self.beta.get(),self.gamma.get(),self.alpha.get()]
                u_0 = N_pop*np.array([1-self.I.get(),0,self.I.get(),0])

                U = methode(u_0,modele,parametres,temps,pas_t) / N_pop
                S = U[:,0]
                E = U[:,1]
                I = U[:,2]
                R = U[:,3]

                self.p1.set_ydata(S)
                self.p2.set_ydata(E)
                self.p3.set_ydata(I)
                self.p4.set_ydata(R)
            else:
                # On récupère les choix de l'utilisateur
                parametres = [N_pop,self.beta.get(),self.gamma.get(),self.alpha.get(),self.D.get()]

                U = methode(u0,modele,parametres,temps,pas_t) / N_pop
                # A MODIFIER : IL FAUT AFFICHER A CHAQUE INSTANT EN UTILISANT PLT.PAUSE()
                # pour obtenir S pour chaque distance à l'instant t on utilise U[t][:,0]
                for t in range(len(temps)):
                    S = U[t][:,0]
                    E = U[t][:,1]
                    I = U[t][:,2]
                    R = U[t][:,3]

                    self.p1.set_ydata(S)
                    self.p2.set_ydata(E)
                    self.p3.set_ydata(I)
                    self.p4.set_ydata(R)
                    #plt.pause(0.1)
            self.canvas.draw()

    #Ouverture de la fenetre
    fenetre = Tk()
    fenetre.title("Modélisation de la propagation d'une épidémie")
    fenetre.geometry("1300x600+30+30")

    boutons = Boutons(fenetre, methode, N_pop, duree, pas_t, temps, diffusion, u0)
    fenetre.mainloop()
