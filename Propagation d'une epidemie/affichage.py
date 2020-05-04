from tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np

import os
os.chdir("C:/Users/Anthony/Desktop/COV_A5/")
import methodes_resolution as res
import modelisation as mod

# Paramètres par défaut
alpha_def = 0.2
beta_def = 1.4
gamma_def = 0.1
S_def = 0.9
E_def = 0.
I_def = 0.1
R_def = 0.
u_def = np.array([S_def,E_def,I_def,R_def])

modeles = mod.modeles

def open(methode,N,duree,pas,temps,Lambda,mu,D):
    class Boutons:
        def __init__(self, master, methode, N0, duree, pas, temps, Lambda, mu, D):
            # Titre fenetre
            champ_label = Label(master, anchor='n', text="Modélisation de l'évolution d'une population de " + str(N0) + " personnes sur " + str(duree) + " jours")
            champ_label.pack()

            # Choix du modèle
            self.liste = Listbox(master, height=4, width=25)
            self.liste.pack(side=RIGHT)
            self.liste.insert(0, "Modèle SIR")
            self.liste.insert(1, "Modèle SIR avec dynamique")
            self.liste.insert(2, "Modèle SIR avec diffusion")
            self.liste.insert(3, "Modèle SEIR sans diffusion")
            self.liste.select_set(0)

            # Choix des paramètres d'évolution
            F = Frame(master)
            F1 = LabelFrame(F, text="Paramètres d'évolution")

            self.alpha=Scale(F1, orient="horizontal", from_=0, to=1, resolution=0.05, label="Inverse de la duree d'incubation", tickinterval=0.25, length=150)
            self.alpha.set(alpha_def)
            self.alpha.pack(side=TOP)

            self.beta=Scale(F1, orient="horizontal", from_=0, to=4, resolution=0.1, label="Contacts par pers. infectee par jour", tickinterval=0.5, length=150)
            self.beta.set(beta_def)
            self.beta.pack()

            self.gamma=Scale(F1, orient="horizontal", from_=0, to=1, resolution=0.05, label="Inverse de la duree de guerison", tickinterval=0.25, length=150)
            self.gamma.set(gamma_def)
            self.gamma.pack(side=BOTTOM)

            F1.pack(side=TOP)

            # Choix des paramètres initiaux
            F2 = LabelFrame(F, text="Paramètres initiaux")

            self.S=Scale(F2, orient='horizontal', from_=0, to=1, resolution=0.01, label="Taux de susceptibles", tickinterval=0.25, length=150, command=self.commande1)
            self.S.set(S_def)
            self.S.pack(side=TOP)

            self.E=Scale(F2, orient='horizontal', from_=0, to=1, resolution=0.01, label="Taux d'exposés", tickinterval=0.25, length=150, command=self.commande1)
            self.E.set(E_def)
            self.E.pack()

            self.I=Scale(F2, orient='horizontal', from_=0, to=1, resolution=0.01, label="Taux d'infectés", tickinterval=0.25, length=150,command=self.commande1)
            self.I.set(I_def)
            self.I.pack(side=BOTTOM)

            self.R=Scale(F2, orient='horizontal', from_=0, to=1, resolution=0.01, label="Taux de guéris", tickinterval=0.25, length=150,command=self.commande1)
            self.R.set(R_def)
            self.R.pack(side=BOTTOM)

            F2.pack(side=BOTTOM)
            F.pack(side=LEFT)

            # Plot initial
            methode = res.Euler_explicite
            modele = mod.f_SIR

            U = methode(u_def,modele,[N0,beta_def,gamma_def,Lambda,mu,alpha_def,D],temps,pas)
            S = U[:,0]
            E = U[:,1]
            I = U[:,2]
            R = U[:,3]

            fig = Figure()
            self.a = fig.add_subplot(111)
            self.p1,=self.a.plot(temps,S,':',color="blue",label="S")
            self.p2,=self.a.plot(temps,E,':',color="orange",label="E")
            self.p3,=self.a.plot(temps,I,':',color="red",label="I")
            self.p4,=self.a.plot(temps,R,':',color="green",label="R")
            self.a.legend()
            self.a.axis([0,duree,0,1])
            self.a.set_xlabel("Nombre de jours")
            self.a.set_ylabel("Population")
            self.canvas = FigureCanvasTkAgg(fig,master=master)
            self.canvas.get_tk_widget().pack()

            # Bouton pour quitter
            quitter = Button(master, text = "Quitter ce super programme", command=master.destroy)
            quitter.pack(side=BOTTOM)

            # Bouton pour actualiser
            actualiser = Button(master, text="Actualiser", command=self.affiche)
            actualiser.pack(side=BOTTOM)

            # Ajout d'un design de qualité
            #master.iconbitmap("enpc_favicon.ico")

        def change_limites_E(self,val):
            self.S.configure(to=1-self.E.get()-self.I.get()-self.R.get(),tickinterval=(1-self.E.get()-self.I.get()-self.R.get())/4)
            self.E.configure(to=1-self.S.get()-self.I.get()-self.R.get(),tickinterval=(1-self.S.get()-self.I.get()-self.R.get())/4)
            self.I.configure(to=1-self.E.get()-self.S.get()-self.R.get(),tickinterval=(1-self.E.get()-self.S.get()-self.R.get())/4)
            self.R.configure(to=1-self.E.get()-self.I.get()-self.S.get(),tickinterval=(1-self.E.get()-self.I.get()-self.S.get())/4)

        def change_limites_noE(self,val):
            self.S.configure(to=1-self.I.get()-self.R.get(),tickinterval=(1-self.I.get()-self.R.get())/4)
            self.I.configure(to=1-self.S.get()-self.R.get(),tickinterval=(1-self.S.get()-self.R.get())/4)
            self.R.configure(to=1-self.I.get()-self.S.get(),tickinterval=(1-self.I.get()-self.S.get())/4)

        def commande1(self,val):
            if (self.liste.curselection()[0]>2):
                self.change_limites_E(val)
            else:
                self.change_limites_noE(val)

        def affiche(self):
            parametres = [N,self.beta.get(),self.gamma.get(),Lambda,mu,self.alpha.get(),D]
            modele = modeles[self.liste.curselection()[0]]
            u0 = N*np.array([self.S.get(),self.E.get(),self.I.get(),self.R.get()])
            cree = True

            U = methode(u0,modele,parametres,temps,pas) / N
            S = U[:,0]
            E = U[:,1]
            I = U[:,2]
            R = U[:,3]

            self.p1.set_ydata(S)
            self.p2.set_ydata(E)
            self.p3.set_ydata(I)
            self.p4.set_ydata(R)
            self.a.legend()
            self.a.axis([0,duree,0,1])
            self.a.set_xlabel("Nombre de jours")
            self.a.set_ylabel("Population")
            self.canvas.draw()

    #Ouverture de la fenetre
    fenetre = Tk()
    fenetre.title("Modélisation de la propagation d'une épidémie")
    fenetre.geometry("1300x600+30+30")

    boutons = Boutons(fenetre,methode,N,duree,pas,temps,Lambda,mu,D)
    fenetre.mainloop()
