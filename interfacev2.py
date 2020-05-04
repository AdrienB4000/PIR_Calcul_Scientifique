from tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

import os
os.chdir("C:/Users/belfe/OneDrive/Documents/Cours/Ponts_1A_2019-2020/PIR/")
import propagation as prop

#Paramètres par défaut
alpha_def=0.2
beta_def=1.4
gamma_def=0.1
S_def=0.9
E_def=0.
I_def=0.1
R_def=0.
Nb=1000000

duree = 80 #jours
nb_pts = 10000
pas = duree/nb_pts
temps = np.linspace(0,duree,nb_pts)

liste_modeles=[prop.f_SIR,prop.f_SEIR,prop.f_SIR_dyn,prop.f_SIS]

class Boutons:
    def __init__(self, master):
        #Titre fenetre
        champ_label = Label(master, anchor=N, text="Modélisation")
        champ_label.pack()

        #Choix du modèle
        self.liste = Listbox(master, height=4, width=25)
        self.liste.pack(side=RIGHT)
        self.liste.insert(0, "Modèle SIR")
        self.liste.insert(1, "Modèle SEIR")
        self.liste.insert(2, "Modèle SIR avec diffusion")
        self.liste.insert(3, "Modèle SEIR avec diffusion")
        self.liste.select_set(0)

        #Choix des paramètres d'évolution
        F=Frame(master)
        F1=LabelFrame(F, text="Paramètres d'évolution")

        self.alpha=Scale(F1, orient='horizontal', from_=0, to=1, resolution=0.05, label="Taux d'incubation", tickinterval=0.25, length=150)
        self.alpha.set(alpha_def)
        self.alpha.pack(side=TOP)

        self.beta=Scale(F1, orient='horizontal', from_=0, to=4, resolution=0.1, label='Taux de contact', tickinterval=0.5, length=150)
        self.beta.set(beta_def)
        self.beta.pack()

        self.gamma=Scale(F1, orient='horizontal', from_=0, to=1, resolution=0.05, label='Taux de guérison', tickinterval=0.25, length=150)
        self.gamma.set(gamma_def)
        self.gamma.pack(side=BOTTOM)

        F1.pack(side=TOP)

        #Choix des paramètres initiaux
        F2=LabelFrame(F, text="Paramètres initiaux")

        self.S=Scale(F2, orient='horizontal', from_=0, to=1, resolution=0.01, label="Taux de susceptibles", tickinterval=0.25, length=150, command=self.commande1)
        self.S.set(S_def)
        self.S.pack(side=TOP)

        self.E=Scale(F2, orient='horizontal', from_=0, to=1, resolution=0.01, label="Taux d'exposés", tickinterval=0.25, length=150, command=self.commande1)
        self.E.set(E_def)
        self.E.pack()

        self.I=Scale(F2, orient='horizontal', from_=0, to=1-self.S.get(), resolution=0.01, label="Taux d'infectés", tickinterval=0.25, length=150,command=self.commande1)
        self.I.set(I_def)
        self.I.pack(side=BOTTOM)

        self.R=Scale(F2, orient='horizontal', from_=0, to=1-self.S.get(), resolution=0.01, label="Taux de guéris", tickinterval=0.25, length=150,command=self.commande1)
        self.R.set(R_def)
        self.R.pack(side=BOTTOM)

        F2.pack(side=BOTTOM)
        F.pack(side=LEFT)

        #Plot initial
        parametres=[Nb,self.beta.get(),self.gamma.get(),0,0,0,self.alpha.get(),0,0]
        methode=prop.Runge_Kutta_4
        modele=liste_modeles[(self.liste.curselection())[0]]
        u0=[0,Nb*self.S.get(),Nb*self.E.get(),Nb*self.I.get(),Nb*self.R.get()]
        print(self.beta.get())
        cree=True

        U=prop.methode(u0,modele,parametres)
        Maternal=U[:,0]
        Susceptible=U[:,1]
        Exposed=U[:,2]
        Infected=U[:,3]
        Recovered=U[:,4]

        fig = Figure()
        self.a = fig.add_subplot(111)
        self.p1,=self.a.plot(temps,Susceptible,':',color="blue",label="S")
        self.p2,=self.a.plot(temps,Infected,':',color="red",label="I")
        self.p3,=self.a.plot(temps,Recovered,':',color="green",label="R")
        self.p4,=self.a.plot(temps,Exposed,':',color="orange",label="E")
        self.a.legend()
        self.canvas = FigureCanvasTkAgg(fig,master=master)
        self.canvas.get_tk_widget().pack()

        #Un bouton pour actualiser
        b=Button(master, text="Actualiser", command=self.affiche)
        b.pack(side=BOTTOM)

    def change_limites_E(self,val):
        self.S.configure(to=1-self.E.get()-self.I.get()-self.R.get(), tickinterval=(1-self.E.get()-self.I.get()-self.R.get())/4)
        self.E.configure(to=1-self.S.get()-self.I.get()-self.R.get(),tickinterval=(1-self.S.get()-self.I.get()-self.R.get())/4)
        self.I.configure(to=1-self.E.get()-self.S.get()-self.R.get(),tickinterval=(1-self.E.get()-self.S.get()-self.R.get())/4)
        self.R.configure(to=1-self.E.get()-self.I.get()-self.S.get(),tickinterval=(1-self.E.get()-self.I.get()-self.S.get())/4)

    def change_limites_noE(self,val):
        self.S.configure(to=1-self.I.get()-self.R.get(), tickinterval=(1-self.I.get()-self.R.get())/4)
        self.I.configure(to=1-self.S.get()-self.R.get(),tickinterval=(1-self.S.get()-self.R.get())/4)
        self.R.configure(to=1-self.I.get()-self.S.get(),tickinterval=(1-self.I.get()-self.S.get())/4)

    def commande1(self,val):
        if ((self.liste.curselection())[0]==1 or (self.liste.curselection())[0]==3):
            self.change_limites_E(val)
        else:
            self.change_limites_noE(val)

    def affiche(self):
        parametres=[Nb,self.beta.get(),self.gamma.get(),0,0,0,self.alpha.get(),0,0]
        methode=prop.Runge_Kutta_4
        modele=liste_modeles[(self.liste.curselection())[0]]
        u0=[0,Nb*self.S.get(),Nb*self.E.get(),Nb*self.I.get(),Nb*self.R.get()]
        cree=True

        U=prop.methode(u0,modele,parametres)
        Maternal=U[:,0]
        Susceptible=U[:,1]
        Exposed=U[:,2]
        Infected=U[:,3]
        Recovered=U[:,4]

        self.p1.set_ydata(Susceptible)
        self.p2.set_ydata(Infected)
        self.p3.set_ydata(Recovered)
        self.p4.set_ydata(Exposed)
        self.a.legend()
        self.canvas.draw()

#Ouverture de la fenetre
fenetre = Tk()
fenetre.title("Modélisation de la propagation d'une épidémie")
fenetre.geometry("1300x600+30+30")

boutons=Boutons(fenetre)
fenetre.mainloop()
