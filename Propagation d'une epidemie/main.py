import numpy as np
from math import *

import os
import affichage as aff

## Parametres

# Parametres d'affichage
duree = 50 # Duree (jour)
nb_pts_t = 10000
pas_t = duree/nb_pts_t

# Parametres de modelistaion
N_pop = 10000 # Population totale
nb_pts_x = 9

# Choix de u0 pour le cas avec diffusion
x = np.linspace(0, 1, nb_pts_x)
sigma = 0.2
I0 = N_pop/2/pi/sigma*np.exp(-((x-1/2)/sigma)**2)
S0 = N_pop*np.ones(nb_pts_x) - I0
u0 = np.zeros((nb_pts_x,5))
u0[:,0] = S0
u0[:,2] = I0

aff.open(N_pop, duree, pas_t, u0)
