# -*- coding: utf-8 -*-
"""
Created on Sat Dec 15 15:35:12 2018

@author: Pedram Tavadze
"""

import pychemia
import numpy as np


def gaussian(x,y,sigma,mu):
    ret = 1/(np.sqrt(2*np.pi*np.linalg.norm(sigma)**2))*np.exp(-1*(((x-mu[0])**2/(2*sigma[0]**2))+(y-mu[1])**2/(2*sigma[1]**2)))
    return ret
    

monolayer_poscar = 'mono-mp-1027525_MoS2.vasp'
supercell_size = (15,15,1)



st = pychemia.code.vasp.read_poscar(monolayer_poscar)
st = st.supercell(supercell_size)
st = st.add_vacuum(10,direction=2)
st.sort_sites()
sigma = (np.max(st.positions[:,:2],axis=0)-np.min(st.positions[:,:2],axis=0))/8
mu = np.average(st.positions[:,:2],axis=0)

height = 100
for iatom in range(st.natom):
    x = st.positions[iatom,0]
    y = st.positions[iatom,1]
    st.positions[iatom,2] += height*gaussian(x,y,sigma,mu)
pychemia.code.vasp.write_poscar(st,filepath="dented_"+monolayer_poscar,direct=False)