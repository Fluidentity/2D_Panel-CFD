#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##########################\/\/\/\/2D_Panel CFD\/\/\/\/##########################
# ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~  ~  ~  #
#   ~   ~   ~   ~   ~   ~   ~______________ ~   ~   ~   ~   ~   ~   ~  ~  ~  ~ #
# ~   ~   ~   ~   ~   ~   ~ /_____________\   ~   ~   ~   ~   ~   ~   ~  ~  ~  #     
#   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~  ~  ~  ~ #
# ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~  ~  ~  #
#   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~   ~  ~  ~  ~ #
##########################\/\/\/\/\/\/\/\/\/\/\/\/\/\/##########################

"""
LICENSE

    This file is a part of 2D_Panel CFD.
    
    2D_Panel CFD is a repository with 2D framework to test new numerical schemes, 
    pressure coupling algorithms, VOF etc. A GUI is intended to be made shortly 
    to make this a user oriented program.
    Copyright (C) <2021>  <Fluidentity>

    2D_Panel CFD is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    2D_Panel CFD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
    
    2D_Panel CFD comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
    
"""

import numpy as np 

def U_face(Uxo, Uyo, D, nx, ny):
    Uue = np.zeros((ny, nx+1))
    Uuw = np.zeros((ny, nx+1))
    Vun = np.zeros((ny, nx+1))
    Vus = np.zeros((ny, nx+1))
    for i in range (0, ny):
        for j in range (0, nx+1):
            if(D[0][i, j]==2 or D[0][i+1, j]==2):
                Uue[i, j] = (Uxo[i, j]+Uxo[i, j+1])/2
                Uuw[i, j] = (Uxo[i, j]+Uxo[i, j-1])/2
                Vun[i, j] = (Uyo[i, j-1]+Uyo[i, j])/2
                Vus[i, j] = (Uyo[i+1, j-1]+Uyo[i+1, j])/2
                
                
    return Uue, Uuw, Vun, Vus
    
    
def V_face(Uxo, Uyo, D, nx, ny):
    Uve = np.zeros((ny+1, nx))
    Uvw = np.zeros((ny+1, nx))
    Vvn = np.zeros((ny+1, nx))
    Vvs = np.zeros((ny+1, nx))
    for i in range (0, ny+1):
        for j in range (0, nx):
            if(D[0][i, j]==2 or D[0][i, j+1]==2):
                Vvs[i, j] = (Uyo[i, j]+Uyo[i+1, j])/2
                Vvn[i, j] = (Uyo[i, j]+Uyo[i-1, j])/2
                Uve[i, j] = (Uxo[i-1, j+1]+Uxo[i, j+1])/2
                Uvw[i, j] = (Uxo[i-1, j]+Uxo[i, j])/2
    return Uve, Uvw, Vvn, Vvs