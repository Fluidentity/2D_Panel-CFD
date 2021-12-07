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

def correctionSIMPLE(Uxs, Uys, Po, Pc, D, Aup, Avp, urf_UV, urf_P, nx, ny, dx, dy):
    P = np.zeros((ny, nx))
    Ux = Uxs.copy()
    Uy = Uys.copy()
    for i in range (0, ny):
        for j in range (0, nx):
            P[i, j] = Po[i, j] + urf_P*Pc[i, j]
            if(D[0][i, j]==2 or D[0][i+1, j]==2):
                Ux[i, j] = Uxs[i, j] + urf_UV*dy*(Pc[i, j-1]-Pc[i, j])/Aup[i, j]
            if(D[0][i, j]==2 or D[0][i, j+1]==2):
                Uy[i, j] = Uys[i, j] + urf_UV*dx*(Pc[i, j]-Pc[i-1, j])/Avp[i, j]             
    for i in range (0, ny):
        for j in range (0, nx+1):
            if(j==nx and D[0][i, j]==-1 and D[0][i+1, j]==-1):
                Ux[i, j] = Ux[i, j-1]

    return Ux, Uy, P
            