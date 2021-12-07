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

def P_coeffSIMPLE(Aup, Avp, Uxs, Uys, D, nx, ny, dx, dy):
    App =np.zeros((ny, nx))
    Ape =np.zeros((ny, nx))
    Apw =np.zeros((ny, nx))
    Apn =np.zeros((ny, nx))
    Aps =np.zeros((ny, nx))
    Bpp =np.zeros((ny, nx))
    for j in range (0, nx):
        for i in range (0, ny):
            if (D[0][i, j+1]!=2 and D[0][i+1, j+1]!=2):
                Ape[i, j] = 0
            else:
                Ape[i, j] = dy**2/Aup[i, j+1]
            if (D[0][i, j]!=2 and D[0][i+1, j]!=2):
                Apw[i, j] = 0
            else:
                Apw[i, j] = dy**2/Aup[i, j]
            if (D[0][i, j]!=2 and D[0][i, j+1]!=2):
                Apn[i, j] = 0
            else:
                Apn[i, j] = dx**2/Avp[i, j]
            if (D[0][i+1, j]!=2 and D[0][i+1, j+1]!=2):
                Aps[i, j] = 0
            else:
                Aps[i, j] = dx**2/Avp[i+1, j]
            App[i, j] = Ape[i, j] + Apw[i, j] + Apn[i, j] + Aps[i, j]
            Bpp[i, j] = dy*(Uxs[i, j]-Uxs[i, j+1]) + \
                        dx*(Uys[i+1, j]-Uys[i, j])
            if(D[0][i, j]==-1 or D[0][i+1, j]==-1 or D[0][i, j+1]==-1 or D[0][i+1, j+1]==-1):
                App[i, j] = 1
                Ape[i, j] = 0
                Apw[i, j] = 0
                Apn[i, j] = 0
                Aps[i, j] = 0
                Bpp[i, j] = 0
            if(D[0][i, j]!=2 and D[0][i+1, j]!=2 and D[0][i, j+1]!=2 and D[0][i+1, j+1]!=2):
                App[i, j] = 1
                Ape[i, j] = 0
                Apw[i, j] = 0
                Apn[i, j] = 0
                Aps[i, j] = 0
                Bpp[i, j] = 0
                
            
    
    M_Bpp = Bpp.flatten()
    
    M_A = np.zeros((ny*nx, ny*nx))
    ite=0
    for i in range (0, ny):
        for j in range (0, nx):
            M_A[ite, ite] = App[i, j]
            if ((ite+1)%(nx)!=0):
                M_A[ite, ite+1] = -Ape[i, j]
            if(ite%(nx)!=0):
                M_A[ite, ite-1] = -Apw[i, j]
            if (ite<(ny-1)*(nx)):
                M_A[ite, ite+nx] = -Aps[i, j]
            if(ite>nx-1):
                M_A[ite, ite-nx] = -Apn[i, j]
            ite+=1

    
    return M_A, M_Bpp, App, Ape, Apw, Apn, Aps, Bpp