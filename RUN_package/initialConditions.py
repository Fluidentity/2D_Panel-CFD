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

def initialCond(D, nx, ny, flag_Outlet):
    Uxo = np.zeros((ny, nx+1))
    Uyo = np.zeros((ny+1, nx))
    Po = np.zeros((ny, nx))
        
        
    for i in range (0, ny):
        for j in range (0, nx+1):
            if(i==0):
                if((D[0][i, j]==0 and D[0][i+1, j]==0) or (D[0][i, j]==1 and D[0][i+1, j]==1)):
                    Uxo[i, j] = D[1][i+1, j]
            else:
                if((D[0][i, j]==0 and D[0][i+1, j]==0) or (D[0][i, j]==1 and D[0][i+1, j]==1)):
                    Uxo[i, j] = D[1][i, j]

    
                
    for i in range (0, ny+1):
        for j in range (0, nx):
            if(j==0):
                if((D[0][i, j]==0 and D[0][i, j+1]==0) or (D[0][i, j]==1 and D[0][i, j+1]==1)):
                    Uyo[i, j] = D[2][i, j+1]
            else:
                if((D[0][i, j]==0 and D[0][i, j+1]==0) or (D[0][i, j]==1 and D[0][i, j+1]==1)):
                    Uyo[i, j] = D[2][i, j]  
    if(flag_Outlet==1):
        for i in range (0, ny):
            for j in range (0, nx):
                if(D[0][i, j]==-1):
                    Po[i, j] = D[3][i, j]
                if(D[0][i+1, j]==-1):
                    Po[i, j] = D[3][i+1, j]
                if(D[0][i, j+1]==-1):
                    Po[i, j] = D[3][i, j+1]
                if(D[0][i+1, j+1]==-1):
                    Po[i, j] = D[3][i+1, j+1]
                
    return Uxo, Uyo, Po