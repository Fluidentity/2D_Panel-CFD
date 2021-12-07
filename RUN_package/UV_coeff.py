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

def U_coeffUpwind(Uxo, Uyo, Po, Uue, Uuw, Vun, Vus, D, rho, visc, urf_UV, nx, ny, dx, dy):
    Aup = np.zeros((ny, nx+1))
    Aue = np.zeros((ny, nx+1))
    Auw = np.zeros((ny, nx+1))
    Aun = np.zeros((ny, nx+1))
    Aus = np.zeros((ny, nx+1))
    Bup = np.zeros((ny, nx+1))
    vyx = visc*dy/dx
    vxy = visc*dx/dy
    for i in range (0, ny):
        for j in range (0, nx+1):
            if(j!=0 and j!=nx):
                Aup[i, j] = rho*dy*(max(Uue[i, j], 0) + max(-Uuw[i, j], 0)) +\
                            rho*dx*(max(Vun[i, j], 0) + max(-Vus[i, j], 0)) +\
                            2*(vyx+vxy)
                Aue[i, j] = rho*dy*max(-Uue[i, j], 0) + vyx
                Auw[i, j] = rho*dy*max(Uuw[i, j], 0) + vyx
                Aun[i, j] = rho*dx*max(-Vun[i, j], 0) + vxy
                Aus[i, j] = rho*dx*max(Vus[i, j], 0) + vxy
                Bup[i, j] = dy*(Po[i, j-1]-Po[i, j])
            
            # South Boundary Outlet
            if(D[0][i, j] == 2 and D[0][i+1, j] == -1):
                Aup[i, j] += -vxy - rho*dx*max(Vus[i, j], 0)
                Aus[i, j] = 0
            # North Boundary Outlet
            if(D[0][i, j] == -1 and D[0][i+1, j] == 2):
                Aup[i, j] += -vxy - rho*dx*max(-Vun[i, j], 0)
                Aue[i, j] = 0
            # South Boundary Wall
            if(D[0][i, j] == 2 and D[0][i+1, j] == 0):
                Aup[i, j] += 2*vxy
                Aun[i, j] += vxy/3
                Aus[i, j] = 0
                Bup[i, j] += vxy*8*D[1][i+1, j]/3
            # North Boundary Wall
            if(D[0][i, j] == 0 and D[0][i+1, j] == 2):
                Aup[i, j] += 2*vxy
                Aus[i, j] += vxy/3
                Aun[i, j] = 0
                Bup[i, j] += vxy*8*D[1][i, j]/3
            # South Boundary Inlet
            if(D[0][i, j] == 2 and D[0][i+1, j] == 1):
                Aup[i, j] += 2*vxy + 2*rho*dx*max(Vus[i, j], 0)
                Aun[i, j] += vxy/3 + rho*dx*max(Vus[i, j], 0)/3
                Aus[i, j] = 0
                Bup[i, j] += vxy*8*D[1][i+1, j]/3 + 8*rho*dx*D[1][i+1, j]*max(Vus[i, j], 0)/3
            # North Boundary Inlet
            if(D[0][i, j] == 1 and D[0][i+1, j] == 2):
                Aup[i, j] += 2*vxy + 2*rho*dx*max(-Vun[i, j], 0)
                Aus[i, j] += vxy/3 + rho*dx*max(-Vun[i, j], 0)/3
                Aun[i, j] = 0
                Bup[i, j] += vxy*8*D[1][i, j]/3 + 8*rho*dx*D[1][i, j]*max(-Vun[i, j], 0)/3
            

            # Data Point on Boundary 
            if(D[0][i, j]!=2 and D[0][i+1, j]!=2):
                Aup[i, j] = 1
                Aue[i, j] = 0
                Auw[i, j]  = 0
                Aun[i, j]  = 0
                Aus[i, j]  = 0
                Bup[i, j]  = Uxo[i, j]
            # Outlet Boundary 
            if(D[0][i, j]==-1 and D[0][i+1, j]==-1 and j==nx):
                Aup[i, j] = 1
                Aue[i, j] = 0
                Auw[i, j]  = 1
                Aun[i, j]  = 0
                Aus[i, j]  = 0
                Bup[i, j]  = 0
            if(D[0][i, j]==-1 and D[0][i+1, j]==-1 and j==0):
                Aup[i, j] = 1
                Aue[i, j] = 1
                Auw[i, j]  = 0
                Aun[i, j]  = 0
                Aus[i, j]  = 0
                Bup[i, j]  = 0
            # Under-Relaxation Factor
            if(D[0][i, j]==2 or D[0][i+1, j]==2):
                Aup[i, j] = Aup[i, j]/urf_UV
                Bup[i, j] += (1-urf_UV)*Aup[i, j]*Uxo[i, j]

                
    # Matrix Creation 
    M_Bup = Bup.flatten()
    M_Au = np.zeros(((ny)*(nx+1), (ny)*(nx+1)))
    ite=0
    for i in range (0, ny):
        for j in range (0, nx+1):
            M_Au[ite, ite] = Aup[i, j]
            if ((ite+1)%(nx+1)!=0):
                M_Au[ite, ite+1] = -Aue[i, j]
            if((ite%(nx+1)!=0) and (ite!=0)):
                M_Au[ite, ite-1] = -Auw[i, j]
            if (ite<(ny-1)*(nx+1)):
                M_Au[ite, ite+nx+1] = -Aus[i, j]
            if(ite>nx):
                M_Au[ite, ite-nx-1] = -Aun[i, j]
            ite+=1
    
    return M_Au, M_Bup, Aup, Aue, Auw, Aun, Aus, Bup
    
    
def V_coeffUpwind(Uxo, Uyo, Po, Uve, Uvw, Vvn, Vvs, D, rho, visc, urf_UV, nx, ny, dx, dy):
    Avp = np.zeros((ny+1, nx))
    Ave = np.zeros((ny+1, nx))
    Avw = np.zeros((ny+1, nx))
    Avn = np.zeros((ny+1, nx))
    Avs = np.zeros((ny+1, nx))
    Bvp = np.zeros((ny+1, nx))
    vyx = visc*dy/dx
    vxy = visc*dx/dy
    for i in range (0, ny+1):
        for j in range (0, nx):
            if(i!=0 and i!=ny):
                Avp[i, j] = rho*dy*(max(Uve[i, j], 0) + max(-Uvw[i, j], 0)) +\
                            rho*dx*(max(Vvn[i, j], 0) + max(-Vvs[i, j], 0)) +\
                            2*(vxy+vyx)
                Ave[i, j] = rho*dy*max(-Uve[i, j], 0) + vyx
                Avw[i, j] = rho*dy*max(Uvw[i, j], 0) + vyx
                Avn[i, j] = rho*dx*max(-Vvn[i, j], 0) + vxy
                Avs[i, j] = rho*dx*max(Vvs[i, j], 0) + vxy
                Bvp[i, j] = dx*(Po[i, j]-Po[i-1, j])
            # West Boundary Outlet
            if(D[0][i, j] == -1 and D[0][i, j+1] == 2):
                Avp[i, j] += -vyx - rho*dy*max(Uvw[i, j], 0)
                Avw[i, j] = 0
            # East Boundary Outlet
            if(D[0][i, j] == 2 and D[0][i, j+1] == -1):
                Avp[i, j] += -vyx - rho*dy*max(-Uve[i, j], 0)
                Ave[i, j] = 0
            # West Boundary Wall
            if(D[0][i, j] == 0 and D[0][i, j+1] == 2):
                Avp[i, j] += 2*vyx
                Ave[i, j] += vyx/3
                Avw[i, j] = 0
                Bvp[i, j] += vyx*8*D[2][i, j]/3
            # East Boundary Wall
            if(D[0][i, j] == 2 and D[0][i, j+1] == 0):
                Avp[i, j] += 2*vyx
                Avw[i, j] += vyx/3
                Ave[i, j] = 0
                Bvp[i, j] += vyx*8*D[2][i, j+1]/3
            # West Boundary Inlet
            if(D[0][i, j] == 1 and D[0][i, j+1] == 2):
                Avp[i, j] += 2*vyx + 2*rho*dy*max(Uvw[i, j], 0)
                Ave[i, j] += vyx/3 + rho*dy*max(Uvw[i, j], 0)/3
                Avw[i, j] = 0
                Bvp[i, j] += vyx*8*D[2][i, j]/3 + 8*rho*dy*D[2][i, j]*max(Uvw[i, j], 0)/3
            # East Boundary Inlet
            if(D[0][i, j] == 2 and D[0][i, j+1] == 1):
                Avp[i, j] += 2*vyx + 2*rho*dy*max(-Uve[i, j], 0)
                Avw[i, j] += vyx/3 + rho*dy*max(-Uve[i, j], 0)/3
                Ave[i, j] = 0
                Bvp[i, j] += vyx*8*D[2][i, j+1]/3 + 8*rho*dy*D[2][i, j+1]*max(-Uve[i, j], 0)/3
            

            # Data Point on Boundary 
            if(D[0][i, j]!=2 and D[0][i, j+1]!=2):
                Avp[i, j] = 1
                Ave[i, j] = 0
                Avw[i, j]  = 0
                Avn[i, j]  = 0
                Avs[i, j]  = 0
                Bvp[i, j]  = Uyo[i, j]
            # Outlet Boundary 
            if(D[0][i, j]==-1 and D[0][i, j+1]==-1 and i==0):
                Avp[i, j] = 1
                Ave[i, j] = 0
                Avw[i, j]  = 0
                Avn[i, j]  = 0
                Avs[i, j]  = 1
                Bvp[i, j]  = 0
            if(D[0][i, j]==-1 and D[0][i, j+1]==-1 and i==ny):
                Avp[i, j] = 1
                Ave[i, j] = 0
                Avw[i, j]  = 0
                Avn[i, j]  = 1
                Avs[i, j]  = 0
                Bvp[i, j]  = 0
            # Under-Relaxation Factor
            if(D[0][i, j]==2 or D[0][i, j+1]==2):
                Avp[i, j] = Avp[i, j]/urf_UV
                Bvp[i, j] += (1-urf_UV)*Avp[i, j]*Uyo[i, j]
    # Matrix Creation
    M_Bvp = Bvp.flatten()
    M_Av = np.zeros(((ny+1)*(nx), (ny+1)*(nx)))
    ite=0
    for i in range (0, ny+1):
        for j in range (0, nx):
            M_Av[ite, ite] = Avp[i, j]
            if ((ite+1)%(nx)!=0):
                M_Av[ite, ite+1] = -Ave[i, j]
            if((ite%(nx)!=0) and (ite!=0)):
                M_Av[ite, ite-1] = -Avw[i, j]
            if (ite<(ny)*(nx)):
                M_Av[ite, ite+nx] = -Avs[i, j]
            if(ite>nx-1):
                M_Av[ite, ite-nx] = -Avn[i, j]
            ite+=1
    
    return M_Av, M_Bvp, Avp, Ave, Avw, Avn, Avs, Bvp