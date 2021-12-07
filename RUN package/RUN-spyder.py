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
#%%
import numpy as np 
import matplotlib as mpl
from matplotlib import cm
from matplotlib import pyplot as plt 
from scipy import linalg 
import sys, os
mpl.rcParams.update({'font.size': 15})
caseFolder = input("\nEnter address of testCase ...\nExample:  /home/user/folder/Case_1\n") 
try:
    os.mkdir(caseFolder)
except: 
    print("\nFolder already there!!\n")
    flag_EmptyFolder = input("\nEmpty Folder ?? [y/n]\n")
    if(flag_EmptyFolder=='y'):
        filelist = [ f for f in os.listdir(caseFolder) ]
        for f in filelist:
            os.remove(os.path.join(caseFolder, f))
    else:
        sys.exit()
from boundaryDomain import domain, obstacle, obstacleMath
from initialConditions import initialCond
from faceInterpolation import U_face, V_face
from UV_coeff import U_coeffUpwind, V_coeffUpwind
from P_coeff import P_coeffSIMPLE
from correction import correctionSIMPLE
from postProcessing import postProcessor
D, \
xMin, xMax, yMin, yMax, \
dx, dy, nx, ny, \
south_Boundary, north_Boundary, west_Boundary, east_Boundary, \
flag_Outlet =\
domain(caseFolder)
flag_obstacleGeometry = int(input("\n1. Use .txt file as Address \n2. Add elliptical obstacle\n"))
if(flag_obstacleGeometry==1):
    D, Points, Points_x, Points_y = obstacle(D, xMin, xMax, yMin, yMax, nx, ny, dx, dy)   
elif(flag_obstacleGeometry==2):
    D, Points, Points_x, Points_y = obstacleMath(D, xMin, xMax, yMin, yMax, nx, ny, dx, dy)
Lx = xMax-xMin
Ly = yMax-yMin

# Grid - Plot 
xg = np.linspace(xMin, xMax, num=nx+2)
yg = np.linspace(yMax, yMin, num=ny+2)
 
Xg, Yg = np.meshgrid(xg, yg)
fig, ax = plt.subplots(figsize =(Lx*6*1.4, Ly*6))
plt.xlabel('x [m]')
plt.ylabel('y [m]')
cmap = plt.get_cmap('tab10')
ylorbr = cm.get_cmap('tab10', 4)
norm = mpl.colors.Normalize(vmin=-1.5, vmax=2.5)
col = plt.pcolormesh(xg, yg, D[0], cmap = ylorbr, norm=norm)
cbar = plt.colorbar(col)
cbar.set_label('Grid\n-1.0 - Outlet    0.0 - Wall    1.0 - Inlet    2.0 - Wall')
plt.axis('scaled')

# show plot
plt.show()

flag_continue = input("\nContinue or Abort ?? [y/n]\n")
if(flag_continue=='n'):
    sys.exit()


urf_UV = float(input("\nEnter Under-relaxation factor for UV\n"))
urf_P = float(input("\nEnter Under-relaxation factor for P\n"))
rho = float(input("\nEnter Density of fluid [kg/m^3]\n"))
visc = float(input("\nEnter Viscosity of Fluid [Pa*s]\n"))
max_iter = int(input("\nEnter Max Iteration for this SIM\n"))
MaxRes_Mass = float(input("\nEnter Min Mass Imbalance Residual for this SIM\n"))
    
Uxo, Uyo, Po = initialCond(D, nx, ny, flag_Outlet)
Ux_values=[]
Ux_values.append(Uxo)
Uy_values=[]
Uy_values.append(Uyo)
P_values=[]
P_values.append(Po)
Res_U_values=[]
Res_V_values=[]
Res_Mass_values=[]
logRes_Mass=0.0
n=0
iter_ = np.linspace(0, max_iter, max_iter)
temp = np.linspace(0, max_iter, max_iter)
# plt.ion()
fig, ax = plt.subplots(figsize =(Lx*6*1.4, Ly*6))
plt.xlabel('x [m]')
plt.ylabel('y [m]')


while n<max_iter and logRes_Mass>MaxRes_Mass:
    Uue, Uuw, Vun, Vus = U_face(Uxo, Uyo, D, nx, ny)
    Uve, Uvw, Vvn, Vvs = V_face(Uxo, Uyo, D, nx, ny)
    
    
    
    M_Au, M_Bup, Aup, Aue, Auw, Aun, Aus, Bup = \
        U_coeffUpwind(Uxo, Uyo, Po, Uue, Uuw, Vun, Vus, D, rho, visc, urf_UV, nx, ny, dx, dy)
    M_Av, M_Bvp, Avp, Ave, Avw, Avn, Avs, Bvp = \
        V_coeffUpwind(Uxo, Uyo, Po, Uve, Uvw, Vvn, Vvs, D, rho, visc, urf_UV, nx, ny, dx, dy)
    
    
    
    # Momentum Predictor 
    M_Uxs = linalg.solve(M_Au, M_Bup)
    M_Uys = linalg.solve(M_Av, M_Bvp)
    Uxs = np.reshape(M_Uxs, (ny, nx+1))
    Uys = np.reshape(M_Uys, (ny+1, nx))

    
    M_Ap, M_Bpp, App, Ape, Apw, Apn, Aps, Bpp = P_coeffSIMPLE(Aup, Avp, Uxs, Uys, D, nx, ny, dx, dy)
    if(flag_Outlet==0):
        M_Ap[int(ny*nx/4), :] = 0
        M_Ap[int(ny*nx/4), int(ny*nx/4)] = 1
        M_Bpp[int(ny*nx/4)] = 0
    
    # Pressure Correction values
    M_Pc = linalg.solve(M_Ap, M_Bpp)
    
    Pc = np.reshape(M_Pc, (ny, nx))
    
    # Residuals 
    Res_U = M_Au.dot(Uxo.flatten())-M_Bup
    Res_V = M_Av.dot(Uyo.flatten())-M_Bvp
    logRes_U = np.log(np.sum(np.abs(Res_U)))
    logRes_V = np.log(np.sum(np.abs(Res_V)))
    logRes_Mass = np.log(np.sum(np.abs(Bpp)))
    
    # Momentum & Pressure Correction step
    Ux, Uy, P = correctionSIMPLE(Uxs, Uys, Po, Pc, D, Aup, Avp, urf_UV, urf_P, nx, ny, dx, dy)
    
    # Writing Results 
    Ux_values.append(Ux.copy())
    Uy_values.append(Uy.copy())
    P_values.append(P.copy())
    Res_U_values.append(logRes_U)
    Res_V_values.append(logRes_V)
    Res_Mass_values.append(logRes_Mass)
    
    # Prepping for next step 
    Uxo = Ux.copy()
    Uyo = Uy.copy()
    Po = P.copy()
    iter_ = np.linspace(0, n+1, n+1)
    ax.set_xlim(-0.5, n+2)
    max_Res=0
    min_Res=0
    for k in range (0, n+1):
        if(max_Res<Res_U_values[k]):
            max_Res=Res_U_values[k]
        if(max_Res<Res_V_values[k]):
            max_Res=Res_U_values[k]
        if(max_Res<Res_Mass_values[k]):
            max_Res=Res_Mass_values[k]
        if(min_Res>Res_U_values[k]):
            min_Res=Res_U_values[k]
        if(min_Res>Res_V_values[k]):
            min_Res=Res_U_values[k]
        if(min_Res>Res_Mass_values[k]):
            min_Res=Res_Mass_values[k]
    
    ax.set_ylim(min_Res-1, max_Res+1)
    plt.plot(iter_, Res_U_values, color='red', label='Ux-Mom Residual')
    plt.plot(iter_, Res_V_values, color='blue', label='Uy-Mom Residual')
    plt.plot(iter_, Res_Mass_values, color='green', label='Mass Residual')
    plt.legend()
    plt.show()
    n+=1



if(n<max_iter):
    max_iter=n

max_iter
uxq_values = []
uyq_values = []
for n in range (0, max_iter):
    ux1 = Ux_values[n].copy()
    uxq = np.zeros((2*ny+1, 2*nx+1))
    for i in range (0, ny):
        for j in range (0, nx+1):
            uxq[2*i+1, 2*j] = ux1[i, j]
    for j in range (0, nx+1):
        uxq[0, 2*j] = D[1][0, j]
        uxq[-1, 2*j] = D[1][-1, j]
    for i in range (0, ny-1):
        for j in range (0, nx+1):
            uxq[2*i+2, 2*j] = (uxq[2*i+1, 2*j] + uxq[2*i+3, 2*j])/2
    for i in range (0, 2*ny):
        for j in range (0, nx):
            uxq[i, 2*j+1] = (uxq[i, 2*j] + uxq[i, 2*j+2])/2

    
    uy1 = Uy_values[n].copy()
    uyq = np.zeros((2*ny+1, 2*nx+1))

    
    for i in range (0, ny+1):
        for j in range (0, nx):
            uyq[2*i, 2*j+1] = uy1[i, j]
    for i in range (0, ny+1):
        uyq[2*i, 0] = D[2][i, 0]
        uyq[2*i, -1] = D[2][i, -1]
    for i in range (0, ny+1):
        for j in range (0, nx-1):
            uyq[2*i, 2*j+2] = (uyq[2*i, 2*j+1] + uyq[2*i, 2*j+3])/2
    for i in range (0, ny):
        for j in range (0, 2*nx):
            uyq[2*i+1, j] = (uyq[2*i, j] + uyq[2*i+2, j])/2

    uxq_values.append(uxq)
    uyq_values.append(uyq)


#%% 
flag_Validation = input("Validation check for Lid-driven Cavity Re-[100, 1000, 5000] ?? [y/n]")
if(flag_Validation=='y'):
    mpl.rcParams.update({'font.size': 10})
    ghia_u100 = np.array([1, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, -0.13641, -0.20581, -0.2109, -0.15662, -0.1015, -0.06434, -0.04775, -0.04192, -0.03717, 0])
    ghia_u1000 = np.array([1, 0.65928, 0.57492, 0.51117, 0.46604, 0.33304, 0.18719, 0.05702, -0.0608, -0.10648, -0.27805, -0.38289, -0.2973, -0.2222, -0.20196, -0.18109, 0])
    ghia_u5000 = np.array([1, 0.48223, 0.4612, 0.45992, 0.46036, 0.33556, 0.20087, 0.08183, -0.03039, -0.07404, -0.22855, -0.3305, -0.40435, -0.43643, -0.42901, -0.41165, 0])
    ghia_y = np.array([1, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0])
    yq = np.linspace(yMax, yMin, num=2*ny+1)  
    Re = rho*D[1][0, int(nx/2)]*Lx/visc
    if(Re==100):
        plt.scatter(ghia_u100, ghia_y, c='green', label='Ghia etal.')
        plt.title('Reynolds No. 100 Mesh -[{}X{}]'.format(nx, ny))
    if(Re==1000):
        plt.scatter(ghia_u1000, ghia_y, c='green', label='Ghia etal.')
        plt.title('Reynolds No. 1000 Mesh -[{}X{}]'.format(nx, ny))
    if(Re==5000):
        plt.scatter(ghia_u5000, ghia_y, c='green', label='Ghia etal.')
        plt.title('Reynolds No. 5000 Mesh -[{}X{}]'.format(nx, ny))
    plt.plot(uxq_values[max_iter-1][:, nx], yq, label ='Centreline X-Velocity')
    plt.legend()
    plt.show()
flag_Validation = input("Validation check for Fully developed flow ?? [y/n]")
if(flag_Validation=='y'):
    mpl.rcParams.update({'font.size': 10})
    Re = rho*D[1][int(ny/2), 0]*Ly/visc
    yq = np.linspace(yMax, yMin, num=2*ny+1)  
    plt.plot(uxq_values[max_iter-1][:, int(8*(2*nx+1)/10)], yq, label ='Centreline X-Velocity - at X={}'.format(xMin+int(8*(2*nx+1)/10)*dx/2))
    # plt.set_xlabel('Y')
    plt.title('Fully Developed Flow Re = {}'.format(Re))
    plt.legend()
    plt.show()
    xq = np.linspace(xMin, xMax, num=2*nx+1)  
    plt.plot(uxq_values[max_iter-1][ny, :], xq, label ='Centreline X-Velocity - at Y={}'.format(yMax-ny*dy/2))
    # plt.set_xlabel('X')    
    plt.title('Fully Developed Flow Re = {}'.format(Re))
    plt.legend()
    plt.show()
#%%

postProcessor(P_values, uxq_values, uyq_values, \
                  xMin, xMax, yMin, yMax, nx, ny, dx, dy, Points_x, Points_y, \
                  caseFolder, Res_U_values, Res_V_values, Res_Mass_values, \
                  max_iter, min_Res, max_Res)
