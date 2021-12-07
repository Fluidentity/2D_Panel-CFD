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
import matplotlib
from matplotlib import pyplot as plt 
import matplotlib.animation as animation
matplotlib.rcParams.update({'font.size': 15})
def postProcessor(P_values, uxq_values, uyq_values, xMin, xMax, yMin, yMax, nx, ny, dx, dy, Points_x, Points_y, caseFolder, Res_U_values, Res_V_values, Res_Mass_values, max_iter, min_Res, max_Res):
    ('plasma')
    flag_generateResults = int(input("\n1. Generate Results\n2. Abort\n"))
    flag_newMap = 'y'
    flag_newMapOK = 2
    if(flag_generateResults==1):
        flag_generateResidualsPlot = input("Generate Residuals plot ?? [y/n]")
        if(flag_generateResidualsPlot=='y'):
            fig, ax = fig, ax = plt.subplots(figsize =(10, 5))
            iter_ = np.linspace(0, max_iter, max_iter)
            ax.plot(iter_, Res_U_values, color = 'red', label = 'Ux-Mom Residual')
            ax.plot(iter_, Res_V_values, color = 'blue', label = 'Uy-Mom Residual')
            ax.plot(iter_, Res_Mass_values, color = 'green', label = 'Mass Residual')
            ax.set_ylim(min_Res, max_Res)
            ax.set_xlim(0, max_iter)
            ax.set_xlabel('Iterations')
            ax.legend()
            plt.savefig(caseFolder+"/Residuals.jpg")

        flag_generateCSV = int(input("\n[1]- Generate CSV files of Residuals \n[2]- Generate CSV files of U, V and P\n[3]- No CSV Results\n"))
        if(flag_generateCSV==1):
            np.savetxt(caseFolder+"/U_Residual.csv", Res_U_values, delimiter=",")
            np.savetxt(caseFolder+"/V-Residual.csv", Res_V_values, delimiter=",")
            np.savetxt(caseFolder+"/Mass_Residual.csv", Res_Mass_values, delimiter=",")
        elif(flag_generateCSV==2):
            np.savetxt(caseFolder+"/U.csv", uxq_values[max_iter-1], delimiter=",")
            np.savetxt(caseFolder+"/V.csv", uyq_values[max_iter-1], delimiter=",")
            np.savetxt(caseFolder+"/P.csv", P_values[max_iter], delimiter=",")
        mapCounter = 1
        matplotlib.rcParams.update({'font.size': 35})
        while(flag_newMap=='y'):
            print("Generating Map-", mapCounter)
            flag_resultColorMap = input("Add ColorMap ?? [y/n]")
            flag_resultGeometry = input("Add Black Obstacle Geometry ?? [y/n]")
            if(flag_resultColorMap=='y'):
                flag_resultColorMapMarker = int(input("\nSelect Variable\n[1]- Pressure\n[2]- Velocity\n[3]- X-Velocity\n[4]- Y-Velocity\n"))
                flag_resultColorMapContourlines = input("\nWith Labelled Contourlines ?? [y/n]\n")
                if(flag_resultColorMapContourlines=='y'):
                    ContourlinesDensity = int(20/int(input("\nDensity of Contourlines ??\n")))
                min_resultColorMap = float(input("Min value of variable"))
                max_resultColorMap = float(input("Max value of variable"))
            flag_resultStreamlines = input("Add Streamlines ?? [y/n]")
            if(flag_resultStreamlines=='y'):
                flag_resultStreamlinesDensity = int(input("Streamlines Density ??"))
            flag_resultVectors = input("Add Vectors ?? [y/n]")
            
            
            xq = np.linspace(xMin, xMax, num=2*nx+1)
            yq = np.linspace(yMax, yMin, num=2*ny+1)
            xqs = np.linspace(xMin, xMax, num=2*nx+1)
            yqs = np.linspace(yMin, yMax, num=2*ny+1)
            xp = np.linspace(xMin+dx/2, xMax-dx/2, num=nx)
            yp = np.linspace(yMax-dy/2, yMin+dy/2, num=ny)
            Xq, Yq = np.meshgrid(xqs, yqs)
            # creating plot
            Lx = xMax-xMin
            Ly = yMax-yMin
            fig, ax = plt.subplots(figsize =(Lx*10, Ly*10))
            if(flag_resultGeometry=='y'):
                ax.fill(Points_x, Points_y, color='black')
            if(flag_resultColorMap=='y'):
                levels = np.linspace(min_resultColorMap, max_resultColorMap, num=20)
                if(flag_resultColorMapMarker==1):
                    contf = ax.contourf(xp, yp, P_values[max_iter-1], levels=levels, cmap=plt.cm.plasma)
                    if(flag_resultColorMapContourlines=='y'):
                        cont = ax.contour(contf, levels=contf.levels[::ContourlinesDensity], colors='r')
                        # ax.clabel(cont, colors='w', fontsize=14)
                    cbar = fig.colorbar(contf)
                    cbar.set_label('P')
                elif(flag_resultColorMapMarker==2):
                    contf = ax.contourf(xq, yq, np.sqrt(uxq_values[max_iter-1]**2+uyq_values[max_iter-1]**2), levels=levels, cmap=plt.cm.plasma)
                    if(flag_resultColorMapContourlines=='y'):
                        cont = ax.contour(contf, levels=contf.levels[::ContourlinesDensity], colors='r')
                        # ax.clabel(cont, colors='w', fontsize=14)
                    cbar = fig.colorbar(contf)
                    cbar.set_label('Vel')
                elif(flag_resultColorMapMarker==3):
                    contf = ax.contourf(xq, yq, uxq_values[max_iter-1], levels=levels, cmap=plt.cm.plasma)
                    if(flag_resultColorMapContourlines=='y'):
                        cont = ax.contour(contf, levels=contf.levels[::ContourlinesDensity], colors='r')
                        # ax.clabel(cont, colors='w', fontsize=14)
                    cbar = fig.colorbar(contf)
                    cbar.set_label('Vel_x')
                elif(flag_resultColorMapMarker==4):
                    contf = ax.contourf(xq, yq, uyq_values[max_iter-1], levels=levels, cmap=plt.cm.plasma)
                    if(flag_resultColorMapContourlines=='y'):
                        cont = ax.contour(contf, levels=contf.levels[::ContourlinesDensity], colors='r')
                        # ax.clabel(cont, colors='w', fontsize=14)
                    cbar = fig.colorbar(contf)
                    cbar.set_label('Vel_y')
            if(flag_resultStreamlines=='y'):
                ax.streamplot(Xq, Yq, np.flip(uxq_values[max_iter-1], 0) ,np.flip(uyq_values[max_iter-1], 0), density=flag_resultStreamlinesDensity, color='black', linewidth=1)
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])
                ax.axis([xMin, xMax, yMin, yMax])
                ax.set_aspect('equal')
            elif(flag_resultVectors=='y'):
                ax.quiver(Xq, Yq, np.flip(uxq_values[max_iter-1], 0) ,np.flip(uyq_values[max_iter-1], 0))
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])
                ax.axis([xMin, xMax, yMin, yMax])
                ax.set_aspect('equal')
            
            # show plot
            plt.show()
            flag_newMapOK = int(input("\n1. Continue with Map \n2. Generate again \n"))
            if(flag_newMapOK==2):
                continue
            flag_plotall = int(input("\n0. Plot Last frame only [n] - [jpg]\n1. Plot every n iterations - [GIF]\n"))
            if(flag_plotall!=0):
                xq = np.linspace(xMin, xMax, num=2*nx+1)
                yq = np.linspace(yMin, yMax, num=2*ny+1)
                xp = np.linspace(xMin+dx/2, xMax-dx/2, num=nx)
                yp = np.linspace(yMax-dy/2, yMin+dy/2, num=ny)
                Xq, Yq = np.meshgrid(xq, yq)
                Lx = xMax-xMin
                Ly = yMax-yMin
                fig, ax = plt.subplots(figsize =(Lx*10, Ly*10))
                if(flag_resultColorMap=='y'):
                    levels = np.linspace(min_resultColorMap, max_resultColorMap, num=20)
                    if(flag_resultColorMapMarker==1):
                        contf = ax.contourf(xp, yp, P_values[0], levels=levels, cmap=plt.cm.plasma)
                        cbar = fig.colorbar(contf)
                        cbar.set_label('P')
                    elif(flag_resultColorMapMarker==2):
                        contf = ax.contourf(xq, yq, np.sqrt(uxq_values[0]**2+uyq_values[0]**2), levels=levels, cmap=plt.cm.plasma)
                        cbar = fig.colorbar(contf)
                        cbar.set_label('Vel')
                    elif(flag_resultColorMapMarker==3):
                        contf = ax.contourf(xq, yq, uxq_values[0], levels=levels, cmap=plt.cm.plasma)
                        cbar = fig.colorbar(contf)
                        cbar.set_label('Vel_x')
                    elif(flag_resultColorMapMarker==4):
                        contf = ax.contourf(xq, yq, uyq_values[0], levels=levels, cmap=plt.cm.plasma)
                        cbar = fig.colorbar(contf)
                        cbar.set_label('Vel_y')
                def animate(i):
                    ax.clear() 
                    if(flag_resultGeometry=='y'):
                        ax.fill(Points_x, Points_y, color='black')
                    if(flag_resultColorMap=='y'):
                        levels = np.linspace(min_resultColorMap, max_resultColorMap, num=20)
                        if(flag_resultColorMapMarker==1):
                            contf = ax.contourf(xp, yp, P_values[i], levels=levels, cmap=plt.cm.plasma)
                            if(flag_resultColorMapContourlines=='y'):
                                cont = ax.contour(contf, levels=contf.levels[::ContourlinesDensity], colors='r')
                                # ax.clabel(cont, colors='w', fontsize=14)
                            # cbar = fig.colorbar(contf)
                            # cbar.set_label('P')
                        elif(flag_resultColorMapMarker==2):
                            contf = ax.contourf(xq, yq, np.sqrt(uxq_values[i]**2+uyq_values[i]**2), levels=levels, cmap=plt.cm.plasma)
                            if(flag_resultColorMapContourlines=='y'):
                                cont = ax.contour(contf, levels=contf.levels[::ContourlinesDensity], colors='r')
                                # ax.clabel(cont, colors='w', fontsize=14)
                            # cbar = fig.colorbar(contf)
                            # cbar.set_label('Vel')
                        elif(flag_resultColorMapMarker==3):
                            contf = ax.contourf(xq, yq, uxq_values[i], levels=levels, cmap=plt.cm.plasma)
                            if(flag_resultColorMapContourlines=='y'):
                                cont = ax.contour(contf, levels=contf.levels[::ContourlinesDensity], colors='r')
                                # ax.clabel(cont, colors='w', fontsize=14)
                            # cbar = fig.colorbar(contf)
                            # cbar.set_label('Vel_x')
                        elif(flag_resultColorMapMarker==4):
                            contf = ax.contourf(xq, yq, uyq_values[i], levels=levels, cmap=plt.cm.plasma)
                            if(flag_resultColorMapContourlines=='y'):
                                cont = ax.contour(contf, levels=contf.levels[::ContourlinesDensity], colors='r')
                                # ax.clabel(cont, colors='w', fontsize=14)
                            # cbar = fig.colorbar(contf)
                            # cbar.set_label('Vel_y')
                    if(flag_resultStreamlines=='y'):
                        ax.streamplot(Xq, Yq, np.flip(uxq_values[i], 0) ,np.flip(uyq_values[i], 0), density=flag_resultStreamlinesDensity, color='black', linewidth=1)
                    elif(flag_resultVectors=='y'):
                        ax.quiver(Xq, Yq, np.flip(uxq_values[i], 0) ,np.flip(uyq_values[i], 0))
                    
                    ax.set_title('Map-{} Step-[{}]'.format(mapCounter, i))
                interval = 0.1
                ani = animation.FuncAnimation(fig ,animate, max_iter, interval=interval*1e+3,blit=False)
                writergif = animation.PillowWriter(fps=10)
                ani.save(caseFolder+"/Map-{}.gif".format(mapCounter), writer = writergif)
            if(flag_plotall==0):
                fig.savefig(caseFolder+"/Map-{} Step-[{}].jpg".format(mapCounter, max_iter)) 
            flag_newMap = input("Generate another Map ?? [y/n]")
            mapCounter+=1