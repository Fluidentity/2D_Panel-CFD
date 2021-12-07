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
from matplotlib import pyplot as plt
import os

def domain(caseFolder):
    # Domain input 
    xMin = float(input("\nEnter X_min value\n"))
    xMax = float(input("\nEnter X_max value\n"))
    dx_temp = float(input("\nEnter Cell size in X-axis [dx]\n"))
    nx = int((xMax-xMin)/dx_temp)
    dx = (xMax-xMin)/nx
    print("\n--- No. of elements in X-axis [{}] ---\n".format(nx))
    yMin = float(input("\nEnter Y_min value\n"))
    yMax = float(input("\nEnter Y_max value\n"))
    dy_temp = float(input("\nEnter Cell size in Y-axis [dy]\n"))
    ny = int((yMax-yMin)/dy_temp)
    dy = (yMax-yMin)/ny
    print("\n--- No. of elements in Y-axis [{}] ---\n".format(ny))
    
    temp = "/boundary.txt"
    boundary_info = open(caseFolder+temp, 'w')
    # Boundary revealing Matrix 
    B1 = 2*np.ones((ny+1, nx+1))
    
    # Scalar Matrix [U, V, P]
    B2 = -1000.5*np.ones((ny+1, nx+1))
    
    # Domain Matrix 
    D = []
    D.append(B1)
    D.append(B2)
    D.append(B2.copy())
    D.append(B2.copy())
    
    # Boundary Definition
    flag_Boundary = 0
    flag_SouthBoundary = 0
    flag_NorthBoundary = 0
    flag_WestBoundary = 0
    flag_EastBoundary = 0
    flag_Outlet = 0
    
    # Boundary Defining Loop
    while (flag_Boundary < 1):
        
        # Southern Boundary Defining Loop
        
        while (flag_SouthBoundary < 1):
            
            flag_Inlet_SouthBoundary = 0
            flag_Outlet_SouthBoundary = 0
            flag_SlidingWall_SouthBoundary = 0
            flag = input("\nInlet in South Boundary ?? [y/n]\n")
            if(flag=='y'):
                flag_Inlet_SouthBoundary = 1
            flag = input("\nSliding Wall in South Boundary ?? [y/n]\n")
            if(flag=='y'):
                flag_SlidingWall_SouthBoundary = 1
            flag = input("\nOutlet in South Boundary ?? [y/n]\n")
            if(flag=='y'):
                flag_Outlet_SouthBoundary = 1         
            
            # Initialising Boundary Markers
            south_Boundary = np.zeros((3, 6))
            south_Boundary[1, 0] = 1
            south_Boundary[2, 0] = -1
            # Temporary Wall
            south_Boundary[0, 1] = xMin
            south_Boundary[0, 2] = xMax
            D[0][-1, :] = 0
            D[1][-1, :] = 0
            D[2][-1, :] = 0
            # South Boundary Sliding U-velocity
            if(flag_SlidingWall_SouthBoundary==1):
                south_Boundary[0, 3] = float(input("\nEnter U-Velocity for Sliding Wall in South Boundary\n"))
                D[1][-1, 1:nx] = south_Boundary[0, 3]
            # Inlet definition in South Boundary
            if(flag_Inlet_SouthBoundary==1):
                south_Boundary[1, 1] = float(input("\nEnter xmin for Inlet in South Boundary\n"))
                southInlet_min_index = int((south_Boundary[1, 1]-xMin)/dx)
                if(southInlet_min_index==0):
                    southInlet_min_index+=1
                south_Boundary[1, 2] = float(input("\nEnter xmax for Inlet in South Boundary\n"))
                southInlet_max_index = int((south_Boundary[1, 2]-xMin)/dx)
                if(southInlet_max_index==nx):
                    southInlet_max_index-=1
                south_Boundary[1, 4] = float(input("\nEnter V-Velocity for Inlet in South Boundary\n"))
                # Boundary Marker
                D[0][-1, southInlet_min_index:southInlet_max_index+1] = 1.0
                # U-vel - Initialised to 0
                D[1][-1, southInlet_min_index:southInlet_max_index+1] = 0
                # V-vel - User Initialised
                D[2][-1, southInlet_min_index:southInlet_max_index+1] = south_Boundary[1, 4]
                
            # Outlet definition in South Boundary
            if(flag_Outlet_SouthBoundary==1):
                flag_Outlet = 1
                south_Boundary[2, 1] = float(input("\nEnter xmin for Outlet in South Boundary\n"))
                southOutlet_min_index = int((south_Boundary[2, 1]-xMin)/dx)
                if(southOutlet_min_index==0):
                    southOutlet_min_index+=1
                south_Boundary[2, 2] = float(input("\nEnter xmax for Outlet in South Boundary\n"))
                southOutlet_max_index = int((south_Boundary[2, 2]-xMin)/dx)
                if(southOutlet_max_index==nx):
                    southOutlet_max_index-=1
                south_Boundary[2, 5] = float(input("\nEnter Pressure for Outlet in South Boundary\n"))
                # Boundary Marker
                D[0][-1, southOutlet_min_index:southOutlet_max_index+1] = -1.0
                # U-vel - Uninitialised
                D[1][-1, southOutlet_min_index:southOutlet_max_index+1] = -1000.5
                # V-vel - Uninitialised
                D[2][-1, southOutlet_min_index:southOutlet_max_index+1] = -1000.5
                # Pressure - User Initialised
                D[3][-1, southOutlet_min_index:southOutlet_max_index+1] = south_Boundary[2, 5]
            
            
            print("\n\n------------------SOUTH BOUNDARY------------------\n")
            boundary_info.write("\n\n------------------SOUTH BOUNDARY------------------\n")
            ite=0
            while ite<nx+1:
                if(D[0][-1, ite]==0):
                    ite1=ite+1
                    if(ite1==nx+1 and D[0][-1, ite-1]!=0):
                        print("\n[Wall]\n")
                        print("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                        print("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1-1)*dx, ite1-1))
                        print("\nNo. of elements = [{}] \n".format(ite1-1-ite+1))
                        boundary_info.write("\n[Wall]\n")
                        boundary_info.write("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                        boundary_info.write("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1-1)*dx, ite1-1))
                        boundary_info.write("\nNo. of elements = [{}] \n".format(ite1-1-ite+1))
                    while ite1<nx+1:
                        if(D[0][-1, ite1]!=0):
                            print("\n[Wall]\n")
                            print("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                            print("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1-1)*dx, ite1-1))
                            print("\nNo. of elements = [{}] \n".format(ite1-1-ite+1))
                            boundary_info.write("\n[Wall]\n")
                            boundary_info.write("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                            boundary_info.write("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1-1)*dx, ite1-1))
                            boundary_info.write("\nNo. of elements = [{}] \n".format(ite1-1-ite+1))
                            if(flag_SlidingWall_SouthBoundary==1):
                                print("\nSliding Wall Velocity {}\n".format(D[1][-1, ite1-1]))
                                boundary_info.write("\nSliding Wall Velocity {}\n".format(D[1][-1, ite1-1]))
                            ite=ite1-1
                            break
                        if(ite1==nx):
                            print("\n[Wall]\n")
                            print("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                            print("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1)*dx, ite1))
                            print("\nNo. of elements = [{}] \n".format(ite1-ite+1))
                            boundary_info.write("\n[Wall]\n")
                            boundary_info.write("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                            boundary_info.write("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1)*dx, ite1))
                            boundary_info.write("\nNo. of elements = [{}] \n".format(ite1-ite+1))
                            if(flag_SlidingWall_SouthBoundary==1):
                                print("\nSliding Wall Velocity {}\n".format(D[1][-1, ite1-1]))
                                boundary_info.write("\nSliding Wall Velocity {}\n".format(D[1][-1, ite1-1]))
                            ite=ite1-1
                            break
                        ite1=ite1+1
                if(D[0][-1, ite]==1):
                    ite1 = ite+1
                    while ite1<nx+1:
                        if(D[0][-1, ite1]!=1):
                            print("\n[Inlet]\n")
                            print("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                            print("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1-1)*dx, ite1-1))
                            print("\nNo. of elements = [{}] \n".format(ite1-1-ite+1))
                            print("\nInlet V-Velocity {}\n".format(D[2][-1, ite1-1]))
                            boundary_info.write("\n[Inlet]\n")
                            boundary_info.write("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                            boundary_info.write("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1-1)*dx, ite1-1))
                            boundary_info.write("\nNo. of elements = [{}] \n".format(ite1-1-ite+1))
                            boundary_info.write("\nInlet V-Velocity {}\n".format(D[2][-1, ite1-1]))
                            ite=ite1-1
                            break
                        if(ite1==nx):
                            print("\n[Inlet]\n")
                            print("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                            print("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1)*dx, ite1))
                            print("\nNo. of elements = [{}] \n".format(ite1-ite+1))
                            print("\nInlet V-Velocity {}\n".format(D[2][-1, ite1-1]))
                            boundary_info.write("\n[Inlet]\n")
                            boundary_info.write("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                            boundary_info.write("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1)*dx, ite1))
                            boundary_info.write("\nNo. of elements = [{}] \n".format(ite1-ite+1))
                            boundary_info.write("\nInlet V-Velocity {}\n".format(D[2][-1, ite1-1]))
                            ite=ite1-1
                            break
                        ite1=ite1+1
                if(D[0][-1, ite]==-1):
                    ite1=ite+1
                    while ite1<nx+1:
                        if(D[0][-1, ite1]!=-1):
                            print("\n[Outlet]\n")
                            print("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                            print("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1-1)*dx, ite1-1))
                            print("\nNo. of elements = [{}] \n".format(ite1-1-ite+1))
                            print("\nOutlet Pressure {}\n".format(D[3][-1, ite1-1]))
                            boundary_info.write("\n[Outlet]\n")
                            boundary_info.write("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                            boundary_info.write("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1-1)*dx, ite1-1))
                            boundary_info.write("\nNo. of elements = [{}] \n".format(ite1-1-ite+1))
                            boundary_info.write("\nOutlet Pressure {}\n".format(D[3][-1, ite1-1]))
                            ite=ite1-1
                            break
                        if(ite1==nx):
                            print("\n[Outlet]\n")
                            print("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                            print("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1)*dx, ite1))
                            print("\nNo. of elements = [{}] \n".format(ite1-ite+1))
                            print("\nOutlet Pressure {}\n".format(D[3][-1, ite1-1]))
                            boundary_info.write("\n[Outlet]\n")
                            boundary_info.write("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                            boundary_info.write("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1)*dx, ite1))
                            boundary_info.write("\nNo. of elements = [{}] \n".format(ite1-ite+1))
                            boundary_info.write("\nOutlet Pressure {}\n".format(D[3][-1, ite1-1]))
                            ite=ite1-1
                            break
                        ite1=ite1+1
                ite=ite+1
                            
            flag = input("\nFinished with South Boundary Condition ?? [y/n]\n")
            if(flag=='y'):
                flag_SouthBoundary = 1
                print("\n--------------------------------------------------\n")
                boundary_info.write("\n--------------------------------------------------\n")
        
        # Northern Boundary Defining Loop    
        
        while (flag_NorthBoundary < 1):
            
            flag_Inlet_NorthBoundary = 0
            flag_Outlet_NorthBoundary = 0
            flag_SlidingWall_NorthBoundary = 0
            flag = input("\nInlet in North Boundary ?? [y/n]\n")
            if(flag=='y'):
                flag_Inlet_NorthBoundary = 1
            flag = input("\nSliding Wall in North Boundary ?? [y/n]\n")
            if(flag=='y'):
                flag_SlidingWall_NorthBoundary = 1
            flag = input("\nOutlet in North Boundary ?? [y/n]\n")
            if(flag=='y'):
                flag_Outlet_NorthBoundary = 1         
            
            # Initialising Boundary Markers
            north_Boundary = np.zeros((3, 6))
            north_Boundary[1, 0] = 1
            north_Boundary[2, 0] = -1
            # Temporary Wall
            north_Boundary[0, 1] = xMin
            north_Boundary[0, 2] = xMax
            D[0][0, :] = 0
            D[1][0, :] = 0
            D[2][0, :] = 0
            # North Boundary Sliding U-velocity
            if(flag_SlidingWall_NorthBoundary==1):
                north_Boundary[0, 3] = float(input("\nEnter U-Velocity for Sliding Wall in North Boundary\n"))
                D[1][0, 1:nx] = north_Boundary[0, 3]
            # Inlet definition in North Boundary
            if(flag_Inlet_NorthBoundary==1):
                north_Boundary[1, 1] = float(input("\nEnter xmin for Inlet in North Boundary\n"))
                northInlet_min_index = int((north_Boundary[1, 1]-xMin)/dx)
                if(northInlet_min_index==0):
                    northInlet_min_index+=1
                north_Boundary[1, 2] = float(input("\nEnter xmax for Inlet in North Boundary\n"))
                northInlet_max_index = int((north_Boundary[1, 2]-xMin)/dx)
                if(northInlet_max_index==nx):
                    northInlet_max_index-=1
                north_Boundary[1, 4] = float(input("\nEnter V-Velocity for Inlet in North Boundary\n"))
                # Boundary Marker
                D[0][0, northInlet_min_index:northInlet_max_index+1] = 1.0
                # U-vel Initialised to 0
                D[1][0, northInlet_min_index:northInlet_max_index+1] = 0.0
                # V-vel User Initialised
                D[2][0, northInlet_min_index:northInlet_max_index+1] = north_Boundary[1, 4]
            # Outlet definition in North Boundary
            if(flag_Outlet_NorthBoundary==1):
                flag_Outlet = 1
                north_Boundary[2, 1] = float(input("\nEnter xmin for Outlet in North Boundary\n"))
                northOutlet_min_index = int((north_Boundary[2, 1]-xMin)/dx)
                if(northOutlet_min_index==0):
                    northOutlet_min_index+=1
                north_Boundary[2, 2] = float(input("\nEnter xmax for Outlet in North Boundary\n"))
                northOutlet_max_index = int((north_Boundary[2, 2]-xMin)/dx)
                if(northOutlet_max_index==nx):
                    northOutlet_max_index-=1
                north_Boundary[2, 5] = float(input("\nEnter Pressure for Outlet in North Boundary\n"))
                # Boundary Marker
                D[0][0, northOutlet_min_index:northOutlet_max_index+1] = -1.0
                # U-vel Uninitialised
                D[1][0, northOutlet_min_index:northOutlet_max_index+1] = -1000.5
                # V-vel Uninitialised
                D[2][0, northOutlet_min_index:northOutlet_max_index+1] = -1000.5
                # Pressure User Initialised
                D[3][0, northOutlet_min_index:northOutlet_max_index+1] = north_Boundary[2, 5]
            
            
            print("\n\n------------------NORTH BOUNDARY------------------\n")
            boundary_info.write("\n\n------------------NORTH BOUNDARY------------------\n")
            ite=0
            while ite<nx+1:
                if(D[0][0, ite]==0):
                    ite1=ite+1
                    if(ite1==nx+1 and D[0][0, ite-1]!=0):
                        print("\n[Wall]\n")
                        print("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                        print("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1-1)*dx, ite1-1))
                        print("\nNo. of elements = [{}] \n".format(ite1-1-ite+1))
                        boundary_info.write("\n[Wall]\n")
                        boundary_info.write("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                        boundary_info.write("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1-1)*dx, ite1-1))
                        boundary_info.write("\nNo. of elements = [{}] \n".format(ite1-1-ite+1))
                    while ite1<nx+1:
                        if(D[0][0, ite1]!=0):
                            print("\n[Wall]\n")
                            print("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                            print("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1-1)*dx, ite1-1))
                            print("\nNo. of elements = [{}] \n".format(ite1-1-ite+1))
                            boundary_info.write("\n[Wall]\n")
                            boundary_info.write("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                            boundary_info.write("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1-1)*dx, ite1-1))
                            boundary_info.write("\nNo. of elements = [{}] \n".format(ite1-1-ite+1))
                            if(flag_SlidingWall_NorthBoundary==1):
                                print("\nSliding Wall Velocity {}\n".format(D[1][0, ite1-1]))
                                boundary_info.write("\nSliding Wall Velocity {}\n".format(D[1][0, ite1-1]))
                            ite=ite1-1
                            break
                        if(ite1==nx):
                            print("\n[Wall]\n")
                            print("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                            print("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1)*dx, ite1))
                            print("\nNo. of elements = [{}] \n".format(ite1-ite+1))
                            boundary_info.write("\n[Wall]\n")
                            boundary_info.write("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                            boundary_info.write("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1)*dx, ite1))
                            boundary_info.write("\nNo. of elements = [{}] \n".format(ite1-ite+1))
                            if(flag_SlidingWall_NorthBoundary==1):
                                print("\nSliding Wall Velocity {}\n".format(D[1][0, ite1-1]))
                                boundary_info.write("\nSliding Wall Velocity {}\n".format(D[1][0, ite1-1]))
                            ite=ite1-1
                            break
                        ite1=ite1+1
                if(D[0][0, ite]==1):
                    ite1 = ite+1
                    while ite1<nx+1:
                        if(D[0][0, ite1]!=1):
                            print("\n[Inlet]\n")
                            print("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                            print("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1-1)*dx, ite1-1))
                            print("\nNo. of elements = [{}] \n".format(ite1-1-ite+1))
                            print("\nInlet V-Velocity {}\n".format(D[2][0, ite1-1]))
                            boundary_info.write("\n[Inlet]\n")
                            boundary_info.write("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                            boundary_info.write("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1-1)*dx, ite1-1))
                            boundary_info.write("\nNo. of elements = [{}] \n".format(ite1-1-ite+1))
                            boundary_info.write("\nInlet V-Velocity {}\n".format(D[2][0, ite1-1]))
                            ite=ite1-1
                            break
                        if(ite1==nx):
                            print("\n[Inlet]\n")
                            print("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                            print("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1)*dx, ite1))
                            print("\nNo. of elements = [{}] \n".format(ite1-ite+1))
                            print("\nInlet V-Velocity {}\n".format(D[2][0, ite1-1]))
                            boundary_info.write("\n[Inlet]\n")
                            boundary_info.write("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                            boundary_info.write("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1)*dx, ite1))
                            boundary_info.write("\nNo. of elements = [{}] \n".format(ite1-ite+1))
                            boundary_info.write("\nInlet V-Velocity {}\n".format(D[2][0, ite1-1]))
                            ite=ite1-1
                            break
                        ite1=ite1+1
                if(D[0][0, ite]==-1):
                    ite1=ite+1
                    while ite1<nx+1:
                        if(D[0][0, ite1]!=-1):
                            print("\n[Outlet]\n")
                            print("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                            print("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1-1)*dx, ite1-1))
                            print("\nNo. of elements = [{}] \n".format(ite1-1-ite+1))
                            print("\nOutlet Pressure {}\n".format(D[3][0, ite1-1]))
                            boundary_info.write("\n[Outlet]\n")
                            boundary_info.write("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                            boundary_info.write("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1-1)*dx, ite1-1))
                            boundary_info.write("\nNo. of elements = [{}] \n".format(ite1-1-ite+1))
                            boundary_info.write("\nOutlet Pressure {}\n".format(D[3][0, ite1-1]))
                            ite=ite1-1
                            break
                        if(ite1==nx):
                            print("\n[Outlet]\n")
                            print("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                            print("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1)*dx, ite1))
                            print("\nNo. of elements = [{}] \n".format(ite1-ite+1))
                            print("\nOutlet Pressure {}\n".format(D[3][0, ite1-1]))
                            boundary_info.write("\n[Outlet]\n")
                            boundary_info.write("\nxMin = [{}] Starting Index = [{}]\n ".format(xMin+ite*dx, ite))
                            boundary_info.write("\nxMax = [{}] Ending Index = [{}]\n ".format(xMin+(ite1)*dx, ite1))
                            boundary_info.write("\nNo. of elements = [{}] \n".format(ite1-ite+1))
                            boundary_info.write("\nOutlet Pressure {}\n".format(D[3][0, ite1-1]))
                            ite=ite1-1
                            break
                        ite1=ite1+1
                ite=ite+1
                            
            flag = input("\nFinished with North Boundary Condition ?? [y/n]\n")
            if(flag=='y'):
                flag_NorthBoundary = 1
                print("\n--------------------------------------------------\n")
                boundary_info.write("\n--------------------------------------------------\n")
                
        # Western Boundary Defining Loop
        
        while (flag_WestBoundary < 1):
            
            flag_Inlet_WestBoundary = 0
            flag_Outlet_WestBoundary = 0
            flag_SlidingWall_WestBoundary = 0
            flag = input("\nInlet in West Boundary ?? [y/n]\n")
            if(flag=='y'):
                flag_Inlet_WestBoundary = 1
            flag = input("\nSliding Wall in West Boundary ?? [y/n]\n")
            if(flag=='y'):
                flag_SlidingWall_WestBoundary = 1
            flag = input("\nOutlet in West Boundary ?? [y/n]\n")
            if(flag=='y'):
                flag_Outlet_WestBoundary = 1         
            
            # Initialising Boundary Markers
            west_Boundary = np.zeros((3, 6))
            west_Boundary[1, 0] = 1
            west_Boundary[2, 0] = -1
            # Temporary Wall
            west_Boundary[0, 1] = yMin
            west_Boundary[0, 2] = yMax
            D[0][:, 0] = 0
            D[1][:, 0] = 0
            D[2][:, 0] = 0
            # West Boundary Sliding U-velocity
            if(flag_SlidingWall_WestBoundary==1):
                west_Boundary[0, 4] = float(input("\nEnter V-Velocity for Sliding Wall in West Boundary\n"))
                D[2][1:ny, 0] = west_Boundary[0, 4]
            # Inlet definition in West Boundary
            if(flag_Inlet_WestBoundary==1):
                west_Boundary[1, 1] = float(input("\nEnter ymin for Inlet in West Boundary\n"))
                westInlet_bottom_index = int((yMax-west_Boundary[1, 1])/dy)
                if(westInlet_bottom_index==ny):
                    westInlet_bottom_index-=1
                west_Boundary[1, 2] = float(input("\nEnter ymax for Inlet in West Boundary\n"))
                westInlet_top_index = int((yMax-west_Boundary[1, 2])/dy)
                if(westInlet_top_index==0):
                    westInlet_top_index+=1
                west_Boundary[1, 3] = float(input("\nEnter U-Velocity for Inlet in West Boundary\n"))
                # Boundary Marker
                D[0][westInlet_top_index:westInlet_bottom_index+1, 0] = 1.0
                # U-vel User Initialised
                D[1][westInlet_top_index:westInlet_bottom_index+1, 0] = west_Boundary[1, 3]
                # V-vel Initialised to 0
                D[2][westInlet_top_index:westInlet_bottom_index+1, 0] = 0.0
            # Outlet definition in West Boundary
            if(flag_Outlet_WestBoundary==1):
                flag_Outlet = 1
                west_Boundary[2, 1] = float(input("\nEnter ymin for Outlet in West Boundary\n"))
                westOutlet_bottom_index = int((yMax-west_Boundary[2, 1])/dy)
                if(westOutlet_bottom_index==ny):
                    westOutlet_bottom_index-=1
                west_Boundary[2, 2] = float(input("\nEnter ymax for Outlet in West Boundary\n"))
                westOutlet_top_index = int((yMax-west_Boundary[2, 2])/dy)
                if(westOutlet_top_index==0):
                    westOutlet_top_index+=1
                west_Boundary[2, 5] = float(input("\nEnter Pressure for Outlet in West Boundary\n"))
                # Boundary Marker
                D[0][westOutlet_top_index:westOutlet_bottom_index+1, 0] = -1.0
                # U-vel Uninitialised
                D[1][westOutlet_top_index:westOutlet_bottom_index+1, 0] = -1000.5
                # V-vel Uninitialised
                D[2][westOutlet_top_index:westOutlet_bottom_index+1, 0] = -1000.5
                # Pressure User Initialised
                D[3][westOutlet_top_index:westOutlet_bottom_index+1, 0] = west_Boundary[2, 5]
            
            
            print("\n\n------------------WEST BOUNDARY-----------------\n")
            boundary_info.write("\n\n------------------WEST BOUNDARY-----------------\n")
            ite=ny
            while ite>-1:
                if(D[0][ite, 0]==0):
                    ite1=ite-1
                    if(ite1==-1 and D[0][ite+1, 0]!=0):
                        print("\n[Wall]\n")
                        print("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                        print("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1+1)*dy, ite1+1))
                        print("\nNo. of elements = [{}] \n".format(ite-ite1))
                        boundary_info.write("\n[Wall]\n")
                        boundary_info.write("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                        boundary_info.write("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1+1)*dy, ite1+1))
                        boundary_info.write("\nNo. of elements = [{}] \n".format(ite-ite1))
                    while ite1>-1:
                        if(D[0][ite1, 0]!=0):
                            print("\n[Wall]\n")
                            print("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                            print("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1+1)*dy, ite1+1))
                            print("\nNo. of elements = [{}] \n".format(ite-ite1))
                            boundary_info.write("\n[Wall]\n")
                            boundary_info.write("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                            boundary_info.write("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1+1)*dy, ite1+1))
                            boundary_info.write("\nNo. of elements = [{}] \n".format(ite-ite1))
                            if(flag_SlidingWall_WestBoundary==1):
                                print("\nSliding Wall Velocity {}\n".format(D[2][ite1+1, 0]))
                                boundary_info.write("\nSliding Wall Velocity {}\n".format(D[2][ite1+1, 0]))
                            ite=ite1+1
                            break
                        if(ite1==0):
                            print("\n[Wall]\n")
                            print("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                            print("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1)*dy, ite1))
                            print("\nNo. of elements = [{}] \n".format(ite-ite1+1))
                            boundary_info.write("\n[Wall]\n")
                            boundary_info.write("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                            boundary_info.write("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1)*dy, ite1))
                            boundary_info.write("\nNo. of elements = [{}] \n".format(ite-ite1+1))
                            if(flag_SlidingWall_WestBoundary==1):
                                print("\nSliding Wall Velocity {}\n".format(D[2][ite1+1, 0]))
                                boundary_info.write("\nSliding Wall Velocity {}\n".format(D[2][ite1+1, 0]))
                            ite=ite1+1
                            break
                        ite1=ite1-1
                if(D[0][ite, 0]==1):
                    ite1 = ite-1
                    while ite1>-1:
                        if(D[0][ite1, 0]!=1):
                            print("\n[Inlet]\n")
                            print("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                            print("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1+1)*dy, ite1+1))
                            print("\nNo. of elements = [{}] \n".format(ite-ite1))
                            print("\nInlet U-Velocity {}\n".format(D[1][ite1+1, 0]))
                            boundary_info.write("\n[Inlet]\n")
                            boundary_info.write("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                            boundary_info.write("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1+1)*dy, ite1+1))
                            boundary_info.write("\nNo. of elements = [{}] \n".format(ite-ite1))
                            boundary_info.write("\nInlet U-Velocity {}\n".format(D[1][ite1+1, 0]))
                            ite=ite1+1
                            break
                        if(ite1==0):
                            print("\n[Inlet]\n")
                            print("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                            print("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1)*dy, ite1))
                            print("\nNo. of elements = [{}] \n".format(ite-ite1+1))
                            print("\nInlet U-Velocity {}\n".format(D[1][ite1+1, 0]))
                            boundary_info.write("\n[Inlet]\n")
                            boundary_info.write("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                            boundary_info.write("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1)*dy, ite1))
                            boundary_info.write("\nNo. of elements = [{}] \n".format(ite-ite1+1))
                            boundary_info.write("\nInlet U-Velocity {}\n".format(D[1][ite1+1, 0]))
                            ite=ite1+1
                            break
                        ite1=ite1-1
                if(D[0][ite, 0]==-1):
                    ite1=ite-1
                    while ite1>-1:
                        if(D[0][ite1, 0]!=-1):
                            print("\n[Outlet]\n")
                            print("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                            print("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1+1)*dy, ite1+1))
                            print("\nNo. of elements = [{}] \n".format(ite-ite1))
                            print("\nOutlet Pressure {}\n".format(D[3][ite1+1, 0]))
                            boundary_info.write("\n[Outlet]\n")
                            boundary_info.write("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                            boundary_info.write("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1+1)*dy, ite1+1))
                            boundary_info.write("\nNo. of elements = [{}] \n".format(ite-ite1))
                            boundary_info.write("\nOutlet Pressure {}\n".format(D[3][ite1+1, 0]))
                            ite=ite1+1
                            break
                        if(ite1==0):
                            print("\n[Outlet]\n")
                            print("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                            print("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1)*dy, ite1))
                            print("\nNo. of elements = [{}] \n".format(ite-ite1+1))
                            print("\nOutlet Pressure {}\n".format(D[3][ite1+1, 0]))
                            boundary_info.write("\n[Outlet]\n")
                            boundary_info.write("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                            boundary_info.write("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1)*dy, ite1))
                            boundary_info.write("\nNo. of elements = [{}] \n".format(ite-ite1+1))
                            boundary_info.write("\nOutlet Pressure {}\n".format(D[3][ite1+1, 0]))
                            ite=ite1+1
                            break
                        ite1=ite1-1
                ite=ite-1
                            
            flag = input("\nFinished with West Boundary Condition ?? [y/n]\n")
            if(flag=='y'):
                flag_WestBoundary = 1
                print("\n-------------------------------------------------\n")
                boundary_info.write("\n-------------------------------------------------\n")
        
        # Eastern Boundary Defining Loop
        
        while (flag_EastBoundary < 1):
            
            flag_Inlet_EastBoundary = 0
            flag_Outlet_EastBoundary = 0
            flag_SlidingWall_EastBoundary = 0
            flag = input("\nInlet in East Boundary ?? [y/n]\n")
            if(flag=='y'):
                flag_Inlet_EastBoundary = 1
            flag = input("\nSliding Wall in East Boundary ?? [y/n]\n")
            if(flag=='y'):
                flag_SlidingWall_EastBoundary = 1
            flag = input("\nOutlet in East Boundary ?? [y/n]\n")
            if(flag=='y'):
                flag_Outlet_EastBoundary = 1         
            
            # Initialising Boundary Markers
            east_Boundary = np.zeros((3, 6))
            east_Boundary[1, 0] = 1
            east_Boundary[2, 0] = -1
            # Temporary Wall
            east_Boundary[0, 1] = yMin
            east_Boundary[0, 2] = yMax
            D[0][:, -1] = 0
            D[1][:, -1] = 0
            D[2][:, -1] = 0
            # East Boundary Sliding U-velocity
            if(flag_SlidingWall_EastBoundary==1):
                east_Boundary[0, 4] = float(input("Enter V-Velocity for Sliding Wall in East Boundary\n"))
                D[2][1:ny, -1] = east_Boundary[0, 4]
            # Inlet definition in East Boundary
            if(flag_Inlet_EastBoundary==1):
                east_Boundary[1, 1] = float(input("Enter ymin for Inlet in East Boundary\n"))
                eastInlet_bottom_index = int((yMax-east_Boundary[1, 1])/dy)
                if(eastInlet_bottom_index==ny):
                    eastInlet_bottom_index-=1
                east_Boundary[1, 2] = float(input("Enter ymax for Inlet in East Boundary\n"))
                eastInlet_top_index = int((yMax-east_Boundary[1, 2])/dy)
                if(eastInlet_top_index==0):
                    eastInlet_top_index+=1
                east_Boundary[1, 3] = float(input("Enter U-Velocity for Inlet in East Boundary\n"))
                # Boundary Marker
                D[0][eastInlet_top_index:eastInlet_bottom_index+1, -1] = 1.0
                # U-vel User Initialised
                D[1][eastInlet_top_index:eastInlet_bottom_index+1, -1] = east_Boundary[1, 3]
                # V-vel Initialised to 0
                D[2][eastInlet_top_index:eastInlet_bottom_index+1, -1] = 0.0
            # Outlet definition in East Boundary
            if(flag_Outlet_EastBoundary==1):
                flag_Outlet = 1
                east_Boundary[2, 1] = float(input("Enter ymin for Outlet in East Boundary\n"))
                eastOutlet_bottom_index = int((yMax-east_Boundary[2, 1])/dy)
                if(eastOutlet_bottom_index==ny):
                    eastOutlet_bottom_index-=1
                east_Boundary[2, 2] = float(input("Enter ymax for Outlet in East Boundary\n"))
                eastOutlet_top_index = int((yMax-east_Boundary[2, 2])/dy)
                if(eastOutlet_top_index==0):
                    eastOutlet_top_index+=1
                east_Boundary[2, 5] = float(input("Enter Pressure for Outlet in East Boundary\n"))
                # Boundary Marker
                D[0][eastOutlet_top_index:eastOutlet_bottom_index+1, -1] = -1.0
                # U-vel Uninitialised
                D[1][eastOutlet_top_index:eastOutlet_bottom_index+1, -1] = -1000.5
                # V-vel Uninitialised
                D[2][eastOutlet_top_index:eastOutlet_bottom_index+1, -1] = -1000.5
                # Pressure User Initialised
                D[3][eastOutlet_top_index:eastOutlet_bottom_index+1, -1] = east_Boundary[2, 5]
            
            
            print("\n\n------------------EAST BOUNDARY-----------------\n")
            boundary_info.write("\n\n------------------EAST BOUNDARY-----------------\n")
            ite=ny
            while ite>-1:
                if(D[0][ite, -1]==0):
                    ite1=ite-1
                    if(ite1==-1 and D[0][ite+1, -1]!=0):
                        print("\n[Wall]\n")
                        print("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                        print("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1+1)*dy, ite1+1))
                        print("\nNo. of elements = [{}] \n".format(ite-ite1))
                        boundary_info.write("\n[Wall]\n")
                        boundary_info.write("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                        boundary_info.write("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1+1)*dy, ite1+1))
                        boundary_info.write("\nNo. of elements = [{}] \n".format(ite-ite1))
                    while ite1>-1:
                        if(D[0][ite1, -1]!=0):
                            print("\n[Wall]\n")
                            print("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                            print("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1+1)*dy, ite1+1))
                            print("\nNo. of elements = [{}]\n".format(ite-ite1))
                            boundary_info.write("\n[Wall]\n")
                            boundary_info.write("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                            boundary_info.write("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1+1)*dy, ite1+1))
                            boundary_info.write("\nNo. of elements = [{}]\n".format(ite-ite1))
                            if(flag_SlidingWall_EastBoundary==1):
                                print("\nSliding Wall Velocity {}\n".format(D[2][ite1+1, -1]))
                                boundary_info.write("\nSliding Wall Velocity {}\n".format(D[2][ite1+1, -1]))
                            ite=ite1+1
                            break
                        if(ite1==0):
                            print("\n[Wall]\n")
                            print("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                            print("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1)*dy, ite1))
                            print("\nNo. of elements = [{}]\n".format(ite-ite1+1))
                            boundary_info.write("\n[Wall]\n")
                            boundary_info.write("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                            boundary_info.write("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1)*dy, ite1))
                            boundary_info.write("\nNo. of elements = [{}]\n".format(ite-ite1+1))
                            if(flag_SlidingWall_EastBoundary==1):
                                print("\nSliding Wall Velocity {}\n".format(D[2][ite1+1, -1]))
                                boundary_info.write("\nSliding Wall Velocity {}\n".format(D[2][ite1+1, -1]))
                            ite=ite1+1
                            break
                        ite1=ite1-1
                if(D[0][ite, -1]==1):
                    ite1 = ite-1
                    while ite1>-1:
                        if(D[0][ite1, -1]!=1):
                            print("\n[Inlet]\n")
                            print("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                            print("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1+1)*dy, ite1+1))
                            print("\nNo. of elements = [{}] \n".format(ite-ite1))
                            print("\nInlet U-Velocity {}\n".format(D[1][ite1+1, -1]))
                            boundary_info.write("\n[Inlet]\n")
                            boundary_info.write("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                            boundary_info.write("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1+1)*dy, ite1+1))
                            boundary_info.write("\nNo. of elements = [{}] \n".format(ite-ite1))
                            boundary_info.write("\nInlet U-Velocity {}\n".format(D[1][ite1+1, -1]))
                            ite=ite1+1
                            break
                        if(ite1==0):
                            print("\n[Inlet]\n")
                            print("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                            print("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1)*dy, ite1))
                            print("\nNo. of elements = [{}] \n".format(ite-ite1+1))
                            print("\nInlet U-Velocity {}\n".format(D[1][ite1+1, -1]))
                            boundary_info.write("\n[Inlet]\n")
                            boundary_info.write("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                            boundary_info.write("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1)*dy, ite1))
                            boundary_info.write("\nNo. of elements = [{}] \n".format(ite-ite1+1))
                            boundary_info.write("\nInlet U-Velocity {}\n".format(D[1][ite1+1, -1]))
                            ite=ite1+1
                            break
                        ite1=ite1-1
                if(D[0][ite, -1]==-1):
                    ite1=ite-1
                    while ite1>-1:
                        if(D[0][ite1, -1]!=-1):
                            print("\n[Outlet]\n")
                            print("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                            print("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1+1)*dy, ite1+1))
                            print("\nNo. of elements = [{}] \n".format(ite-ite1))
                            print("\nOutlet Pressure {}\n".format(D[3][ite1+1, -1]))
                            boundary_info.write("\n[Outlet]\n")
                            boundary_info.write("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                            boundary_info.write("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1+1)*dy, ite1+1))
                            boundary_info.write("\nNo. of elements = [{}] \n".format(ite-ite1))
                            boundary_info.write("\nOutlet Pressure {}\n".format(D[3][ite1+1, -1]))
                            ite=ite1+1
                            break
                        if(ite1==0):
                            print("\n[Outlet]\n")
                            print("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                            print("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1)*dy, ite1))
                            print("\nNo. of elements = [{}] \n".format(ite-ite1+1))
                            print("\nOutlet Pressure {}\n".format(D[3][ite1+1, -1]))
                            boundary_info.write("\n[Outlet]\n")
                            boundary_info.write("\nyMin = [{:.2f}] Bottom Index = [{}]\n ".format(yMax-(ite)*dy, ite))
                            boundary_info.write("\nyMax = [{:.2f}] Top Index = [{}]\n ".format(yMax-(ite1)*dy, ite1))
                            boundary_info.write("\nNo. of elements = [{}] \n".format(ite-ite1+1))
                            boundary_info.write("\nOutlet Pressure {}\n".format(D[3][ite1+1, -1]))
                            ite=ite1+1
                            break
                        ite1=ite1-1
                ite=ite-1
                            
            flag = input("Finished with East Boundary Condition ?? [y/n]")
            if(flag=='y'):
                flag_EastBoundary = 1
                print("\n-------------------------------------------------\n")
                boundary_info.write("\n-------------------------------------------------\n")
        flag_Boundary = 1
        
    return D, xMin, xMax, yMin, yMax, dx, dy, nx, ny, south_Boundary, north_Boundary, west_Boundary, east_Boundary, flag_Outlet

def obstacle(D, xMin, xMax, yMin, yMax, nx, ny, dx, dy):
    package_dir = os.path.dirname(os.path.realpath(__file__))
    flag_obstacle = 'n'
    while flag_obstacle!='y':
        scal_x = float(input("\nEnter Scaling factor for X-direction\n"))
        scal_y = float(input("\nEnter Scaling factor for Y-direction\n"))
        traf_x = float(input("\nOffset in X\n"))
        traf_y = float(input("\nOffset in Y\n")) + dy/10
        rot = float(input("\nEnter anti-clockwise rotation angle\n"))
        flag_geometryCase = input("\n1. Enter address of .txt file\n2. Use default address\n")
        if(flag_geometryCase==1):
            geometryCase = input("\nEnter address of .txt file\nExample: /home/user/folder/Case_1/Points.txt\n")
        else:
            geometryCase = package_dir+"/Points.txt"
        file1 = open(geometryCase, 'r')
        Lines = file1.readlines()
        P=[]
        N=0
        for line in Lines:
            N+=1
            P.append(np.fromstring(line.strip(), dtype=float, sep=' '))
        N+=1
        P.append(np.fromstring(Lines[0].strip(), dtype=float, sep=' '))
        Lx = xMax-xMin
        Ly = yMax-yMin
        Points_x=[]
        Points_y=[]
        for i in range (0, N):
            Points_x.append(np.sqrt((P[i][0]*scal_x)**2 + (P[i][1]*scal_y)**2)*np.cos(rot*np.pi/180+np.arctan((P[i][1])/(P[i][0])))+traf_x)
            Points_y.append(np.sqrt((P[i][0]*scal_x)**2 + (P[i][1]*scal_y)**2)*np.sin(rot*np.pi/180+np.arctan(P[i][1]/P[i][0]))+traf_y)
        fig = plt.subplots(figsize =(Lx*5, Ly*5))
        plt.xlim(xMin, xMax)
        plt.ylim(yMin, yMax)
        plt.fill(Points_x, Points_y, color = 'black')
        plt.show()
        
        flag_obstacle = (input("\nConitnue [Y/N] ??\n"))
        if(flag_obstacle!='y'):
            continue
        for i in range(0, ny+1):
            for j in range (0, nx+1):
                X = xMin + j*dx
                Y = yMax - i*dy
                count=0
                for ite in range (0, N-1):
                    X1 = Points_x[ite]
                    X2 = Points_x[ite+1]
                    Y1 = Points_y[ite]
                    Y2 = Points_y[ite+1]
                    if(Y1<=Y<=Y2 or Y2<=Y<=Y1):
                        if(X1<=X2):
                            if(X<=(X2+(Y-Y2)*(X1-X2)/(Y1-Y2))):
                                count +=1
                        if(X2<=X1):
                            if(X<=(X1+(Y-Y1)*(X2-X1)/(Y2-Y1))):
                                count +=1
                if(count==1):
                    D[0][i, j] = 0
                    D[1][i, j] = 0
                    D[2][i, j] = 0
                    D[3][i, j] = -100000
    
    return D, P, Points_x, Points_y

def obstacleMath(D, xMin, xMax, yMin, yMax, nx, ny, dx, dy):
    flag_obstacle = 'n'
    while flag_obstacle!='y':
        scal_x = float(input("\nEnter Scaling factor for X-direction\n"))
        scal_y = float(input("\nEnter Scaling factor for Y-direction\n"))
        traf_x = float(input("\nOffset in X\n"))
        traf_y = float(input("\nOffset in Y\n")) + dy/10
        rot = int(input("\nEnter anti-clockwise rotation angle\n"))
            
        P=[]
        asq = float(input("\nEnter Horizontal axis length for Elliptical Geometry\n"))
        bsq = float(input("\nEnter Vertical axis length for Elliptical Geometry\n"))
        a=asq**2
        b=bsq**2
        N=0
        for theta in range (0, 360):
            if(theta%2==0):
                temp1=a*b/(b+a*np.tan(theta*np.pi/180)**2)
                if(0<=theta<=90):
                    temp = [np.abs(np.sqrt(temp1)), np.abs(np.tan(theta*np.pi/180))*(np.abs(np.sqrt(temp1)))]
                elif(90<=theta<=180):
                    temp = [-np.abs(np.sqrt(temp1)), np.abs(np.tan(theta*np.pi/180))*(np.abs(np.sqrt(temp1)))]
                elif(180<=theta<=270):
                    temp = [-np.abs(np.sqrt(temp1)), -np.abs(np.tan(theta*np.pi/180))*(np.abs(np.sqrt(temp1)))]
                elif(270<=theta<=360):
                    temp = [np.abs(np.sqrt(temp1)), -np.abs(np.tan(theta*np.pi/180))*(np.abs(np.sqrt(temp1)))]
                P.append(temp)
                N+=1
        P.append(P[0])
        Lx = xMax-xMin
        Ly = yMax-yMin
        Points_x=[]
        Points_y=[]
        for i in range(0, 181):
            if(0<=(i+rot/2)<=45):
                Points_x.append(np.abs(np.sqrt((P[i][0]*scal_x)**2 + (P[i][1]*scal_y)**2)*np.cos(rot*np.pi/180+np.arctan(P[i][1]/P[i][0])))+traf_x)
                Points_y.append(np.abs(np.sqrt((P[i][0]*scal_x)**2 + (P[i][1]*scal_y)**2)*np.sin(rot*np.pi/180+np.arctan(P[i][1]/P[i][0])))+traf_y)
            elif(45<=(i+rot/2)<=90):
                Points_x.append(-np.abs(np.sqrt((P[i][0]*scal_x)**2 + (P[i][1]*scal_y)**2)*np.cos(rot*np.pi/180+np.arctan(P[i][1]/P[i][0])))+traf_x)
                Points_y.append(np.abs(np.sqrt((P[i][0]*scal_x)**2 + (P[i][1]*scal_y)**2)*np.sin(rot*np.pi/180+np.arctan(P[i][1]/P[i][0])))+traf_y)
            elif(90<=(i+rot/2)<=135):
                Points_x.append(-np.abs(np.sqrt((P[i][0]*scal_x)**2 + (P[i][1]*scal_y)**2)*np.cos(rot*np.pi/180+np.arctan(P[i][1]/P[i][0])))+traf_x)
                Points_y.append(-np.abs(np.sqrt((P[i][0]*scal_x)**2 + (P[i][1]*scal_y)**2)*np.sin(rot*np.pi/180+np.arctan(P[i][1]/P[i][0])))+traf_y)
            elif(135<=(i+rot/2)<=180):
                Points_x.append(np.abs(np.sqrt((P[i][0]*scal_x)**2 + (P[i][1]*scal_y)**2)*np.cos(rot*np.pi/180+np.arctan(P[i][1]/P[i][0])))+traf_x)
                Points_y.append(-np.abs(np.sqrt((P[i][0]*scal_x)**2 + (P[i][1]*scal_y)**2)*np.sin(rot*np.pi/180+np.arctan(P[i][1]/P[i][0])))+traf_y)
            elif(180<=(i+rot/2)<=225):
                Points_x.append(np.abs(np.sqrt((P[i][0]*scal_x)**2 + (P[i][1]*scal_y)**2)*np.cos(rot*np.pi/180+np.arctan(P[i][1]/P[i][0])))+traf_x)
                Points_y.append(np.abs(np.sqrt((P[i][0]*scal_x)**2 + (P[i][1]*scal_y)**2)*np.sin(rot*np.pi/180+np.arctan(P[i][1]/P[i][0])))+traf_y)
            elif(-45<=(i+rot/2)<=0):
                Points_x.append(np.abs(np.sqrt((P[i][0]*scal_x)**2 + (P[i][1]*scal_y)**2)*np.cos(rot*np.pi/180+np.arctan(P[i][1]/P[i][0])))+traf_x)
                Points_y.append(-np.abs(np.sqrt((P[i][0]*scal_x)**2 + (P[i][1]*scal_y)**2)*np.sin(rot*np.pi/180+np.arctan(P[i][1]/P[i][0])))+traf_y)
            elif(-90<=(i+rot/2)<=-45):
                Points_x.append(-np.abs(np.sqrt((P[i][0]*scal_x)**2 + (P[i][1]*scal_y)**2)*np.cos(rot*np.pi/180+np.arctan(P[i][1]/P[i][0])))+traf_x)
                Points_y.append(-np.abs(np.sqrt((P[i][0]*scal_x)**2 + (P[i][1]*scal_y)**2)*np.sin(rot*np.pi/180+np.arctan(P[i][1]/P[i][0])))+traf_y)
            elif(-135<=(i+rot/2)<=-90):
                Points_x.append(-np.abs(np.sqrt((P[i][0]*scal_x)**2 + (P[i][1]*scal_y)**2)*np.cos(rot*np.pi/180+np.arctan(P[i][1]/P[i][0])))+traf_x)
                Points_y.append(np.abs(np.sqrt((P[i][0]*scal_x)**2 + (P[i][1]*scal_y)**2)*np.sin(rot*np.pi/180+np.arctan(P[i][1]/P[i][0])))+traf_y)
            elif(-180<=(i+rot/2)<=-135):
                Points_x.append(np.abs(np.sqrt((P[i][0]*scal_x)**2 + (P[i][1]*scal_y)**2)*np.cos(rot*np.pi/180+np.arctan(P[i][1]/P[i][0])))+traf_x)
                Points_y.append(np.abs(np.sqrt((P[i][0]*scal_x)**2 + (P[i][1]*scal_y)**2)*np.sin(rot*np.pi/180+np.arctan(P[i][1]/P[i][0])))+traf_y)
        fig = plt.subplots(figsize =(Lx*5, Ly*5))
        plt.xlim(xMin, xMax)
        plt.ylim(yMin, yMax)
        plt.fill(Points_x, Points_y, color='black')
        plt.show()
        
        flag_obstacle = (input("\nConitnue [Y/N] ??\n"))
        if(flag_obstacle!='y'):
            continue
        for i in range(0, ny+1):
            for j in range (0, nx+1):
                X = xMin + j*dx
                Y = yMax - i*dy
                count=0
                for ite in range (0, 180):
                    X1 = Points_x[ite]
                    X2 = Points_x[ite+1]
                    Y1 = Points_y[ite]
                    Y2 = Points_y[ite+1]
                    if(Y1<Y<=Y2 or Y2<=Y<Y1):
                        if(X<=X1 or X<=X2):
                            count +=1
                if(count==1):
                    D[0][i, j] = 0
                    D[1][i, j] = 0
                    D[2][i, j] = 0
                    D[3][i, j] = -100000
        

    return D, P, Points_x, Points_y