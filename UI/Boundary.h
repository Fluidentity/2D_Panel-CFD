#ifndef Boundary_H
#define Boundary_H
#include<iostream>
#include <Eigen/Dense>
#include <QMainWindow>


using namespace std;
using namespace Eigen;

class BoundaryStaggered
{

// Initialising Boundary Marker Matrix for South Wall
// Wall Edge   [ 0 xMax xMin U_vel V_vel Pressure ]
// Inlet Edge  | 1 xMax xMin U_vel V_vel Pressure |
// Outlet Edge [-1 xMax xMin U_vel V_vel Pressure ]
public:

Eigen::MatrixXd southWall;
Eigen::MatrixXd westWall;
Eigen::MatrixXd northWall;
Eigen::MatrixXd eastWall;

int southWallctr;
int westWallctr;
int northWallctr;
int eastWallctr;

BoundaryStaggered();






};

#endif // DIALOG_H
