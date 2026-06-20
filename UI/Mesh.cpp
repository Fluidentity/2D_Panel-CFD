#include <QTextBrowser>
#include <QMainWindow>
#include <string>
#include <QString>
#include <QTextStream>
#include "Mesh.h"
#include "Boundary.h"
#include "Domain.h"

MeshStaggered::MeshStaggered() {
  meshdomainindex=0;
  meshboundaryindex=0;

}

MeshStaggered::MeshStaggered(DomainStaggered d) {

  U_Matrix = MatrixXd::Constant(d.ny, d.nx+1, 0);
  V_Matrix = MatrixXd::Constant(d.ny+1, d.nx, 0);
  P_Matrix = MatrixXd::Constant(d.ny, d.nx, 0);
  D_Marker = d.D_boundaryMarker;
}

void MeshStaggered::updateMeshBoundary(DomainStaggered d, BoundaryStaggered b)
{
  for (int i=0; i<b.southWallctr; i++) {
    int xMin = b.southWall(i, 0);
    int xMax = b.southWall(i, 1);
    int V_Vel = b.southWall(i, 3);

    int nx_min = (xMin - d.xMin)/d.dx;
    int nx_max = (xMax - d.xMin)/d.dx;
    for (int j=nx_min; j<nx_max+1; j++)
    {
      D_Marker(d.ny+1, j) = 1;
      V_Matrix(d.ny+1, j) = V_Vel;

    }
  }

  for (int i=0; i<b.northWallctr; i++) {
    int xMin = b.northWall(i, 0);
    int xMax = b.northWall(i, 1);
    int V_Vel = b.northWall(i, 3);

    int nx_min = (xMin - d.xMin)/d.dx;
    int nx_max = (xMax - d.xMin)/d.dx;
    for (int j=nx_min; j<nx_max+1; j++)
    {
      D_Marker(0, j) = 1;
      V_Matrix(0, j) = V_Vel;

    }
  }

  for (int i=0; i<b.westWallctr; i++) {
    int yMin = b.westWall(i, 0);
    int yMax = b.westWall(i, 1);
    int U_Vel = b.westWall(i, 2);

    int ny_min = (d.yMax - yMax)/d.dy;
    int ny_max = (d.yMax - yMin)/d.dy;
    for (int j=ny_max; j<ny_min+1; j++)
    {
      D_Marker(j, 0) = 1;
      U_Matrix(j, 0) = U_Vel;

    }
  }

  for (int i=0; i<b.eastWallctr; i++) {
    int yMin = b.eastWall(i, 0);
    int yMax = b.eastWall(i, 1);
    int U_Vel = b.eastWall(i, 2);

    int ny_min = (d.yMax - yMax)/d.dy;
    int ny_max = (d.yMax - yMin)/d.dy;
    for (int j=ny_max; j<ny_min+1; j++)
    {
      D_Marker(j, d.nx+1) = 1;
      U_Matrix(j, d.nx+1) = U_Vel;

    }
  }

}

