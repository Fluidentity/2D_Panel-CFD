#include "Domain.h"


DomainStaggered::DomainStaggered() {

    dx_temp=0.2;
    dy_temp=0.2;
    int nx, ny;

    xMax=1.0;
    xMin=0.0;
    yMax=1.0;
    yMin=0.0;
    nx = (int)((xMax-xMin)/dx_temp);
    ny = (int)((yMax-yMin)/dy_temp);

    D_boundaryMarker = MatrixXd::Constant(ny+1, nx+1, 0);
    for (int i = 0; i<ny+1; i++) {
        for(int j=0; j<nx+1; j++){
          if(i==0 || j==0 || i==ny || j==nx) {
            D_boundaryMarker(i, j) = 0;
          }

        }
    }
    D_UMarker = MatrixXd::Constant(ny+1, nx+1, 0);
    D_VMarker = MatrixXd::Constant(ny+1, nx+1, 0);
    D_PMarker = MatrixXd::Constant(ny+1, nx+1, 0);

    dx = (xMax-xMin)/nx;
    dy = (yMax-yMin)/ny;





}

void DomainStaggered::calcDomain() {
  nx = (int)((xMax-xMin)/dx_temp);
  ny = (int)((yMax-yMin)/dy_temp);

  dx = (xMax-xMin)/nx;
  dy = (yMax-yMin)/ny;
  D_boundaryMarker = MatrixXd::Constant(ny+1, nx+1, 2);
  for (int i = 0; i<ny+1; i++) {
      for(int j=0; j<nx+1; j++){
        if(i==0 || j==0 || i==ny || j==nx) {
          D_boundaryMarker(i, j) = 0;
        }

      }
  }
  D_UMarker = MatrixXd::Constant(ny+1, nx+1, 0);
  D_VMarker = MatrixXd::Constant(ny+1, nx+1, 0);
  D_PMarker = MatrixXd::Constant(ny+1, nx+1, 0);
}


