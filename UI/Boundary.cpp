#include "Boundary.h"


BoundaryStaggered::BoundaryStaggered() {


    southWall = MatrixXd::Constant(1, 4, 0);
    westWall = MatrixXd::Constant(1, 4, 0);
    northWall = MatrixXd::Constant(1, 4, 0);
    eastWall = MatrixXd::Constant(1, 4, 0);
    southWallctr =1;
    westWallctr = 1;
    northWallctr =1;
    eastWallctr = 1;



}




