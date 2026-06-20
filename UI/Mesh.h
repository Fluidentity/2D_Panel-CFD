#ifndef MESH_H
#define MESH_H
#include<iostream>
#include <Eigen/Dense>
#include "Domain.h"
#include "Boundary.h"



using namespace std;
using namespace Eigen;

class MeshStaggered
{
public:

int meshdomainindex;
int meshboundaryindex;

Eigen::MatrixXd U_Matrix;
Eigen::MatrixXd V_Matrix;
Eigen::MatrixXd P_Matrix;
Eigen::MatrixXd D_Marker;


MeshStaggered();


MeshStaggered(DomainStaggered d);
void updateMeshBoundary(DomainStaggered d, BoundaryStaggered b);

void displayMesh();
private:


};

#endif // MESH_H
