#ifndef Domain_H
#define Domain_H
#include<iostream>
#include <Eigen/Dense>



using namespace std;
using namespace Eigen;

class DomainStaggered
{
public:

Eigen::MatrixXd D_boundaryMarker;
Eigen::MatrixXd D_UMarker;
Eigen::MatrixXd D_VMarker;
Eigen::MatrixXd D_PMarker;
float xMax;
float xMin;
float yMax;
float yMin;
float dx_temp;
float dy_temp;

float dx;
float dy;

int nx;
int ny;



DomainStaggered();
void calcDomain();





};

#endif // DIALOG_H
