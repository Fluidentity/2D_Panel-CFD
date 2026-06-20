#ifndef BOUNDARYLIST_H
#define BOUNDARYLIST_H

#include <QMainWindow>
#include <QStandardItem>
#include <QStandardItemModel>
#include "Boundary.h"

class boundaryList
{
public:
  boundaryList();
  QStandardItem *treeNode;
  BoundaryStaggered *matrix;

};

#endif // BOUNDARYLIST_H
