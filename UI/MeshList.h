#ifndef MESHLIST_H
#define MESHLIST_H

#include <QMainWindow>
#include <QStandardItem>
#include <QStandardItemModel>
#include "Mesh.h"

class meshList
{
public:
  meshList();
  QStandardItem *treeNode;
  MeshStaggered *matrix;

};

#endif // MESHLIST_H
