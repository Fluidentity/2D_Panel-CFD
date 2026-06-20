#ifndef DOMAINLIST_H
#define DOMAINLIST_H

#include <QMainWindow>
#include <QStandardItem>
#include <QStandardItemModel>
#include "Domain.h"

class domainList
{
public:
  domainList();
  QStandardItem *treeNode;
  DomainStaggered *matrix;

};

#endif // DOMAINLIST_H
