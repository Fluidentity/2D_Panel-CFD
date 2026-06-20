#ifndef LINKMODEL_H
#define LINKMODEL_H
#include <QStandardItem>
#include <QStandardItemModel>
#include "Domain.h"

class LinkModel
{
public:
  LinkModel();
  QStandardItem tree_0;
  DomainStaggered *domain_0;
};

#endif // LINKMODEL_H
