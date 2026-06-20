#ifndef CUSTOMSTANDARDITEMMODEL_H
#define CUSTOMSTANDARDITEMMODEL_H

#include <QStandardItemModel>
#include <QComboBox>

class CustomStandardItemModel : public QStandardItemModel
{
public:
  explicit CustomStandardItemModel(QObject *parent = nullptr);
};

#endif // CUSTOMSTANDARDITEMMODEL_H
