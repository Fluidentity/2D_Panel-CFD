#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QMouseEvent>
#include <QEvent>
#include <QDebug>
#include <QModelIndex>
#include <QTreeView>
#include <QTableView>
#include <Domain.h>
#include <Eigen/Dense>
#include <QTableWidget>
#include <QStandardItemModel>
#include <QStandardItem>
#include "Domain.h"
#include "domainlist.h"
#include "Boundary.h"
#include "boundaryList.h"
#include "MeshList.h"
#include <string>
#include  <QPushButton>
#include <list>


QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
  Q_OBJECT


public:
  MainWindow(QWidget *parent = 0);  
  ~MainWindow();
  DomainStaggered domainTable;
  BoundaryStaggered boundaryTable;
  QStandardItemModel domainTableGUI;
  QStandardItemModel defaultTree;
  MeshStaggered meshInputTable;
  int domainNo;
  int boundaryNo;
  int meshNo;
  int domainGUIIndex;
  int boundaryGUIIndex;
  int edgeGUIIndex;
  int meshGUIIndex;
  list <domainList> Dlist;
  list <boundaryList> BList;
  list <meshList> MList;



public slots:
  void saveCurrentIndex1(MeshStaggered &m, int index);
  void saveCurrentIndex2(MeshStaggered &m, int index);

private slots:

  void on_BoundaryMenu(const QPoint &pos);
  void on_domainsave_released();
  void on_pushButton_clicked();
  void on_treeView_clicked(const QModelIndex &index);
  void on_boundarysave_clicked();

  void on_pushButton_4_clicked();


  void on_btnMeshgeneration_clicked();

private:
  Ui::MainWindow *ui;

  QStandardItemModel *model;



  //void initTableView();
  void prepareDefaultTree();
  void updateDomainTableGUI(DomainStaggered d);
  void setDomainTable(DomainStaggered &d);
  void addInstanceinBoundaryWall(QModelIndex&index);



  void createDomainNode();
  void addBoundaryInstance();
  void createMesh();
  void updatesouthBoundaryTableGUI(BoundaryStaggered b);
  void updatewestBoundaryTableGUI(BoundaryStaggered b);
  void updatenorthBoundaryTableGUI(BoundaryStaggered b);
  void updateeastBoundaryTableGUI(BoundaryStaggered b);

  void updateMeshLinkTableGUI(MeshStaggered &m);

  void createBoundaryNode();
  void setsouthBoundaryTable(BoundaryStaggered &b);
  void setwestBoundaryTable(BoundaryStaggered &b);
  void setnorthBoundaryTable(BoundaryStaggered &b);
  void seteastBoundaryTable(BoundaryStaggered &b);
};
#endif // MAINWINDOW_H
