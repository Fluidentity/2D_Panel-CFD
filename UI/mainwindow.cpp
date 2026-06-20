#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <iostream>
#include <fstream>
#include "Domain.h"
#include "Boundary.h"
#include "linkmodel.h"
#include "domainlist.h"
#include "boundaryList.h"
#include "MeshList.h"
#include <Eigen/Dense>
#include <QTreeView>
#include <QTableView>
#include <QModelIndex>
#include <QStandardItemModel>
#include <QStandardItem>
#include <string>
#include <QPushButton>
#include <QDebug>
#include <QAbstractItemModel>
#include <QTextStream>
#include <QContextMenuEvent>
#include <QMenu>
#include <QComboBox>
#include <QAbstractItemDelegate>
#include <QStyledItemDelegate>
#include <QItemDelegate>
#include "comboboxdelegate.h"
#include "dialog.h"
#include <QMessageBox>

MainWindow::MainWindow(QWidget *parent)
  : QMainWindow(parent)
  , ui(new Ui::MainWindow)
{

  domainNo = 0;
  boundaryNo=0;
  meshNo = 0;
  domainGUIIndex = 1;
  edgeGUIIndex = 1;
  boundaryGUIIndex = 1;
  meshGUIIndex=1;
  ui->setupUi(this);
  prepareDefaultTree();
  updateDomainTableGUI(domainTable);

  ui->treeView->setContextMenuPolicy(Qt::CustomContextMenu);
  connect(ui->treeView, SIGNAL(customContextMenuRequested(const QPoint &)), this, SLOT(on_BoundaryMenu(const QPoint &)));

}

MainWindow::~MainWindow()
{
  delete ui;
}

void MainWindow::updateDomainTableGUI(DomainStaggered d)
{
  ui->tableView->setModel(nullptr);
  domainTableGUI.setRowCount(6);
  domainTableGUI.setColumnCount(2);
  QStringList verticalheaders;
  QStringList horizontalheaders;
  verticalheaders << "xMin" << "xMax" << "yMin" << "ymax" << "dx" << "dy";
  horizontalheaders << "Values" << "Unit";
  domainTableGUI.setVerticalHeaderLabels(verticalheaders);
  domainTableGUI.setHorizontalHeaderLabels(horizontalheaders);


  QStandardItem *value00 = new QStandardItem;
  value00->setText(QString::number(d.xMin));
  domainTableGUI.setItem(0, 0, value00);

  QStandardItem *value01 = new QStandardItem;
  value01->setText("m");
  domainTableGUI.setItem(0, 1, value01);

  QStandardItem *value10 = new QStandardItem;
  value10->setText(QString::number(d.xMax));
  domainTableGUI.setItem(1, 0, value10);

  QStandardItem *value11 = new QStandardItem;
  value11->setText("m");
  domainTableGUI.setItem(1, 1, value11);

  QStandardItem *value20 = new QStandardItem;
  value20->setText(QString::number(d.yMin));
  domainTableGUI.setItem(2, 0, value20);

  QStandardItem *value21 = new QStandardItem;
  value21->setText("m");
  domainTableGUI.setItem(2, 1, value21);

  QStandardItem *value30 = new QStandardItem;
  value30->setText(QString::number(d.yMax));
  domainTableGUI.setItem(3, 0, value30);

  QStandardItem *value31 = new QStandardItem;
  value31->setText("m");
  domainTableGUI.setItem(3, 1, value31);

  QStandardItem *value40 = new QStandardItem;
  value40->setText(QString::number(d.dx_temp));
  domainTableGUI.setItem(4, 0, value40);

  QStandardItem *value41 = new QStandardItem;
  value41->setText("m");
  domainTableGUI.setItem(4, 1, value41);

  QStandardItem *value50 = new QStandardItem;
  value50->setText(QString::number(d.dy_temp));
  domainTableGUI.setItem(5, 0, value50);

  QStandardItem *value51 = new QStandardItem;
  value51->setText("m");
  domainTableGUI.setItem(5, 1, value51);


  ui->tableView->setModel(&domainTableGUI);
  ui->textBrowser->append("Updated GUI Table No." + QString::number(domainGUIIndex));



}

void MainWindow::updatesouthBoundaryTableGUI(BoundaryStaggered b) {
  ui->tableView->setModel(nullptr);
  domainTableGUI.setRowCount(b.southWallctr);
  domainTableGUI.setColumnCount(4);
  QStringList horizontalheaders;
  horizontalheaders << "Marker" << "X_Min" << "X_Max" << "Velocity\nPressure";
  QStringList verticalheaders;
  for (int i=0; i<b.southWallctr; i++) {
      verticalheaders << QString::number(i);
  }
  domainTableGUI.setVerticalHeaderLabels(verticalheaders);
  domainTableGUI.setHorizontalHeaderLabels(horizontalheaders);
  QStandardItem *value00;
  for (int i=0; i<b.southWallctr; i++) {
      for (int j=0; j<4; j++) {
          value00 = new QStandardItem;
          value00->setText(QString::number(b.southWall(i, j)));
          domainTableGUI.setItem(i, j, value00);
      }
    }




  ui->tableView->setModel(&domainTableGUI);
  ui->textBrowser->append("Updated GUI Table No." + QString::number(domainGUIIndex));

}

void MainWindow::updatewestBoundaryTableGUI(BoundaryStaggered b) {
  ui->tableView->setModel(nullptr);
  domainTableGUI.setRowCount(b.westWallctr);
  domainTableGUI.setColumnCount(4);
  QStringList horizontalheaders;
  horizontalheaders << "Boundary_Marker" << "X_Min" << "X_Max" << "Velocity\nPressure";
  domainTableGUI.setHorizontalHeaderLabels(horizontalheaders);
  QStringList verticalheaders;
  for (int i=0; i<b.westWallctr; i++) {
      verticalheaders << QString::number(i);
  }
  domainTableGUI.setVerticalHeaderLabels(verticalheaders);
  QStandardItem *value00;
  for (int i=0; i<b.westWallctr; i++) {
      for (int j=0; j<4; j++) {
          value00 = new QStandardItem;
          value00->setText(QString::number(b.westWall(i, j)));
          domainTableGUI.setItem(i, j, value00);
      }
    }




  ui->tableView->setModel(&domainTableGUI);
  ui->textBrowser->append("Updated GUI Table No." + QString::number(domainGUIIndex));

}

void MainWindow::updatenorthBoundaryTableGUI(BoundaryStaggered b) {
  ui->tableView->setModel(nullptr);
  domainTableGUI.setRowCount(b.northWallctr);
  domainTableGUI.setColumnCount(4);
  QStringList horizontalheaders;
  horizontalheaders << "Boundary_Marker" << "X_Min" << "X_Max" << "Velocity\nPressure";
  domainTableGUI.setHorizontalHeaderLabels(horizontalheaders);
  QStringList verticalheaders;
  for (int i=0; i<b.northWallctr; i++) {
      verticalheaders << QString::number(i);
  }
  domainTableGUI.setVerticalHeaderLabels(verticalheaders);
  QStandardItem *value00;
  for (int i=0; i<b.northWallctr; i++) {
      for (int j=0; j<4; j++) {
          value00 = new QStandardItem;
          value00->setText(QString::number(b.northWall(i, j)));
          domainTableGUI.setItem(i, j, value00);
      }
    }




  ui->tableView->setModel(&domainTableGUI);
  ui->textBrowser->append("Updated GUI Table No." + QString::number(domainGUIIndex));

}

void MainWindow::updateeastBoundaryTableGUI(BoundaryStaggered b) {
  ui->tableView->setModel(nullptr);
  domainTableGUI.setRowCount(b.eastWallctr);
  domainTableGUI.setColumnCount(4);
  QStringList horizontalheaders;
  horizontalheaders << "Boundary_Marker" << "X_Min" << "X_Max" << "Velocity\nPressure";
  domainTableGUI.setHorizontalHeaderLabels(horizontalheaders);
  QStringList verticalheaders;
  for (int i=0; i<b.eastWallctr; i++) {
      verticalheaders << QString::number(i);
  }
  domainTableGUI.setVerticalHeaderLabels(verticalheaders);
  QStandardItem *value00;
  for (int i=0; i<b.eastWallctr; i++) {
      for (int j=0; j<4; j++) {
          value00 = new QStandardItem;
          value00->setText(QString::number(b.eastWall(i, j)));
          domainTableGUI.setItem(i, j, value00);
      }
    }




  ui->tableView->setModel(&domainTableGUI);
  ui->textBrowser->append("Updated GUI Table No." + QString::number(domainGUIIndex));

}

void MainWindow::updateMeshLinkTableGUI(MeshStaggered &m) {
  ui->tableView->setModel(nullptr);
  domainTableGUI.setRowCount(2);
  domainTableGUI.setColumnCount(1);
  QStringList verticalheaders;
  QStringList horizontalheaders;
  verticalheaders << "Domain" << "Boundary";
  horizontalheaders << "Mesh No. ";
  domainTableGUI.setVerticalHeaderLabels(verticalheaders);
  domainTableGUI.setHorizontalHeaderLabels(horizontalheaders);
  int domainsize = Dlist.size();
  int boundarysize = BList.size();

  ui->tableView->setModel(&domainTableGUI);


  QComboBox *comboBox = new QComboBox(this);
  for (int i=0;i<domainsize;i++) {
    comboBox->addItem("Domain_" + QString::number(i));
  }

  ui->tableView->setIndexWidget(domainTableGUI.index(0, 0), comboBox);

  QComboBox *comboBox2 = new QComboBox(this);
  for (int i=0;i<boundarysize;i++) {
    comboBox2->addItem("Boundary_" + QString::number(i));
  }

  ui->tableView->setIndexWidget(domainTableGUI.index(1, 0), comboBox2);
  comboBox->setCurrentIndex(m.meshdomainindex);
  comboBox2->setCurrentIndex(m.meshboundaryindex);
  connect(comboBox, &QComboBox::currentIndexChanged, this, [this, &m](int index) {saveCurrentIndex1(m, index);});
  connect(comboBox2, &QComboBox::currentIndexChanged, this, [this, &m](int index) {saveCurrentIndex2(m, index);});


  ui->textBrowser->append("Updated Mesh GUI Table No." + QString::number(meshGUIIndex));
}

void MainWindow::prepareDefaultTree() {
  QStandardItem *parent = defaultTree.invisibleRootItem();
  QStandardItem *geometryNode = new QStandardItem;
  geometryNode->setText("Geometry");
  parent->appendRow(geometryNode);

  QStandardItem *boundaryNode = new QStandardItem;
  boundaryNode->setText("Boundary");
  parent->appendRow(boundaryNode);

  QStandardItem *meshNode = new QStandardItem;
  meshNode->setText("Mesh");
  parent->appendRow(meshNode);

  ui->treeView->setModel(&defaultTree);

}

void MainWindow::createDomainNode() {

  QStandardItem *itm = new QStandardItem;
  DomainStaggered *table = new DomainStaggered;

  updateDomainTableGUI(*table);
  setDomainTable(*table);

  QStandardItem *root = defaultTree.invisibleRootItem();

  itm->setText("Domain_" + QString::number(domainNo));
  root->child(0, 0)->appendRow(itm);
  domainNo++;
  domainList *XX = new domainList;
  XX->matrix = table;
  XX->treeNode = itm;
  Dlist.push_back(*XX);
}

void MainWindow::createBoundaryNode() {

  QStandardItem *itm = new QStandardItem;
  BoundaryStaggered *table = new BoundaryStaggered();


  QStandardItem *root = defaultTree.invisibleRootItem();

  itm->setText("Boundary_" + QString::number(domainNo));
  root->child(1, 0)->appendRow(itm);
  QStandardItem *child_1 = new QStandardItem;
  child_1->setText(("southWall_Table"));
  QStandardItem *child_2 = new QStandardItem;
  child_2->setText(("westWall_Table"));

  QStandardItem *child_3 = new QStandardItem;
  child_3->setText(("northWall_Table"));

  QStandardItem *child_4 = new QStandardItem;
  child_4->setText(("eastWall_Table"));

  itm->appendRow(child_1);
  itm->appendRow(child_2);
  itm->appendRow(child_3);
  itm->appendRow(child_4);

  boundaryList *XX = new boundaryList;
  XX->matrix = table;
  XX->treeNode = itm;
  BList.push_back(*XX);
}

void MainWindow::createMesh() {

  QStandardItem *itm = new QStandardItem;
  MeshStaggered *table = new MeshStaggered;


  updateMeshLinkTableGUI(*table);

  QStandardItem *root = defaultTree.invisibleRootItem();

  itm->setText("Mesh_" + QString::number(meshNo));
  root->child(2, 0)->appendRow(itm);
  meshNo++;
  meshList *XX = new meshList;
  XX->matrix = table;
  XX->treeNode = itm;
  MList.push_back(*XX);
}

void MainWindow::setDomainTable(DomainStaggered &d)
{

  QStandardItem *itm =  domainTableGUI.item(0, 0);
  d.xMin = itm->text().toFloat();
  itm =  domainTableGUI.item(1, 0);
  d.xMax = itm->text().toFloat();
  itm =  domainTableGUI.item(2, 0);
  d.yMin = itm->text().toFloat();
  itm =  domainTableGUI.item(3, 0);
  d.yMax = itm->text().toFloat();
  itm =  domainTableGUI.item(4, 0);
  d.dx_temp = itm->text().toFloat();
  itm =  domainTableGUI.item(5, 0);
  d.dy_temp = itm->text().toFloat();
  d.calcDomain();
  ui->textBrowser->append("BackEnd Updated");
}

void MainWindow::setsouthBoundaryTable(BoundaryStaggered &b)
{

  QStandardItem *itm;

  for (int i=0; i<b.southWallctr; i++) {
      for (int j=0; j<4; j++) {
          itm =  domainTableGUI.item(i, j);
          b.southWall(i, j) = itm->text().toFloat();
        }
    }

  ui->textBrowser->append("BackEnd Updated");
}

void MainWindow::setwestBoundaryTable(BoundaryStaggered &b)
{

  QStandardItem *itm;

  for (int i=0; i<b.westWallctr; i++) {
      for (int j=0; j<4; j++) {
          itm =  domainTableGUI.item(i, j);
          b.westWall(i, j) = itm->text().toFloat();
        }
    }

  ui->textBrowser->append("BackEnd Updated");
}

void MainWindow::setnorthBoundaryTable(BoundaryStaggered &b)
{

  QStandardItem *itm;

  for (int i=0; i<b.northWallctr; i++) {
      for (int j=0; j<4; j++) {
          itm =  domainTableGUI.item(i, j);
          b.northWall(i, j) = itm->text().toFloat();
        }
    }

  ui->textBrowser->append("BackEnd Updated");
}

void MainWindow::seteastBoundaryTable(BoundaryStaggered &b)
{

  QStandardItem *itm;

  for (int i=0; i<b.eastWallctr; i++) {
      for (int j=0; j<4; j++) {
          itm =  domainTableGUI.item(i, j);
          b.eastWall(i, j) = itm->text().toFloat();
        }
    }

  ui->textBrowser->append("BackEnd Updated");
}

void MainWindow::on_domainsave_released()
{

    auto iter1 = Dlist.begin();
    std::advance(iter1, domainGUIIndex);
    setDomainTable(*iter1->matrix);
    ui->textBrowser->append("Saved Domain Table" + QString::number(domainGUIIndex));
}

void MainWindow::on_pushButton_clicked()
{
  createDomainNode();
  createBoundaryNode();

}

void MainWindow::on_treeView_clicked(const QModelIndex &index)
{
  if (defaultTree.itemFromIndex(index)->parent()) {
    QStandardItem *itm = defaultTree.itemFromIndex(index)->parent();
    if(itm->text()== "Geometry") {
      int row = index.row();
      auto iter = Dlist.begin();
      std::advance(iter, row);
      ui->textBrowser->append("TreeView Clicked"+QString::number(row));
      domainGUIIndex = row;
      updateDomainTableGUI(*iter->matrix);
      QTextStream ms1;
      stringstream ms2;
      string str;
      QString qstr;
      ms2 << iter->matrix->D_boundaryMarker;
      while(getline(ms2, str, '\n')) {
        qstr = QString::fromStdString(str);
        ui->textBrowser->append(qstr);
      }
    }
    if(itm->text()== "Mesh") {
      int row = index.row();
      meshNo = row;
      auto iter = MList.begin();
      std::advance(iter, row);
      ui->textBrowser->append("TreeView Clicked"+QString::number(row));
      meshGUIIndex = row;
      updateMeshLinkTableGUI(*iter->matrix);
      QTextStream ms1;
      stringstream ms2;
      string str;
      QString qstr;
      ms2 << iter->matrix->U_Matrix;
      while(getline(ms2, str, '\n')) {
        qstr = QString::fromStdString(str);
        ui->textBrowser->append(qstr);
      }
      ui->textBrowser->append("\n\n\n");
      ms2 << iter->matrix->V_Matrix;
      while(getline(ms2, str, '\n')) {
        qstr = QString::fromStdString(str);
        ui->textBrowser->append(qstr);
      }
      ui->textBrowser->append("\n\n\n");
      ms2 << iter->matrix->P_Matrix;
      while(getline(ms2, str, '\n')) {
        qstr = QString::fromStdString(str);
        ui->textBrowser->append(qstr);
      }
      ui->textBrowser->append("\n\n\n");
    }


    if (defaultTree.itemFromIndex(index)->parent()->parent()) {
      QStandardItem *itm = defaultTree.itemFromIndex(index)->parent();

      if(itm->parent()->text()== "Boundary") {
          int row_parent = defaultTree.itemFromIndex(index)->parent()->row();
          boundaryGUIIndex = row_parent;
          int row = index.row();
          edgeGUIIndex = row;
          auto iter = BList.begin();
          std::advance(iter, row_parent);
          ui->textBrowser->append("TreeView Clicked"+QString::number(row));
          QTextStream ms1;
          stringstream ms2;
          string str;
          QString qstr;
          if(row==0) {
            updatesouthBoundaryTableGUI(*iter->matrix);
            ms2 << iter->matrix->southWall;
            }
          if(row==1) {
            updatewestBoundaryTableGUI(*iter->matrix);
            ms2 << iter->matrix->westWall;
            }
          if(row==2) {
            updatenorthBoundaryTableGUI(*iter->matrix);
            ms2 << iter->matrix->northWall;
            }
          if(row==3) {
            updateeastBoundaryTableGUI(*iter->matrix);
            ms2 << iter->matrix->eastWall;
            }


          while(getline(ms2, str, '\n')) {
            qstr = QString::fromStdString(str);
            ui->textBrowser->append(qstr);
          }
      }

    }

    }

}

void MainWindow::on_BoundaryMenu(const QPoint &pos){
  QMenu *contextMenu = new QMenu(ui->treeView);
  QAction *newAction = new QAction("Add_Instance", this);
  contextMenu->addAction(newAction);
  contextMenu->exec(ui->treeView->viewport()->mapToGlobal(pos));
  QModelIndex index = ui->treeView->indexAt(pos);
  addInstanceinBoundaryWall(index);
}

void MainWindow::addInstanceinBoundaryWall(QModelIndex &index){

  if (defaultTree.itemFromIndex(index)->parent()->parent()->text()=="Boundary") {
      int row_boundary = index.parent().row();
      int row = index.row();
      auto iter = BList.begin();
      std::advance(iter, row_boundary);
      if (row==0) {
        iter->matrix->southWallctr++;
        iter->matrix->southWall.conservativeResize(iter->matrix->southWallctr, 4);
        for (int i=0; i<4; i++) {
          iter->matrix->southWall(iter->matrix->southWallctr-1, i) = 0;

          }
        updatesouthBoundaryTableGUI(*iter->matrix);
        QString s = "Updated South Wall " + QString::number(iter->matrix->southWallctr) +  " of Boundary Domain " + iter->treeNode->text();
        ui->textBrowser->append(s);
        }
      if (row==1) {
        iter->matrix->westWallctr++;
        iter->matrix->westWall.conservativeResize(iter->matrix->westWallctr, 4);
        for (int i=0; i<4; i++) {
          iter->matrix->westWall(iter->matrix->westWallctr-1, i) = 0;
          }
        updatewestBoundaryTableGUI(*iter->matrix);
        QString s = "Updated West Wall " + QString::number(iter->matrix->westWallctr) +  " of Boundary Domain " + iter->treeNode->text();
        ui->textBrowser->append(s);
        }
      if (row==2) {
        iter->matrix->northWallctr++;
        iter->matrix->northWall.resize(iter->matrix->northWallctr, 4);
        for (int i=0; i<4; i++){
          iter->matrix->northWall(iter->matrix->northWallctr-1, i) = 0;
        }
        updatenorthBoundaryTableGUI(*iter->matrix);
        QString s = "Updated North Wall " + QString::number(iter->matrix->northWallctr) +  " of Boundary Domain " + iter->treeNode->text();
        ui->textBrowser->append(s);
        }
      if (row==3) {
        iter->matrix->eastWallctr++;
        iter->matrix->eastWall.resize(iter->matrix->eastWallctr, 4);
        for (int i=0; i<4; i++) {
          iter->matrix->eastWall(iter->matrix->eastWallctr-1, i) = 0;
          }
        updateeastBoundaryTableGUI(*iter->matrix);
        QString s = "Updated East Wall " + QString::number(iter->matrix->eastWallctr) +  " of Boundary Domain " + iter->treeNode->text();
        ui->textBrowser->append(s);
        }
  }
}

void MainWindow::on_boundarysave_clicked()
{
  auto iter2 = BList.begin();
  std::advance(iter2, boundaryGUIIndex);
  if (edgeGUIIndex==0){
    setsouthBoundaryTable(*iter2->matrix);
    ui->textBrowser->append("Saved Boundary Table" + QString::number(boundaryGUIIndex) + "Edge" + QString::number(edgeGUIIndex));
  }
  if (edgeGUIIndex==1) {
    setwestBoundaryTable(*iter2->matrix);
    ui->textBrowser->append("Saved Boundary Table" + QString::number(boundaryGUIIndex) + "Edge" + QString::number(edgeGUIIndex));
  }
  if (edgeGUIIndex==2) {
    setnorthBoundaryTable(*iter2->matrix);
    ui->textBrowser->append("Saved Boundary Table" + QString::number(boundaryGUIIndex) + "Edge" + QString::number(edgeGUIIndex));
  }
  if (edgeGUIIndex==3) {
    seteastBoundaryTable(*iter2->matrix);
    ui->textBrowser->append("Saved Boundary Table" + QString::number(boundaryGUIIndex) + "Edge" + QString::number(edgeGUIIndex));
  }
  auto iter3 = MList.begin();
  std::advance(iter3, 0);


}

void MainWindow::on_pushButton_4_clicked()
{
    createMesh();
}

void MainWindow::saveCurrentIndex1(MeshStaggered &m, int index)
{
    m.meshdomainindex = index;
    ui->textBrowser->append("set boundary at Row_1 to " + QString::number(index));

}

void MainWindow::saveCurrentIndex2(MeshStaggered &m, int index)
{
    m.meshboundaryindex = index;
    ui->textBrowser->append("set boundary at Row_2 to " + QString::number(index));

}

void MainWindow::on_btnMeshgeneration_clicked()
{
  auto iter1 = MList.begin();
  std::advance(iter1, meshNo);
  int dno = iter1->matrix->meshdomainindex;
  int bno = iter1->matrix->meshboundaryindex;

  auto iter2 = Dlist.begin();
  std::advance(iter2, dno);
  auto iter3 = BList.begin();
  std::advance(iter3, bno);

  float xMin, xMax, yMin, yMax;
  float dx = iter2->matrix->dx;
  float dy = iter2->matrix->dy;
  int nx = iter2->matrix->nx;
  int ny = iter2->matrix->ny;
  float D_marker;
  float U_marker;
  float V_marker;
  float P_marker;
  for(int ctr=0; ctr<iter3->matrix->southWallctr; ctr++) {
    if(iter3->matrix->southWallctr==1) {
        D_marker = iter3->matrix->southWall(0);
        if (D_marker==1)
          U_marker = iter3->matrix->southWall(3);
        if (D_marker==-1)
          P_marker = iter3->matrix->southWall(3);
        xMin = iter3->matrix->southWall(1);
        xMax = iter3->matrix->southWall(2);
      }
    else {
      D_marker = iter3->matrix->southWall(ctr, 0);
      if (D_marker==1)
        U_marker = iter3->matrix->southWall(ctr, 3);
      if (D_marker==-1)
        P_marker = iter3->matrix->southWall(ctr, 3);

      xMin = iter3->matrix->southWall(ctr, 1);
      xMax = iter3->matrix->southWall(ctr, 2);
    }
    int Xmin_iter = (int)((xMin - iter2->matrix->xMin)/dx);
    int Xmax_iter = (int)((xMax - iter2->matrix->xMin)/dx);



    for (int i=Xmin_iter; i<=Xmax_iter; i++) {
      iter1->matrix->D_Marker(ny, i) = D_marker;
      if(D_marker==1)
        iter1->matrix->U_Matrix(ny, i) = U_marker;
      if(D_marker==1)
        iter1->matrix->P_Matrix(ny, i) = P_marker;
    }
  }
  for(int ctr=0; ctr<iter3->matrix->northWallctr; ctr++) {
    if (iter3->matrix->northWallctr==1) {
        D_marker = iter3->matrix->northWall(0);
        if(D_marker==1)
          U_marker = iter3->matrix->northWall(3);
        if (D_marker==-1)
          P_marker = iter3->matrix->northWall(3);

        xMin = iter3->matrix->northWall(1);
        xMax = iter3->matrix->northWall(2);

      }
    else {
      D_marker = iter3->matrix->northWall(ctr, 0);
      if (D_marker==1)
        U_marker = iter3->matrix->northWall(ctr, 3);
      if (D_marker==-1)
        P_marker = iter3->matrix->northWall(ctr, 3);

      xMin = iter3->matrix->northWall(ctr, 1);
      xMax = iter3->matrix->northWall(ctr, 2);
      }
    int Xmin_iter = (int)((xMin - iter2->matrix->xMin)/dx);
    int Xmax_iter = (int)((xMax - iter2->matrix->xMin)/dx);



    for (int i=Xmin_iter; i<=Xmax_iter; i++) {
      iter1->matrix->D_Marker(0, i) = D_marker;
      if(D_marker==1)
        iter1->matrix->U_Matrix(0, i) = U_marker;
      if(D_marker==-1)
        iter1->matrix->P_Matrix(0, i) = P_marker;
    }
  }
  for(int ctr=0; ctr<iter3->matrix->westWallctr; ctr++) {
    if(iter3->matrix->westWallctr==1) {
        D_marker = iter3->matrix->westWall(0);
        if (D_marker==1)
          V_marker = iter3->matrix->westWall(3);
        if (D_marker==-1)
          P_marker = iter3->matrix->westWall(3);

        yMin = iter3->matrix->westWall(1);
        yMax = iter3->matrix->westWall(2);

      }
    else {
      D_marker = iter3->matrix->westWall(ctr, 0);
      if (D_marker==1)
        V_marker = iter3->matrix->westWall(ctr, 3);
      if(D_marker==-1)
        P_marker = iter3->matrix->westWall(ctr, 3);

      yMin = iter3->matrix->westWall(ctr, 1);
      yMax = iter3->matrix->westWall(ctr, 2);
      }
    int Ymin_iter = (int)((iter2->matrix->yMax - yMax)/dy);
    int Ymax_iter = (int)((iter2->matrix->yMax - yMin)/dx);



    for (int i=Ymin_iter; i<=Ymax_iter; i++) {
      iter1->matrix->D_Marker(i, 0) = D_marker;
      if(D_marker==1)
        iter1->matrix->U_Matrix(i, 0) = V_marker;
      if (D_marker==-1)
        iter1->matrix->P_Matrix(i, 0) = P_marker;
    }
  }
  for(int ctr=0; ctr<iter3->matrix->eastWallctr; ctr++) {
    if(iter3->matrix->eastWallctr==1) {
        D_marker = iter3->matrix->eastWall(0);
        if(D_marker==1)
          V_marker = iter3->matrix->eastWall(3);
        if(D_marker==-1)
          P_marker = iter3->matrix->eastWall(3);

        yMin = iter3->matrix->eastWall(1);
        yMax = iter3->matrix->eastWall(2);

      }
    else {
      D_marker = iter3->matrix->eastWall(ctr, 0);
      if(D_marker==1)
        V_marker = iter3->matrix->eastWall(ctr, 3);
      if(D_marker==-1)
        P_marker = iter3->matrix->eastWall(ctr, 3);

      yMin = iter3->matrix->eastWall(ctr, 1);
      yMax = iter3->matrix->eastWall(ctr, 2);
      }
    int Ymin_iter = (int)((iter2->matrix->yMax - yMax)/dy);
    int Ymax_iter = (int)((iter2->matrix->yMax - yMin)/dx);



    for (int i=Ymin_iter; i<=Ymax_iter; i++) {
      iter1->matrix->D_Marker(i, nx) = D_marker;
      if (D_marker==1)
        iter1->matrix->U_Matrix(i, nx) = V_marker;
      if (D_marker==-1)
        iter1->matrix->P_Matrix(i, nx) = P_marker;
    }
  }
}

