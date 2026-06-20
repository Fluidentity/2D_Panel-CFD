#include "dialog.h"
#include "ui_dialog.h"

Dialog::Dialog(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::Dialog)
{
  ui->setupUi(this);
  for (int i=0; i<5; i++)
    ui->comboBox->addItem("Domain_" + QString::number(i));
  for (int i=0; i<5; i++)
    ui->comboBox_2->addItem("Boundary_" + QString::number(i));
}

Dialog::~Dialog()
{
  delete ui;
}

void Dialog::on_buttonBox_rejected()
{
  reject();

}


void Dialog::on_buttonBox_accepted()
{
  accept();

}

