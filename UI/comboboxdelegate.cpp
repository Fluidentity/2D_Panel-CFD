#include "comboboxdelegate.h"
#include <QComboBox>
#include <QModelIndex>
#include <QApplication>
#include <QStyledItemDelegate>

ComboBoxDelegate::ComboBoxDelegate(QObject *parent)
    : QStyledItemDelegate(parent)
{
}

QWidget *ComboBoxDelegate::createEditor(QWidget *parent,
    const QStyleOptionViewItem &/* option */,
    const QModelIndex &index) const
{
    QComboBox *editor = new QComboBox(parent);
    editor->addItems(QStringList() << "Option 1" << "Option 2" << "Option 3");

    // Connect the currentIndexChanged signal of the combo box to the custom signal of the delegate
    connect(editor, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged), [this, index](int currentIndex) {
        emit currentIndexChanged(currentIndex);
    });

    return editor;
}

void ComboBoxDelegate::setEditorData(QWidget *editor,
                                    const QModelIndex &index) const
{
    int currentIndex = index.model()->data(index, Qt::EditRole).toInt();
    QComboBox *comboBox = static_cast<QComboBox*>(editor);
    comboBox->setCurrentIndex(currentIndex);
}


void ComboBoxDelegate::setModelData(QWidget *editor, QAbstractItemModel *model,
                                    const QModelIndex &index) const
{
    QComboBox *comboBox = static_cast<QComboBox*>(editor);
    int currentIndex = comboBox->currentIndex();
    model->setData(index, currentIndex, Qt::EditRole);
}

void ComboBoxDelegate::updateEditorGeometry(QWidget *editor,
    const QStyleOptionViewItem &option, const QModelIndex &/* index */) const
{
    editor->setGeometry(option.rect);
}


