#ifndef COMBOXDELEGATE_H
#define COMBOXDELEGATE_H

#include <QItemDelegate>
#include <QStyledItemDelegate>
#include <QComboBox>

class ComboBoxDelegate : public QStyledItemDelegate
{
    Q_OBJECT

public:
    ComboBoxDelegate(QObject* parent = nullptr);

    QWidget* createEditor(QWidget* parent, const QStyleOptionViewItem& option, const QModelIndex& index) const override;
    void setEditorData(QWidget* editor, const QModelIndex& index) const override;
    void setModelData(QWidget* editor, QAbstractItemModel* model, const QModelIndex& index) const override;
    void updateEditorGeometry(QWidget* editor, const QStyleOptionViewItem& option, const QModelIndex& index) const override;
    bool commitData(QWidget* editor);

signals:
    void currentIndexChanged(int index) const;
};

#endif // COMBOXDELEGATE_H
