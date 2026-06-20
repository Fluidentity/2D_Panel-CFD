#include <QStandardItemModel>
#include <QComboBox>

class CustomStandardItemModel : public QStandardItemModel
{
public:
    explicit CustomStandardItemModel(int rows, int columns, QObject* parent = nullptr)
        : QStandardItemModel(rows, columns, parent)
    {
    }

    QVariant data(const QModelIndex& index, int role = Qt::DisplayRole) const override
    {
        if (role == Qt::UserRole && index.column() == 0) {
            QComboBox* editor = new QComboBox();
            // Populate the combobox with items
            editor->addItem("Option 1");
            editor->addItem("Option 2");
            editor->addItem("Option 3");
            editor->setCurrentIndex(m_items.value(index.row(), 0));
            connect(editor, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index){
                m_items[index.row()] = index;
            });
            return qVariantFromValue((QWidget*)editor);
        }
        return QStandardItemModel::data(index, role);
    }

    bool setData(const QModelIndex& index, const QVariant& value, int role = Qt::EditRole) override
    {
        if (role == Qt::UserRole && index.column() == 0) {
            QComboBox* editor = qvariant_cast<QComboBox*>(value);
            m_items[index.row()] = editor->currentIndex();
            emit dataChanged(index, index);
            return true;
        }
        return QStandardItemModel::setData(index, value, role);
    }

private:
    QMap<int, int> m_items;
};
