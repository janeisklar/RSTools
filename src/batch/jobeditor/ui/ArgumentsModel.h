#ifndef rstools_rsbatch_jobeditor_ui_argumentmodel_h
#define rstools_rsbatch_jobeditor_ui_argumentmodel_h

#include <QAbstractTableModel>
#include "batch/util/rsjob.hpp"

using namespace std;

namespace rstools {
namespace batch {
namespace util {

class ArgumentsModel : public QAbstractTableModel
{
    Q_OBJECT
public:
    explicit ArgumentsModel(RSJob* job, QObject *parent = 0);
    ~ArgumentsModel();
    
    int rowCount(const QModelIndex &parent = QModelIndex()) const;
    int columnCount(const QModelIndex &parent = QModelIndex()) const;
    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const;
    bool setData(const QModelIndex & index, const QVariant & value, int role = Qt::EditRole);
    QVariant headerData(int section, Qt::Orientation orientation, int role) const;
    Qt::ItemFlags flags(const QModelIndex & /*index*/) const;
    
protected:
    RSJob *job;
    
signals:
    void editCompleted(const QString &);
};

}}} // namespace rstools::batch::util

#endif
