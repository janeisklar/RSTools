#include <QSize>
#include "ArgumentsModel.h"

namespace rstools {
namespace batch {
namespace util {

ArgumentsModel::ArgumentsModel(RSJob* job, QObject *parent) : QAbstractTableModel(parent)
{
    this->job = job;
}

ArgumentsModel::~ArgumentsModel()
{}

int ArgumentsModel::rowCount(const QModelIndex & /*parent*/) const
{
   return job->getArguments().size() + 1;
}

int ArgumentsModel::columnCount(const QModelIndex & /*parent*/) const
{
    return 2;
}

QVariant ArgumentsModel::data(const QModelIndex &index, int role) const
{    
    if (role == Qt::DisplayRole || role == Qt::EditRole ) {

        vector<rsArgument*> args = job->getArguments();
        
        if ( index.row() >= (int)args.size() ) {
            return QVariant();
        }
        
        rsArgument* arg = args[index.row()];
        
        switch ( index.column() ) {
            case 0:
                return QString(arg->key);
            case 1:
                return QString(arg->value);
        }
    }
    return QVariant();
}

bool ArgumentsModel::setData(const QModelIndex & index, const QVariant & value, int role)
{
    if (role == Qt::EditRole) {
        QString result = value.toString();
        QByteArray result2 = result.toLatin1();
        const char *result3 = result2.data(); 
        char *v = (char*)rsMalloc(sizeof(char)*(strlen(result3)+1));
        sprintf(v, "%s", result3);
        
        vector<rsArgument*> args = job->getArguments();
        
        if ( index.row() >= (int)args.size() ) {
            rsArgument* arg = (rsArgument*)rsMalloc(sizeof(rsArgument));
            arg->key = (char*)"<empty>";
            arg->value = (char*)"";
            job->addArgument(arg);
        }
        
        rsArgument* arg = args[index.row()];
        
        if ( index.column() == 0 ) {
            arg->key = v;
        } else {
            arg->value = v;
        }
        
        emit editCompleted(result);
    }
    return true;
}

QVariant ArgumentsModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (role == Qt::DisplayRole) {
        if (orientation == Qt::Horizontal) {
            switch (section) {
                case 0:
                    return QString("Key");
                case 1:
                    return QString("Value");
            }
        }
     }
     return QVariant();
}

Qt::ItemFlags ArgumentsModel::flags(const QModelIndex & /*index*/) const
{
    return Qt::ItemIsSelectable | Qt::ItemIsEditable | Qt::ItemIsEnabled ;
}

}}} // namespace rstools::batch::util
