#ifndef rstools_rsbatch_execution_timecourse_h
#define rstools_rsbatch_execution_timecourse_h

#include <QApplication>
#include <QWidget>
#include <QMenuBar>
#include <QSignalMapper>
#include "ui/jobeditor.ui.h"
#include "ui/TaskWidget.h"
#include "src/batch/util/rstool.hpp"
#include "src/batch/util/rstask.hpp"
#include "src/batch/util/rsjob.hpp"
#include "src/batch/util/pluginmanager.hpp"
#include "src/batch/util/rsjobparser.hpp"

class JobEditorWindow : public QMainWindow
{
    Q_OBJECT
public:
    explicit JobEditorWindow(QMainWindow *parent = 0);
    ~JobEditorWindow();
    
    void openJob(char* job);

protected slots:
    void newFile();
    void open();
    void save();
    void insertNewTask(int taskIndex);
    
protected:
    void createActions();
    void createMenus();
    void createInsertTaskMenuItems();
    void insertTask(RSTask* task);
    
    
    Ui::JobEditor ui;
    
    QMenuBar *_menuBar;
    
    QMenu *fileMenu;
    QAction *newAct;
    QAction *openAct;
    QAction *saveAct;
    QAction *exitAct;
    
    QMenu *insertMenu;
    
    RSJob *currentJob;
};

#endif
