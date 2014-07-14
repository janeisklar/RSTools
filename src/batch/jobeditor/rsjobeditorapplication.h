#ifndef rstools_rsbatch_execution_timecourse_h
#define rstools_rsbatch_execution_timecourse_h

#include <QApplication>
#include <QWidget>
#include <QMenuBar>
#include <QSignalMapper>
#include "ui/jobeditor.ui.h"
#include "src/batch/execution/tool.hpp"

class JobEditorWindow : public QMainWindow
{
    Q_OBJECT
public:
    explicit JobEditorWindow(QMainWindow *parent = 0);
    ~JobEditorWindow();

protected slots:
    void newFile();
    void open();
    void save();
    void insertTask(int code);
    
protected:
    void createActions();
    void createMenus();
    void createInsertTaskMenuItems();
    
    Ui::JobEditor ui;
    
    QMenuBar *_menuBar;
    
    QMenu *fileMenu;
    QAction *newAct;
    QAction *openAct;
    QAction *saveAct;
    QAction *exitAct;
    
    QMenu *insertMenu;
};

#endif
