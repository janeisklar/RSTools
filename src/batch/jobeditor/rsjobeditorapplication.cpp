#include "rsjobeditorapplication.h"

using namespace rstools::batch::util;

void JobEditorWindow::createActions()
{
    newAct = new QAction(tr("&New"), this);
    newAct->setShortcuts(QKeySequence::New);
    newAct->setStatusTip(tr("New"));
    connect(newAct, SIGNAL(triggered()), this, SLOT(newFile()));

    openAct = new QAction(tr("&Open..."), this);
    openAct->setShortcuts(QKeySequence::Open);
    openAct->setStatusTip(tr("Open"));
    connect(openAct, SIGNAL(triggered()), this, SLOT(open()));

    saveAct = new QAction(tr("&Save"), this);
    saveAct->setShortcuts(QKeySequence::Save);
    saveAct->setStatusTip(tr("Save"));
    connect(saveAct, SIGNAL(triggered()), this, SLOT(save()));

    exitAct = new QAction(tr("E&xit"), this);
    exitAct->setShortcuts(QKeySequence::Quit);
    exitAct->setStatusTip(tr("Exit"));
    connect(exitAct, SIGNAL(triggered()), this, SLOT(close()));
}

void JobEditorWindow::createMenus()
{
    #if defined(Q_WS_MAC)
    //_menuBar = new QMenuBar(0);
    //_menuBar = menuBar();
    //_menuBar->setNativeMenuBar(false);
    #else
    //_menuBar = menuBar();
    #endif
    
    _menuBar = menuBar();
    _menuBar->setNativeMenuBar(false);
    
    fileMenu = _menuBar->addMenu(tr("&File"));
    fileMenu->addAction(newAct);
    fileMenu->addAction(openAct);
    fileMenu->addAction(saveAct);
    fileMenu->addAction(exitAct);
    
    insertMenu = _menuBar->addMenu(tr("&Insert"));
    
    createInsertTaskMenuItems();
}

void JobEditorWindow::createInsertTaskMenuItems()
{
    // ensure that plugins are loaded
    PluginManager::getInstance().loadPlugins();
    
    QSignalMapper* signalMapper = new QSignalMapper(this);
    
    vector<const char*> tools = RSTool::getTools();
    size_t i=0;
    
    for(vector<const char*>::iterator it = tools.begin(); it != tools.end(); ++it) {

        const char* code   = (char*)*it;
        RSTask* task = RSTask::taskFactory(code);
        const char* name   = task->getName();
        
        fprintf(stdout, "Tool: %s\n", name);

        QAction *action = new QAction(tr(name), this);
        connect(action, SIGNAL(triggered()), signalMapper, SLOT(map()));
        signalMapper->setMapping(action, i);
        
        insertMenu->addAction(action);
        
        i++;
    }
    
    connect(signalMapper, SIGNAL(mapped(int)), this, SLOT(insertTask(int))) ;
}

void JobEditorWindow::newFile()
{
}

void JobEditorWindow::open()
{
}

void JobEditorWindow::save()
{
}

void JobEditorWindow::insertTask(int taskIndex)
{
    const char* code = RSTool::getTools().at(taskIndex);
    RSTask* task = RSTask::taskFactory(code);
    const char* name = task->getName();
    fprintf(stdout, "Insert tool: %s\n", name);
}

JobEditorWindow::JobEditorWindow(QMainWindow *parent) : QMainWindow(parent)
{    
    ui.setupUi(this);
    createActions();
    createMenus();
}

JobEditorWindow::~JobEditorWindow()
{
    
}