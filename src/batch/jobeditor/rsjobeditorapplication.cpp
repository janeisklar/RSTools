#include "rsjobeditorapplication.h"
#include "ui/ArgumentsModel.h"
#include <QFileDialog>
#include <QErrorMessage>

using namespace rstools::batch::util;

void JobEditorWindow::createActions()
{
    newAct = new QAction(tr("&New"), this);
    newAct->setShortcuts(QKeySequence::New);
    newAct->setStatusTip(tr("New"));
    newAct->setEnabled(true);
    newAct->setAutoRepeat(false);
    addAction(newAct);
    connect(newAct, SIGNAL(triggered()), this, SLOT(newFile()));

    openAct = new QAction(tr("&Open..."), this);
    openAct->setShortcuts(QKeySequence::Open);
    openAct->setStatusTip(tr("Open"));
    openAct->setShortcutContext(Qt::ApplicationShortcut);
    openAct->setEnabled(true);
    openAct->setAutoRepeat(false);
    addAction(openAct);
    connect(openAct, SIGNAL(triggered()), this, SLOT(open()));

    saveAct = new QAction(tr("&Save"), this);
    saveAct->setShortcuts(QKeySequence::Save);
    saveAct->setStatusTip(tr("Save"));
    saveAct->setShortcutContext(Qt::ApplicationShortcut);
    saveAct->setEnabled(true);
    saveAct->setAutoRepeat(false);
    addAction(saveAct);
    connect(saveAct, SIGNAL(triggered()), this, SLOT(save()));

    exitAct = new QAction(tr("E&xit"), this);
    exitAct->setShortcuts(QKeySequence::Quit);
    exitAct->setStatusTip(tr("Exit"));
    exitAct->setShortcutContext(Qt::ApplicationShortcut);
    exitAct->setEnabled(true);
    exitAct->setAutoRepeat(false);
    addAction(exitAct);
    connect(exitAct, SIGNAL(triggered()), this, SLOT(close()));
}

void JobEditorWindow::createMenus()
{
    #pragma message "TODO: Fix Mac OS X's native menu bar in " __FILE__

    // Weird QT bug that prevents us from using the Mac's native menu bar..
    /*
    #if defined(Q_WS_MAC)
    _menuBar = new QMenuBar(0);
    #else
    _menuBar = menuBar();
    #endif
    */
    
    _menuBar = menuBar();
    //#if defined(Q_WS_MAC)
    //_menuBar->setNativeMenuBar(false);
    //#endif
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
    
    for (vector<const char*>::iterator it = tools.begin(); it != tools.end(); ++it) {

        const char* code   = (char*)*it;
        RSTask* task = RSTask::taskFactory(code);
        const char* name   = task->getName();

        QAction *action = new QAction(tr(name), this);
        connect(action, SIGNAL(triggered()), signalMapper, SLOT(map()));
        signalMapper->setMapping(action, i);
        
        insertMenu->addAction(action);
        
        i++;
    }
    
    connect(signalMapper, SIGNAL(mapped(int)), this, SLOT(insertNewTask(int))) ;
}

void JobEditorWindow::newFile()
{
}

void JobEditorWindow::open()
{
    try {
        QString fileName = QFileDialog::getOpenFileName(this, tr("Open Job"), "", tr("Job (*.job)"));
        if ( fileName != NULL ) {
            openJob(fileName.toUtf8().data());
        }
    } catch (const exception& e) {
    	QErrorMessage errorMessage;
    	errorMessage.showMessage(e.what());
    	errorMessage.exec();        
    } catch (...) {
    	QErrorMessage errorMessage;
    	errorMessage.showMessage("Unknown error while opening the job file");
    	errorMessage.exec();
    }
}

void JobEditorWindow::openJob(char* jobFile)
{
    
    RSJobParser *parser = new RSJobParser(jobFile);
    parser->parse();
    currentJob = parser->getJob();
    
    ui.argumentsTable->setModel(new ArgumentsModel(currentJob));
    ui.argumentsTable->setSortingEnabled(true);
#if QT_VERSION >= 0x050000
    ui.argumentsTable->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
#else
    ui.argumentsTable->horizontalHeader()->setResizeMode(QHeaderView::ResizeToContents);
#endif
    
    vector<RSTask*> tasks = currentJob->getTasks();
    
    for(vector<RSTask*>::iterator it = tasks.begin(); it != tasks.end(); ++it) {
        RSTask* task = (RSTask*)*it;
        insertTask(task);
    }
}

void JobEditorWindow::save()
{
}

void JobEditorWindow::insertNewTask(int taskIndex)
{
    const char* code = RSTool::getTools().at(taskIndex);
    RSTask* task = RSTask::taskFactory(code);
    const char *name = task->getName();
    char *description = (char*)malloc(sizeof(char)*(strlen(name)+1));
    sprintf(description, "%s", name);
    task->setDescription(description);
    insertTask(task);
}

void JobEditorWindow::insertTask(RSTask* task)
{
    const char* code = task->getCode();
    RSTool* tool = RSTool::toolFactory(code);
    tool->setTask(task);
    const char* name = task->getDescription();

    TaskWidget *widget  = new TaskWidget(tool, ui.pipelineWidget);
    const QString title = QString(name);
    
    ui.pipelineWidget->addPage(widget, QIcon(), title);
}

JobEditorWindow::JobEditorWindow(QMainWindow *parent) : QMainWindow(parent)
{    
    ui.setupUi(this);
    ui.pipelineWidget->removePage(0);
    
    try {
        createActions();
        createMenus();
    } catch (const exception& e) {
    	QErrorMessage errorMessage;
    	errorMessage.showMessage(e.what());
    	errorMessage.exec();        
    } catch (...) {
    	QErrorMessage errorMessage;
    	errorMessage.showMessage("Unknown error while intializing the application");
    	errorMessage.exec();
    }
}

JobEditorWindow::~JobEditorWindow()
{
    
}