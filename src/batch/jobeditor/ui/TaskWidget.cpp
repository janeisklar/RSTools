#include "TaskWidget.h"
#include <QPushButton>
#include <QBoxLayout>
#include <QSpacerItem>
#include <QLabel>
#include <QTextEdit>
#include <QScrollArea>
#include <QTabWidget>

TaskWidget::TaskWidget(RSTool *tool, QWidget *parent) : QWidget(parent, 0)
{
    this->tool = tool;
    setupLayout();
}

TaskWidget::~TaskWidget()
{
    
}

RSTask* TaskWidget::getTask()
{
    return tool->getTask();
}

RSTool* TaskWidget::getTool()
{
    return tool;
}

void TaskWidget::setupLayout()
{   
    QTabWidget *tabWidget = new QTabWidget();
    QWidget *mainContent = new QWidget();
    QWidget *extendedContent = new QWidget();
    QScrollArea *mainScrollArea = new QScrollArea();
    QScrollArea *extendedScrollArea = new QScrollArea();
    
    QBoxLayout *mainLayout = new QBoxLayout(QBoxLayout::TopToBottom);
    QBoxLayout *extendedLayout = new QBoxLayout(QBoxLayout::TopToBottom);
    rsUIInterface* I = tool->createUI();
    
    nWidgets = I->nOptions;
    widgets = (SettingWidget**)malloc(sizeof(SettingWidget*)*nWidgets);
    
    for ( size_t i=0; i<I->nOptions; i++ ) {
        rsUIOption* o = I->options[i];
        
        if ( ! o->showInGUI ) {
            continue;
        }
        
        SettingWidget *setting = new SettingWidget(getTask(), o);
        widgets[i] = setting;    
        
        if ( o->group == RS_UI_GROUP_EXTENDED ) {
            extendedLayout->addWidget(setting);
        } else {
            mainLayout->addWidget(setting);
        }
    }
    
    mainLayout->addStretch(1);
    mainLayout->setSpacing(20);
    extendedLayout->addStretch(1);
    extendedLayout->setSpacing(20);
    
    mainContent->setLayout(mainLayout);
    extendedContent->setLayout(extendedLayout);
    
    mainScrollArea->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    extendedScrollArea->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    
    mainScrollArea->setWidgetResizable(true);
    extendedScrollArea->setWidgetResizable(true);
    
    mainScrollArea->setWidget(mainContent);
    extendedScrollArea->setWidget(extendedContent);
    
    tabWidget->setTabPosition(QTabWidget::East);
    
    tabWidget->addTab(mainScrollArea, QString("Main Settings"));
    tabWidget->addTab(extendedScrollArea, QString("Advanced Settings"));
    
    QLabel* headline = new QLabel(getTask()->getName());
    QFont f("Arial", 12, QFont::Bold);
    headline->setFont(f);
    
    QLabel* description = new QLabel(I->gui_description==NULL ? I->description : I->gui_description);
    QFont f2("Arial", 11, QFont::Normal);
    description->setFont(f2);
    description->setWordWrap(true);
    
    QBoxLayout *taskLayout = new QBoxLayout(QBoxLayout::TopToBottom);
    taskLayout->addWidget(headline);
    taskLayout->addWidget(description);
    taskLayout->addWidget(tabWidget);
    setLayout(taskLayout);
}
