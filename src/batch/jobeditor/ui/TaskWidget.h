#ifndef rstools_rsbatch_jobeditor_ui_taskwidget_h
#define rstools_rsbatch_jobeditor_ui_taskwidget_h

#include <QScrollArea>
#include <QTabWidget>
#include "src/batch/util/rstool.hpp"
#include "src/batch/util/pluginmanager.hpp"
#include "src/batch/util/rstask.hpp"
#include "src/batch/util/rstool.hpp"
#include "src/utils/rsui.h"
#include "SettingWidget.h"

using namespace std;
using namespace rstools::batch::util;

class TaskWidget : public QTabWidget
{
    Q_OBJECT
public:
    explicit TaskWidget(RSTool* tool, QWidget * parent = 0);
    ~TaskWidget();
    
    RSTask* getTask();
    RSTool* getTool();
    
    void setupLayout();
    
protected:
    RSTool *tool;
    SettingWidget **widgets;
    size_t nWidgets;
};

#endif
