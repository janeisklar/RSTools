#ifndef rstools_rsbatch_jobeditor_ui_taskwidget_h
#define rstools_rsbatch_jobeditor_ui_taskwidget_h

#include "batch/util/rstool.hpp"
#include "batch/util/pluginmanager.hpp"
#include "batch/util/rstask.hpp"
#include "batch/util/rstool.hpp"
#include "utils/rsui.h"
#include "SettingWidget.h"

using namespace std;
using namespace rstools::batch::util;

class TaskWidget : public QWidget
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
