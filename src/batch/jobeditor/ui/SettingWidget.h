#ifndef rstools_rsbatch_jobeditor_ui_settingwidget_h
#define rstools_rsbatch_jobeditor_ui_settingwidget_h

#include <stdexcept>
#include <QGroupBox>
#include "utils/rsui.h"
#include "batch/util/rstask.hpp"

using namespace rstools::batch::util;

class SettingWidget : public QGroupBox
{
    Q_OBJECT
public:
    explicit SettingWidget(RSTask* task, rsUIOption* option, QWidget * parent = 0);
    ~SettingWidget();
    
    rsUIOption* getSetting();
    
protected:
    void createValueWidget();
    void setupLayout();
    
    rsUIOption *option;
    QWidget *valueWidget;
    RSTask* task;
};

#endif
