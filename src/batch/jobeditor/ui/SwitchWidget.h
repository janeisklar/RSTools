#ifndef rstools_rsbatch_jobeditor_ui_switchwidget_h
#define rstools_rsbatch_jobeditor_ui_switchwidget_h

#include <QAbstractButton>
#include "src/batch/util/rstool.hpp"
#include "src/batch/util/pluginmanager.hpp"
#include "src/batch/util/rstask.hpp"
#include "src/batch/util/rstool.hpp"
#include "src/utils/rsui.h"

using namespace rstools::batch::util;

class SwitchWidget : public QAbstractButton
{
    Q_OBJECT
public:
    explicit SwitchWidget(QWidget * parent = 0);
    ~SwitchWidget();
    
protected:
    void paintEvent(QPaintEvent* event);
};

#endif
