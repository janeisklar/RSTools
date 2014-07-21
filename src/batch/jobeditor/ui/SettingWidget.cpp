#include "SettingWidget.h"
#include <QPushButton>
#include <QBoxLayout>
#include <QSpacerItem>
#include <QLabel>
#include <QLineEdit>
#include <QSpacerItem>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QCheckBox>
#include <glib.h>

using namespace rstools::batch::util;

SettingWidget::SettingWidget(RSTask* task, rsUIOption *option, QWidget *parent) : QGroupBox(parent)
{
    this->task   = task;
    this->option = option;
    setupLayout();
}

SettingWidget::~SettingWidget()
{
    
}

rsUIOption* SettingWidget::getSetting()
{
    return option;
}

void SettingWidget::setupLayout()
{       
    QBoxLayout *layout = new QBoxLayout(QBoxLayout::TopToBottom);
    
    setTitle(QString(option->name));

    QLabel *description = new QLabel();
    description->setText(QString(
        option->gui_description == NULL
        ? option->cli_description
        : option->gui_description 
    ));
    description->setWordWrap(true);
    layout->addWidget(description);
        
    createValueWidget();
    layout->addWidget(valueWidget);
    
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Minimum);
    setLayout(layout);
}

void SettingWidget::createValueWidget()
{
    rsArgument* argument = task->getArgument(option->name);
    
    switch(option->type) {
        case G_OPTION_ARG_FILENAME:
        case G_OPTION_ARG_STRING:
        case G_OPTION_ARG_CALLBACK:
        case G_OPTION_ARG_INT:
        case G_OPTION_ARG_INT64:
        case G_OPTION_ARG_DOUBLE:
            {
                QLineEdit *w = new QLineEdit();
                w->setPlaceholderText(option->cli_arg_description);
                if ( argument != NULL ) {
                    w->setText(argument->value);
                }
                valueWidget = w;
            }
            break;
        /*
        case G_OPTION_ARG_INT:
        case G_OPTION_ARG_INT64:
            valueWidget = new QSpinBox();
            break;
        case G_OPTION_ARG_DOUBLE:
            valueWidget = new QDoubleSpinBox();
            break;
        */
        case G_OPTION_ARG_NONE:
            {
                QCheckBox *w = new QCheckBox("Enabled"); // new SwitchWidget();
                if ( argument != NULL ) {
                    w->setCheckState(Qt::Checked);
                } else {
                    w->setCheckState(Qt::Unchecked);
                }
                valueWidget = w;
            }
            break;
        default:
            throw std::invalid_argument("UI argument type unknown");
    }
}