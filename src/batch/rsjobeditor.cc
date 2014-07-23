#include <QApplication>
#include <QWidget>
#include "jobeditor/rsjobeditorapplication.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    JobEditorWindow widget;
    if ( argc > 1 ) {
        widget.openJob(argv[1]);
    }
    widget.show();
    return app.exec();
}
