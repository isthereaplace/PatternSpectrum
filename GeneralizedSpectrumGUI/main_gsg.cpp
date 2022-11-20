#include "generalizedspectrumgui.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	GeneralizedSpectrumGUI w;
	w.show();
	return a.exec();
}
