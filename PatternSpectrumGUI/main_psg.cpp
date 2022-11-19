#include "patternspectrumgui.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	PatternSpectrumGUI w;
	w.show();
	return a.exec();
}
