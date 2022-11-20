#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_generalizedspectrumgui.h"

#include <QFileDialog>

#include <fstream>
using std::ofstream;

#include "../SkeletonLib/BSTrans.h"
#include "../GeneralizedSpectrumConsole/GeneralizedSpectrum.h"

class GeneralizedSpectrumGUI : public QMainWindow
{
	Q_OBJECT

public:
	GeneralizedSpectrumGUI(QWidget *parent = 0);
	~GeneralizedSpectrumGUI();

private:
	Ui::GeneralizedSpectrumGUIClass ui;
	GeneralizedSpectrum* spectrum;
	QString filepath;

	void CalcSpectrum(double);
	void DrawPlot();
	
private slots:
	void openImageButtonClicked();
	void checkBoxesChanged(int);
	void radiusChanged(double);
	void scaleChanged(int);
	void updateSkeleton();
	void saveImage();
    void savePolygons();

signals:
	void newImageLoaded(QString filepath_);
};
