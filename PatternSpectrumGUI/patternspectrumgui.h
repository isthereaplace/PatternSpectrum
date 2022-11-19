#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_patternspectrumgui.h"

#include <QFileDialog>

#include <fstream>
using std::ofstream;

#include "../SkeletonLib/BSTrans.h"
#include "../PatternSpectrumConsole/PatternSpectrum.h"

class PatternSpectrumGUI : public QMainWindow
{
	Q_OBJECT

public:
	PatternSpectrumGUI(QWidget *parent = 0);
	~PatternSpectrumGUI();

private:
	Ui::PatternSpectrumGUIClass ui;
	PatternSpectrum* spectrum;
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
