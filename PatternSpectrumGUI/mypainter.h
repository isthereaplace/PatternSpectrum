#pragma once

#include <QWidget>
#include <QImage>
#include <QGraphicsView>
#include <QGraphicsEllipseItem>

#include "../SkeletonLib/BSTrans.h"

class MyPainter : public QGraphicsView
{
	Q_OBJECT

public:
	MyPainter(QWidget *parent);
	~MyPainter();
	void makeVisualization();

	TPolFigure *skeleton;
	bool ready;
	bool drawBones;
	bool drawCircles;
	bool drawContours;
	bool drawImage;
	double radius;
	double scope;
	bool makeLacunas;
	QGraphicsScene* scene;

private:
	QString filepath;
	QImage image;

public slots:
	void imageChanged(QString filepath_);
	
};
