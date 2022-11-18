#pragma once

/*Модуль содержит процедуры построения многоугольной фигуры
для бинарного изображения, заданного в виде матрицы из
чёрных и белых точек.*/

#include "LinkedList.h"
#include "BitRaster.h"
#include <cmath>

class RasterPoint : public LinkedListElement < RasterPoint >
{
public:
	int x, y; // coordinates
	bool color; // color
	unsigned char kind; // type of point
	RasterPoint(int _x, int _y, bool _color);
	RasterPoint();
};

// Следящий мост
class Grip
{
public:
	RasterPoint *left, *right;
	Grip();
};

class ClosedPath : public LinkedListElement < ClosedPath > 
{
public:
	LinkedListTail<RasterPoint>* path;  /* кусочно-линейное представление */
	double Area;
	bool Internal;  /*Внутренний или внешний контур*/
	ClosedPath();
    ~ClosedPath();
	RasterPoint* initialPoint();
};

class ContourTracer
{
	BitRaster* imageRaster;
	BitRaster* markRaster;
	float areaIgnore;
	int xCurr, yCurr;

	inline bool isOutsideRaster(int X, int Y) 
	{
		return (X < 0) || (X >= w) || (Y < 0) || (Y >= h);
	}

	inline bool isBlack(int X, int Y)
	{
		if (isOutsideRaster(X, Y)) {
			return false;
		} 
		else {
			return imageRaster->getBit(X, Y);
		}
	}

	inline bool getMark(int X, int Y)
	{
		if (isOutsideRaster(X, Y)) {
			return true;
		} 
		else {
			return markRaster->getBit(X, Y);
		}
	}

	inline void setMark(int X, int Y)
	{
		if (!isOutsideRaster(X, Y))
		{
			markRaster->setBit(X, Y, true);
		}
	}

	// На выходе получаем координаты точки (X, Y), у которой:
	// 1. цвет с предыдущей различен
	// 2. либо сама точка, либо предыдущая - непомеченные
	// При этом считается, что точки вне картинки все белые и помеченные
	// Возвращаемое значение функции - false, если картинка исчерпана
	bool getNextBorderPoint(int& X, int& Y);

	bool setGrip(Grip* S);

	// Шаг белой (правой) точкой 
	RasterPoint* gripStep(Grip* S);

	// Вычисление положения точки rprob относительно прямой r0 r1
	// r0-нач.точка, r1-конечная, rprob-пробная
	// результат >0 - справа, <0 - слева, =0 - на линии
	int scal(RasterPoint* r0, RasterPoint* r1, RasterPoint* rprob);

	int traceBoundaryCorridor(ClosedPath* CP, Grip* S);

	int extractPolygonOfMinimalPerimeter(ClosedPath* CP, int n);

	void fillPolygonArea(ClosedPath* CP);

public:
	int w;
	int h;
	LinkedListTail<ClosedPath>* contourList;

	ContourTracer(BitRaster* _imageRaster, float _areaIgnore);
	
	~ContourTracer();

	ClosedPath* initialContour();

	void traceContours();
};
