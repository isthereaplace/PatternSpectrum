# pragma once

/*Модуль содержит описание графа смежности сайтов многоугольной
  фигуры в виде обобщённой триангуляции Делоне, состоящей из
  вершин, рёбер и граней.*/


#include "LinkedList.h"
#include "ContourTracer.h"
#include <QDebug>

#include <iostream>
using std::cout;
using std::endl;


class Couple;
class Domain;
class Edge;
class Element;
class Point;
class TContour;
class TDisc;
class TriangleMap;
class Triplet;
class Vertex;


class Point : public LinkedListElement < Point >    /* Точка на плоскости */ {
	friend class Couple;
	friend class Domain;
	friend class Edge;
	friend class Element;
	friend class TContour;
	friend class TDisc;
	friend class TriangleMap;
	friend class Triplet;
	friend class Vertex;
public:
	double X, Y;  /* Координаты */
public:
	int Number;  /*Порядковый номер*/
	Point(double X, double Y);
	Point();
};

/*Здесь термин Элемент используется в смысле Сайт - старая терминология*/


class Element : public LinkedListElement < Element >  /* Элемент коллекции (сайт) - вершина или ребро */ {
	friend class Couple;
	friend class Domain;
	friend class Edge;
	friend class Point;
	friend class TContour;
	friend class TDisc;
	friend class TriangleMap;
	friend class Triplet;
	friend class Vertex;
public:
	bool isVertex;
	int NEl;   /*Номер*/
	TContour* Cont;  /*Контур, в котором находится сайт*/
	bool IsVertex();
};


class TDisc : public LinkedListElement < TDisc >     /* Окружность */ {
	friend class Couple;
	friend class Domain;
	friend class Edge;
	friend class Element;
	friend class Point;
	friend class TContour;
	friend class TriangleMap;
	friend class Triplet;
	friend class Vertex;
public:
	double X, Y, Rad;      /* Координаты центра и радиус */
	bool HalfPlane;   /* Признак вырожденного круга - полуплоскости*/
	TDisc(double XC, double YC, double Radius);
	TDisc(double a, double B, double C, bool halfPlane);
    ~TDisc() {}
    TDisc() {}
};


class Edge : public Element  /* Сайт-сегмент */ {
	friend class Couple;
	friend class Domain;
	friend class Element;
	friend class Point;
	friend class TContour;
	friend class TDisc;
	friend class TriangleMap;
	friend class Triplet;
	friend class Vertex;
public:
	Point* org, *dest;      /* Точки - концы сегмента */
	Edge(Point* P1, Point* P2);
	bool WestDirect();
	virtual ~Edge()
		/*западная направленность*/;
};


class Vertex : public Element/* Сайт-точка */ {
	friend class Couple;
	friend class Domain;
	friend class Edge;
	friend class Element;
	friend class Point;
	friend class TContour;
	friend class TDisc;
	friend class TriangleMap;
	friend class Triplet;
public:
	Point* p;             /* Точка */
	bool Reflex; /*TRUE - угол при вершине больше 180 град*/
	Vertex(Point* OwnPoint);
	bool ConformWithDisc(TDisc* Cr);
	/* Согласованность с диском: он не "за спиной" */
};


enum PlanePosition {
	LEFT_POS,
	RIGHT_POS,
	Beyond,
	Behind,
	Between,
	Origin,
	Destination,
	Cross
};

//  TDisc = Disc;


class Couple : public LinkedListElement < Couple >   /* Пара - дуга в триангуляции */ {
	friend class Domain;
	friend class Edge;
	friend class Element;
	friend class Point;
	friend class TContour;
	friend class TDisc;
	friend class TriangleMap;
	friend class Triplet;
	friend class Vertex;
public:
	Element* Frst, *Scnd;   /* Первый и второй элементы */
	Triplet* right;       /* Правый симплекс */
	Couple* PasteCouple;  /*Пара для склейки - рабочая переменная*/
	Couple(Element* El1, Element* El2);
	bool EquCoples(Couple* p); /*TRUE - вершины совпадают*/
};


class Triplet : public LinkedListElement < Triplet >  /* Тройка - симплекс в триангуляции */ {
	friend class Couple;
	friend class Domain;
	friend class Edge;
	friend class Element;
	friend class Point;
	friend class TContour;
	friend class TDisc;
	friend class TriangleMap;
	friend class Vertex;
public:
	int Numb;
	Element* E1, *E2, *e3;    /* Элементы - вершины симплекса против ЧС */
	Triplet* t1, *t2, *t3;    /* Смежные симплексы против ЧС: t1=>e1=>t2=>e2=>t3=>e3 */
	TDisc* Circ;           /* Пустой круг симплекса */
	double Depth;        /* Глубина центра (для очистки скелета)*/
	Couple* ExtCouple;    /* Смежная пара из оболочки, если такая есть*/
	Triplet(Element* El1, Element* El2, Element* El3, TDisc* Cr);
	virtual ~Triplet();
	Element* FollowingElement(Element* E);
	/* Следующий за E элемент против ЧС */
	Triplet* AdjacentForElement(Element* E);
	/* Смежный симплекс против ЧС относительно элемента E */
	Element* AdjacentForTriplet(Triplet* t);
	/* Смежный элемент против ЧС относительно симплекса T */
	void PasteWithTriplet(Element* E, Triplet* t);
	/* Приклеить симплекс T следом за элементом E */
	void BreakAdjacence(Triplet* t);
	/* Разорвать связь с симплексом t */
	bool ConsistsElements(Element* El1, Element* El2);
	/* Инцидентность дуге El1-El2 против ЧС */
};


class TriangleMap : public LinkedListElement < TriangleMap > /* Карта смежности коллекции точек и сегментов */ {
	friend class Couple;
	friend class Domain;
	friend class Edge;
	friend class Element;
	friend class Point;
	friend class TContour;
	friend class TDisc;
	friend class Triplet;
	friend class Vertex;
public:
	LinkedListTail<Couple>* MapHull;        /* Список внешних дуг карты */
	LinkedListTail<Triplet>* MapTriplet;    /* Список граней-симплексов карты */
	int NSite;      /*Количество сайтов в карте*/
	TriangleMap();
	virtual ~TriangleMap();
};

typedef Element* AdArrElem;

class TContour : public LinkedListElement < TContour >  /* Замкнутый контур - граница области. Область слева. */ {
	friend class Couple;
	friend class Domain;
	friend class Edge;
	friend class Element;
	friend class Point;
	friend class TDisc;
	friend class TriangleMap;
	friend class Triplet;
	friend class Vertex;
public:
	int Number;  /*Порядковый номер*/
	LinkedListTail<Point>* ListPoints;     /* Список точек */
	LinkedListTail<Element>* ListElements;   /* Список элементов*/
	AdArrElem* Elements;   /* Массив элементов, составляющих контур */
	int NumbElem;    /* Количество элементов */
	bool Internal; /* TRUE - внутренний (по ЧС), FALSE - внешний (против ЧС)*/
	bool ClosedPath; /* Признак замкнутого контура */
	TContour* Container;  /* Контур, объемлющий SELF */
	LinkedListTail<TContour>* MySons;           /* Контура, содержащиеся внутри SELF */
	TriangleMap* Map;    /* Карта триангуляции контура */
	Element* WestElement; /*самая левая точка в контуре*/
	Element* ClosestSite;   /*ближайший сайт слева вне контура*/
	bool Fiction;   /*Фиктивный внешний контур для внешнего скелета*/
	TContour();
	virtual ~TContour();
	void AddPoint(double X, double Y);
	void ShiftHead();  /* Установка самой левой точки в начало */
	bool ConterClockWise();   /* Контур против ЧС */
	void Invert();                   /* Сменить направление обхода */
	void CreateElements(int& Number);
	/* Построение массива элементов с номерами, начиная с Number */
};


class Domain {
	friend class Couple;
	friend class Edge;
	friend class Element;
	friend class Point;
	friend class TContour;
	friend class TDisc;
	friend class TriangleMap;
	friend class Triplet;
	friend class Vertex;
public:
	LinkedListTail<TContour>* Boundary;       /* Список замкнутых контуров - границ объекта */
	bool ElementsExist; /* TRUE - Списки элементов составлены */
	bool MapExist;      /* TRUE - Триангуляция построена */
	Domain();
	TContour* AddContour(); /* открыть новый контур в объекте */
	virtual ~Domain();
};
