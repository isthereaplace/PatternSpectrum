#include "StructureTD.h"

Point::Point() :
X(0), Y(0), Number(0)
{
}

Point::Point(double _X, double _Y)
	: X(_X),
	Y(_Y),
	Number(0)
{
}


Edge::Edge(Point* P1, Point* P2)
	: org(P1),
	dest(P2)
{
	isVertex = false;
}


Edge::~Edge()
{
	// todo check:  inherited::Destroy;
}


bool Edge::WestDirect()
{
	/*западная направленность*/
	return (org->X > dest->X) || ((org->X == dest->X) && (org->Y > dest->Y));
}


bool Element::IsVertex()
{
	return this->isVertex;
}


Vertex::Vertex(Point* OwnPoint)
	: p(OwnPoint)
{
	isVertex = true;
	Point* PP = NULL, *PS = NULL;
	PP = p->getPrevLooped();
	PS = p->getNextLooped();
	Reflex = ((PS == PP) || ((p->X - PP->X)*(PS->Y - p->Y) < (p->Y - PP->Y)*(PS->X - p->X)));
}


bool Vertex::ConformWithDisc(TDisc* Cr)
{
	bool result = false;
	Point* Pr = NULL, *Sc = NULL;
	double a = 0.0, B = 0.0, C = 0.0, D = 0.0; /* Согласованность с диском */
	Pr = p->getPrevLooped();
	Sc = p->getNextLooped();
	if (Cr->Rad == 0)
		result = true;
	else
	{
		a = Cr->X - p->X;
		B = Cr->Y - p->Y;
		C = (Pr->X - p->X) * a + (Pr->Y - p->Y) * B;
		D = (Sc->X - p->X) * a + (Sc->Y - p->Y) * B;
		result = (C < 0.0001) && (D < 0.0001);
	}
	return result;
}


Couple::Couple(Element* El1, Element* El2)
	: Frst(El1),
	Scnd(El2),
	right(NULL),
	PasteCouple(NULL)
{
}


bool Couple::EquCoples(Couple* p) /*TRUE - вершины совпадают*/
{
	return (p != NULL) && (Scnd == p->Frst) && (Frst == p->Scnd);
}


Triplet::Triplet(Element* El1, Element* El2, Element* El3, TDisc* Cr)
	: Numb(0),
	E1(El1),
	E2(El2),
	e3(El3),
	t1(NULL),
	t2(NULL),
	t3(NULL),
	Circ(Cr),
	Depth(0.0),
	ExtCouple(NULL)
{
}


void Triplet::BreakAdjacence(Triplet* t)
{
    if (t != NULL) {
        if (t->t1 == this) {
			t->t1 = NULL;
        } else {
			if (t->t2 == this)
				t->t2 = NULL;
			else
				if (t->t3 == this)
					t->t3 = NULL;
        }
    }

	if (t1 == t)
		t1 = NULL;
	else
		if (t2 == t)
			t2 = NULL;
		else
			if (t3 == t)
				t3 = NULL;
}


Triplet::~Triplet()
{
	if (Circ != NULL)
		delete Circ;
	BreakAdjacence(t1);
	BreakAdjacence(t2);
	BreakAdjacence(t3);
	removeFromCurrentList();
	// todo check:  inherited::Destroy;
}


Element* Triplet::FollowingElement(Element* E)
{
	Element* result = NULL;  /* Следующий за E элемент против ЧС */
	if (E == E1)
		result = E2;
	else
		if (E == E2)
			result = e3;
		else
			result = E1;
	return result;
}


Triplet* Triplet::AdjacentForElement(Element* E)
{
	Triplet* result = NULL; /* Смежный симплекс против ЧС относительно элемента E */
	if (E == E1)
		result = t2;
	else
		if (E == E2)
			result = t3;
		else
			if (E == e3)
				result = t1;
			else
				result = NULL;
	return result;
}


Element* Triplet::AdjacentForTriplet(Triplet* t)
{
	Element* result = NULL; /* Смежный элемент против ЧС относительно симплекса T */
	if (t == t1)
		result = E1;
	else
		if (t == t2)
			result = E2;
		else
			if (t == t3)
				result = e3;
			else
				result = NULL;
	return result;
}


void Triplet::PasteWithTriplet(Element* E, Triplet* t)
{ /* Приклеить симплекс T следом за элементом E */
	if (E == E1)
		t2 = t;
	else
		if (E == E2)
			t3 = t;
		else
			if (E == e3)
				t1 = t;
}


bool Triplet::ConsistsElements(Element* El1, Element* El2)
{
	/* Инцидентность дуге El1-El2 против ЧС */
	return ((E1 == El1) && (E2 == El2)) || ((E2 == El1) && (e3 == El2)) || ((e3 == El1) && (E1 == El2));
}


TDisc::TDisc(double XC, double YC, double Radius)
	: X(XC),
	Y(YC),
	Rad(Radius),
	HalfPlane(false)
{
}


TDisc::TDisc(double a, double B, double C, bool)
	: X(a),
	Y(B),
	Rad(C),
	HalfPlane(true)
{
}


TContour::TContour()
	: Number(0),
	ListPoints(new LinkedListTail<Point>),
	ListElements(new LinkedListTail<Element>),
	Elements(NULL),
	NumbElem(0),
	Internal(false),
	ClosedPath(false),
	Container(NULL),
	MySons(NULL),
	Map(NULL),
	WestElement(NULL),
	ClosestSite(NULL),
	Fiction(false)
{
}


TContour::~TContour()
{
	Point* p = NULL;
    Vertex* V = NULL;
	Edge* Ed = NULL;
	Element* E = NULL;
	if (NumbElem > 0)
	{
		for (int stop = NumbElem - 1, i = 0; i <= stop; i++)
		{
			E = Elements[i];
			if (E->IsVertex())
			{
				V = (Vertex*)E;
				delete V;
			}
			else
				if (!E->IsVertex())
				{
				Ed = (Edge*)E;
				delete Ed;
				}
				else
                    qDebug() << "Ошибка: неверный элемент";
		}
		free(Elements);
	}
	while (!ListPoints->isEmpty())
	{
		p = ListPoints->first();
		p->removeFromCurrentList();
		delete p;
	}
	delete ListPoints;
	if (!ListElements->isEmpty())
        qDebug() << "Ошибка NOT ListElements.empty";
	delete ListElements;
	if ((MySons != NULL) && !(MySons->isEmpty()))
        qDebug() << "Ошибка sons<>NIL and NOT empty";
	if (MySons != NULL)
		delete MySons;
	if (Map != NULL)
		delete Map;
	// todo check:  inherited::Destroy;
}


void TContour::AddPoint(double X, double Y)
{
	(new Point(X, Y))->moveIntoTail(ListPoints);
}


void TContour::ShiftHead()
{
	Point* PFirst = NULL, *p = NULL;  /* Установка самой левой точки в начало */
	PFirst = ListPoints->first();
	p = PFirst->getNext();
	while (p != NULL)
	{
		if ((p->X < PFirst->X) || ((p->X == PFirst->X) && (p->Y < PFirst->Y)))
			PFirst = p;
		p = p->getNext();
	}
	p = ListPoints->first();
	while (p != PFirst)
	{
		p->moveIntoTail(ListPoints);
		p = ListPoints->first();
	}
}


bool TContour::ConterClockWise() /* Контур против ЧС */
{
	bool result = false;
	double X0 = 0.0, Y0 = 0.0, A1 = 0.0, A2 = 0.0, b1 = 0.0, b2 = 0.0, V = 0.0;
	Point* p = NULL;
	p = ListPoints->first();
	X0 = p->X;
	Y0 = p->Y;
	A1 = 0;
	A2 = 0;
	V = 0;
	while (p != NULL)
	{
		b1 = p->X - X0;
		b2 = p->Y - Y0;
		V = V + A1 * b2 - A2 * b1;
		A1 = b1;
		A2 = b2;
		p = p->getNext();
	}
	result = V > 0;
	return result;
}


void TContour::Invert()
{
	Point* p = NULL, *q = NULL;
	p = ListPoints->first();
	while (p->getNext() != NULL)
	{
		q = p->getNext();
		q->moveAsPrevFor(ListPoints->first());
	}
	p->moveAsPrevFor(ListPoints->first());
}


Domain::Domain()
	: Boundary(new LinkedListTail<TContour>),
	ElementsExist(false),
	MapExist(false)
{
	// inherited::Create();
}


TContour* Domain::AddContour() /* открыть новый контур в объекте */
{
	TContour* result = NULL;
	result = new TContour;
	result->moveIntoTail(Boundary);
	return result;
}


void TContour::CreateElements(int& Number)
/* Построение массива элементов по списку точек */
{
	Point* p = NULL, *q = NULL, *PFirst = NULL;
	int n = 0;
	Vertex* V = NULL;
	Edge* Ed = NULL;
	if (NumbElem != 0)
	{
		for (int stop = NumbElem - 1, n = 0; n <= stop; n++)
			delete Elements[n];
		free(Elements);
		NumbElem = 0;
	}
	n = ListPoints->cardinal() * 2;
	Elements = (AdArrElem*)malloc(n * sizeof(AdArrElem));

	/*Образование элементов*/
	PFirst = ListPoints->first();
	p = PFirst->getNext();
	while (p != NULL)
	{
		if ((p->X < PFirst->X) || ((p->X == PFirst->X) && (p->Y < PFirst->Y)))
			PFirst = p;
		p = p->getNext();
	}
	while (true)
	{
		p = ListPoints->first();
		if (p != PFirst)
			p->moveIntoTail(ListPoints);
		else
			break;
	}
	p = ListPoints->first();
	while (p != NULL)
	{
		q = p->getNextLooped();
		V = new Vertex(p);
		V->Cont = this;
		Elements[NumbElem] = V;
		V->NEl = Number;
		Number++;
		NumbElem++;
		V->moveIntoTail(ListElements);
		Ed = new Edge(p, q);
		Ed->Cont = this;
		Ed->NEl = Number;
		Number++;
		Elements[NumbElem] = Ed;
		NumbElem++;
		Ed->moveIntoTail(ListElements);
		p = p->getNext();
	}
	WestElement = Elements[0];
}


Domain::~Domain()
{
	TContour* C = NULL;
	if (Boundary != NULL)
	{
		while (!Boundary->isEmpty())
		{
			C = Boundary->first();
			C->removeFromCurrentList();
			delete C;
		}
		delete Boundary;
	}
	// todo check:  inherited::Destroy;
}


TriangleMap::TriangleMap()
	: MapHull(new LinkedListTail<Couple>),
	MapTriplet(new LinkedListTail<Triplet>),
	NSite(0)
{
}


TriangleMap::~TriangleMap()
{
	Triplet* t = NULL;
	Couple* C = NULL;
	while (!MapTriplet->isEmpty())
	{
		t = MapTriplet->first();
		t->removeFromCurrentList();
		delete t;
	}
	delete MapTriplet;
	while (!MapHull->isEmpty())
	{
		C = MapHull->first();
		C->removeFromCurrentList();
		delete C;
	}
	delete MapHull;
	// todo check:  inherited::Destroy;
}
