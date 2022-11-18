#include "ContourTracer.h"

RasterPoint::RasterPoint(int _x, int _y, bool _color)
{
	x = _x;
	y = _y;
	color = _color;
	kind = 0;
}

RasterPoint::RasterPoint()
{
	x = y = 0;
	color = false;
	kind = 0;
}

Grip::Grip()
{
	left = right = NULL;
}

ClosedPath::ClosedPath()
{
	path = new LinkedListTail < RasterPoint > ;
	Area = 0;
	Internal = false;
}

ClosedPath::~ClosedPath()
{
	RasterPoint* r0 = NULL;
	if (path != NULL)
	{
		r0 = path->first();
		while (r0 != NULL)
		{
			delete r0;
			r0 = path->first();
		}
		delete path;
	}
}

RasterPoint* ClosedPath::initialPoint()
{
	return (path == NULL) ? NULL : path->first();
}

bool ContourTracer::getNextBorderPoint(int& X, int& Y)
{
	bool CurrColor, PredColor, CurrMark, PredMark;
	CurrColor = isBlack(xCurr, yCurr);
	CurrMark = getMark(xCurr, yCurr);
	while (yCurr < h)
	{
		// если стоим в начале байтика и он весь белый, можно сразу перейти на его конец
		if (!((imageRaster->BITS_IN_BYTE - 1) & xCurr) && (imageRaster->byteIsFalse(xCurr, yCurr) || imageRaster->byteIsTrue(xCurr, yCurr))) {
			xCurr = xCurr + imageRaster->BITS_IN_BYTE - 1;
			if (xCurr >= w)
			{
				xCurr = -1;
				yCurr = yCurr + 1;
			}
			CurrColor = isBlack(xCurr, yCurr);
			CurrMark = getMark(xCurr, yCurr);
		}
		if (yCurr >= h) {
			break;
		}
		// пробуем на шаг вперед
		xCurr = xCurr + 1;
		PredColor = CurrColor;
		PredMark = CurrMark;
		if (xCurr >= w)
		{
			xCurr = 0;
			yCurr = yCurr + 1;
			PredColor = false;
			PredMark = true;
		}
		if (yCurr >= h) {
			break;
		}
		CurrColor = isBlack(xCurr, yCurr);
		CurrMark = getMark(xCurr, yCurr);
		if ((CurrColor != PredColor) && !(PredMark && CurrMark)) {
			break;
		}
	}
	X = xCurr;
	Y = yCurr;
	return yCurr < h;
}

bool ContourTracer::setGrip(Grip* S)
{
	int X0, Y0;

	if (!getNextBorderPoint(X0, Y0))
	{
		return false;
	}

	RasterPoint* P1 = new RasterPoint(X0, Y0, isBlack(X0, Y0));
	RasterPoint* P2 = new RasterPoint(X0 - 1, Y0, isBlack(X0 - 1, Y0));

	S->right = P1->color ? P2 : P1;
	S->left = P1->color ? P1 : P2;

	setMark(S->right->x, S->right->y);
	setMark(S->left->x, S->left->y);

	return true;
}

RasterPoint* ContourTracer::gripStep(Grip* S)
{
	RasterPoint* p = new RasterPoint;
	int du = S->right->x - S->left->x;
	int dv = S->right->y - S->left->y;
	if ((du == 0) || (dv == 0)) {
		/* следящая пара вертикальна или горизонтальна */
		p->x = S->right->x - dv;
		p->y = S->right->y + du;
	}
	else {
		/* следящая пара диагональна */
		p->x = S->left->x + (du - dv) / 2;
		p->y = S->left->y + (du + dv) / 2;
	}
	p->color = isBlack(p->x, p->y);
	if (p->color == S->right->color) {
		S->right = p;
	}
	else {
		S->left = p;
	}
	return p;
}

int ContourTracer::scal(RasterPoint* r0, RasterPoint* r1, RasterPoint* rprob)
{
	int C = 0;
	int A1 = r1->x - r0->x;
	int A2 = r1->y - r0->y;
	int b1 = rprob->x - r0->x;
	int b2 = rprob->y - r0->y;
	A1 = -(A1 * b2 - A2 * b1);
	if (A1 > 0)
		C = 1;
	else
		if (A1 < 0)
			C = -1;
		else
			C = 0;
	return C;
}

int ContourTracer::traceBoundaryCorridor(ClosedPath* CP, Grip* S)
{
	S->right->moveIntoTail(CP->path);
	if (S->right->x > S->left->x)
		S->left->moveIntoTail(CP->path);
	else
		S->left->moveAsPrevFor(S->right);
	RasterPoint* r = CP->path->first();
	bool Color = r->color;
	int n = 1;
	r = S->right;
	int u = S->right->x + S->left->x;
	int V = S->right->y + S->left->y;
	while (r != NULL)
	{
		r = gripStep(S);
		if (r != NULL)
		{
			int u1 = S->right->x + S->left->x;
			int V1 = S->right->y + S->left->y;
			if ((u != u1) || (V != V1))
			{
				r->moveIntoTail(CP->path);
				setMark(r->x, r->y);
				if (r->color == Color)
					n++;
			}
			else
			{
				delete r;
				r = NULL;
			}
		}
	}
	return n;
}

int ContourTracer::extractPolygonOfMinimalPerimeter(ClosedPath* CP, int n)
{
	bool B;
	RasterPoint* test, *FirstPoint, *left, *right, *Corner, *rt;
	int i = 0;
	if (n > 3)
	{
		test = CP->path->first();
		FirstPoint = test;
		Corner = test;
		do
		{
			{
				Corner->kind = 1;
				left = Corner->getNextLooped();
				right = left;
				while (!left->color)  /* ищем первую внутреннюю за Corner'ом */
					left = left->getNextLooped();
				while (right->color)  /* ищем первую внешнюю за Corner'ом */
					right = right->getNextLooped();
				if (right->color == Corner->color)
					while (true)   /* перебираем внутренние до первой внешней */
					{
					rt = left->getNextLooped();
					if (rt->color)
						left = rt;
					else
						break;
					}
				else
					while (true)   /* перебираем внешние до первой внутренней */
					{
					rt = right->getNextLooped();
					if (!rt->color)
						right = rt;
					else
						break;
					}

				/* Теперь Right внешняя белая, Left -внутренняя черная */
				test = left;        /* test - пробная */
				B = true;
				while (B)
				{
					if (!test->color) /*test - внешняя белая точка*/
					{
						if (scal(Corner, left, test) < 0)
						{ /* test левее левого створа */
							Corner = left;
							B = false;
						}
						else /* test левее правого створа */
							if (scal(Corner, right, test) <= 0)
								right = test;
					}
					else             /*test - внутрення черная точка*/
					{
						if (scal(Corner, right, test) > 0)
						{ /* test правее правого створа */
							Corner = right;
							B = false;
						}
						else    /* test правее левого створа */
							if (scal(Corner, left, test) >= 0)
								left = test;
					}
					test = test->getNextLooped();
				}
			}
		} while (!((Corner->x == FirstPoint->x) && (Corner->y == FirstPoint->y)));
	}
	else  /* n<=3 */
	{
		test = CP->path->first();
		FirstPoint = test;
		FirstPoint->kind = 1;
		test = test->getNext();
		while (test != NULL)
		{
			if ((test->color == FirstPoint->color) && ((test->x != FirstPoint->x) || (test->y != FirstPoint->y)))
				test->kind = 1;
			test = test->getNext();
		}
	}
	Corner = CP->path->first();
	i = 0;
	while (Corner != NULL)
	{
		test = Corner->getNext();
		if (Corner->kind == 1)
			i++;
		else
			delete Corner;
		Corner = test;
	}
	return i;
}

void ContourTracer::fillPolygonArea(ClosedPath* CP)
{
	float Area = 0;
	RasterPoint* r;
	RasterPoint* p = CP->path->first();
	RasterPoint* q = CP->path->last();
	r = p;
	while (r != q)
	{
		Area += (p->x - r->x) * (q->y - r->y) - (p->y - r->y) * (q->x - r->x);
		p = r;
		r = r->getNext();
	}
	Area = abs(Area / 2);
	CP->Area = Area;
}

ContourTracer::ContourTracer(BitRaster* _imageRaster, float _areaIgnore)
{
	imageRaster = _imageRaster;
	w = imageRaster->w;
	h = imageRaster->h;

	areaIgnore = _areaIgnore;

	// starting position
	yCurr = 0;
	xCurr = -1;

	// allocate mark raster
	markRaster = new BitRaster(w, h);

	contourList = new LinkedListTail < ClosedPath > ;
}

ContourTracer::~ContourTracer()
{
	delete markRaster;

	// delete contours
	ClosedPath* CP = contourList->first();
	while (CP != NULL)
	{
		delete CP;
		CP = contourList->first();
	}
	delete contourList;
}

ClosedPath* ContourTracer::initialContour()
{
	return contourList->first();
}

void ContourTracer::traceContours()
{
	ClosedPath* CP;
	int n, M;
	Grip* S = new Grip;
	while (setGrip(S))
	{
		CP = new ClosedPath;
		CP->Internal = (S->left->x < S->right->x);
		n = traceBoundaryCorridor(CP, S);
		M = extractPolygonOfMinimalPerimeter(CP, n);
		CP->moveIntoTail(contourList);
		fillPolygonArea(CP);
		if (CP->Area < areaIgnore)
			delete CP;
	}
	delete S;
}
