#include "Geometry.h"

const int MAXDOUBLE = 100000000;
const double Mnull = 1e-4;  /*Константа для машинного нуля*/
const double Mnull2 = 1e-6;  /*Константа для машинного нуля 2*/


PlanePosition Classify(Point *P1, Point *P2, Point *q)
/*Положение точки q относительно вектора p1,p2 */
{
	PlanePosition result;
	double ax = 0.0, ay = 0.0, bx = 0.0, by = 0.0, S = 0.0, C = 0.0, ma = 0.0, mb = 0.0;
	ax = P2->X - P1->X;
	ay = P2->Y - P1->Y;
	bx = q->X - P1->X;
	by = q->Y - P1->Y;
	S = ax * by - ay * bx; /* векторное произведение */
	if (S > 0)
		result = LEFT_POS;
	else
		if (S < 0)
			result = RIGHT_POS;
		else
		{
			C = ax * bx + ay * by; /* скалярное произведение */
			if (C < 0)
				result = Behind;
			else
				if (C == 0)
					result = Origin;
				else
				{
					ma = ax * ax + ay * ay; /* квадрат длины вектора */
					mb = bx * bx + by * by; /* квадрат длины вектора */
					if (ma < mb)
						result = Beyond;
					else
						if (ma == mb)
							result = Destination;
						else
							result = Between;
				}
		}
	return result;
}


bool Codirect(Point *r0, Point *r1, Point *r2)
/*монотонность трех собственных точек*/
{
	bool result = false;
	result = ((r2->X - r1->X) * (r1->X - r0->X) + (r2->Y - r1->Y) * (r1->Y - r0->Y) > 0);
	return result;
}


bool Collinear(Point *r0, Point *r1, Point *r2)
/*коллинеарность трех собственных точек*/
{
	bool result = false;
	result = ((r1->X - r0->X) * (r2->Y - r0->Y) == (r1->Y - r0->Y) * (r2->X - r0->X));
	return result;
}


void LinSyst(double A1, double b1, double C1, double A2, double b2, double C2, double* X, double* Y, double* D)
{
	double d1 = 0.0, d2 = 0.0;
	*D = A1 * b2 - A2 * b1;
    if (abs(*D) < Mnull2) {
        *D = 0.0;
    }
	else
	{
		d1 = C1 * b2 - C2 * b1;
		d2 = A1 * C2 - A2 * C1;
		*X = d1 / *D;
		*Y = d2 / *D;
	}
} /*LinSyst*/


void Circle(double r0x, double r0y, double r1x, double r1y, double r2x, double r2y, double* X, double* Y, double* Rad)
{
	double A1 = 0.0, b1 = 0.0, A2 = 0.0, b2 = 0.0, C1 = 0.0, C2 = 0.0, D = 0.0;
	A1 = r1x - r0x;
	b1 = r1y - r0y;
	A2 = r2x - r0x;
	b2 = r2y - r0y;
	C2 = r0x * r0x + r0y * r0y;
	C1 = (r1x * r1x + r1y * r1y - C2) * 0.5;
	C2 = (r2x * r2x + r2y * r2y - C2) * 0.5;
	LinSyst(A1, b1, C1, A2, b2, C2, X, Y, &D);
    if (abs(D) < Mnull2)
		*Rad = 0.0;
	else
	{
		A1 = *X - r0x;
		b1 = *Y - r0y;
		*Rad = sqrt(A1 * A1 + b1 * b1);
	}
} /*Circle*/


bool Tangent(Point *r0, Point *r1, double X, double Y)
{
	bool result = false;
	double zx = 0.0, zy = 0.0, ax = 0.0, ay = 0.0, bx = 0.0, by = 0.0;
	double S1 = 0.0, S2 = 0.0;
	bool B = false;
	/* Лежит ли проекция точки X,Y на прямую r0,r1 между r0 и r1 */
	zx = r1->X - r0->X;
	zy = r1->Y - r0->Y;
	ax = X - r0->X;
	ay = Y - r0->Y;
	bx = X - r1->X;
	by = Y - r1->Y;
	S1 = zx * ax + zy * ay;
	S2 = zx * bx + zy * by;
	B = (S1 >= -Mnull) && (S2 <= Mnull);
	result = B;
	return result;
} /*Tangent*/


bool RightOrientation(Point *r0, Point *r1, Point *r2, Point *r3, Point *r4, Point *r5)
/*Проверка правой ориентации трёх сегментов*/
{
	bool result = false;
	double V12 = 0.0, V23 = 0.0, V31 = 0.0;
	int k = 0;
	k = 0;
	V12 = (r1->X - r0->X) * (r3->Y - r2->Y) - (r1->Y - r0->Y) * (r3->X - r2->X);
	if (V12 > 0)
		k++;
	V23 = (r3->X - r2->X) * (r5->Y - r4->Y) - (r3->Y - r2->Y) * (r5->X - r4->X);
	if (V23 > 0)
		k++;
	V31 = (r5->X - r4->X) * (r1->Y - r0->Y) - (r5->Y - r4->Y) * (r1->X - r0->X);
	if (V31 > 0)
		k++;
	result = (k > 1);
	return result;
}


TDisc* NewDisc(double X, double Y, double Rad)
{
	return new TDisc(X, Y, Rad);
} /*NewDisc*/


TDisc* P3(Point *r0, Point *r1, Point *r2)
{
	TDisc* result;
	double X = 0.0, Y = 0.0, u = 0.0, V = 0.0, r = 0.0; /* Построение круга, касающегося трех точек против ЧС! */
	X = (r1->X - r0->X);
	Y = (r1->Y - r0->Y);
	u = (r2->X - r0->X);
	V = (r2->Y - r0->Y);
	if (X * V - Y * u > 0.0)
		Circle(r0->X, r0->Y, r1->X, r1->Y, r2->X, r2->Y, &X, &Y, &r);
	else
		r = 0.0;
	if (r > 0.0)
		result = NewDisc(X, Y, r);
	else
		result = NULL;
	return result;
} /*P3*/


TDisc* S1_P1(Point *r0, Point *r1, Point *r2, Point *r3)
/* Построение круга, касающегося сегмента в его конце и еще одной точки */
/* r0,r1 - сегмент, r2 и r3 - точки, одна из которых совпадает с r0 или r1 */
/* r0,r1,r2,r3 выпуклый 4-угольник, ориентированный против ЧС */
{
	TDisc* result;
	double sx = 0.0, sy = 0.0, t = 0.0, rx = 0.0, ry = 0.0, XC = 0.0, YC = 0.0;
	Point *RR2, *RR3;
	if ((r1 == r3) || (r2 == r0))
	{
		result = NULL;
		return result;
	}
	if (r2 == r1)
	{
		RR2 = r2;
		RR3 = r3;
	}
	else
	{
		RR2 = r3;
		RR3 = r2;
	}
	sx = -(r1->Y - r0->Y);
	sy = (r1->X - r0->X);
	t = sqrt(sx * sx + sy * sy);
	sx = sx / t;
	sy = sy / t;
	rx = (RR3->X - RR2->X);
	ry = (RR3->Y - RR2->Y);
	t = rx * sx + ry * sy;
	if (abs(t) < Mnull)
	{
		result = NULL;
		return result;
	}
	t = (rx * rx + ry * ry) / t / 2.0;
	if (t <= 0.0)
	{
		result = NULL;
		return result;
	}
	else
	{
		XC = (RR2->X) + sx * t;
		YC = (RR2->Y) + sy * t;
		result = NewDisc(XC, YC, t);
	}
	return result;
} /*S1_P1*/


TDisc* S3(Point *r0, Point *r1, Point *r2, Point *r3, Point *r4, Point *r5)
/* Построение круга, касающегося трех сегментов (r0,r1),(r2,r3),(r4,r5) */
{
	TDisc* result;
	double s1x = 0.0, s1y = 0.0, s2x = 0.0, s2y = 0.0, s3x = 0.0, s3y = 0.0, A1 = 0.0, A2 = 0.0, b1 = 0.0, b2 = 0.0, C1 = 0.0, C2 = 0.0, P1 = 0.0, P2 = 0.0, P3 = 0.0, D = 0.0, X = 0.0, Y = 0.0;
	if (!RightOrientation(r0, r1, r2, r3, r4, r5))
	{
		result = NULL;
		return result;
	}
	s1x = (r1->Y - r0->Y);
	s1y = -(r1->X - r0->X);
	s2x = (r3->Y - r2->Y);
	s2y = -(r3->X - r2->X);
	s3x = (r5->Y - r4->Y);
	s3y = -(r5->X - r4->X);
	D = sqrt(s1x * s1x + s1y * s1y);
	s1x = s1x / D;
	s1y = s1y / D;
	D = sqrt(s2x * s2x + s2y * s2y);
	s2x = s2x / D;
	s2y = s2y / D;
	D = sqrt(s3x * s3x + s3y * s3y);
	s3x = s3x / D;
	s3y = s3y / D;
	A1 = s1x - s2x;
	b1 = s1y - s2y;
	A2 = s1x - s3x;
	b2 = s1y - s3y;
	P1 = s1x * (r0->X) + s1y * (r0->Y);
	P2 = s2x * (r2->X) + s2y * (r2->Y);
	P3 = s3x * (r4->X) + s3y * (r4->Y);
	D = sqrt(Sqr(A1) + Sqr(b1)) + sqrt(Sqr(A2) + Sqr(b2));
	C1 = P1 - P2;
	C2 = P1 - P3;
	A1 = A1 / D;
	b1 = b1 / D;
	A2 = A2 / D;
	b2 = b2 / D;
	C1 = C1 / D;
	C2 = C2 / D;
	LinSyst(A1, b1, C1, A2, b2, C2, &X, &Y, &D);
	if (D == 0.0)
	{
		result = NULL;
		return result;
	}
	A1 = s1x * X + s1y * Y - P1;
	b1 = s2x * X + s2y * Y - P2;
	C1 = s3x * X + s3y * Y - P3;
	if ((A1 > 0.0) || (b1 > 0.0) || (C1 > 0.0))
	{
		result = NULL;
		return result;
	}
	if (Tangent(r0, r1, X, Y) && Tangent(r2, r3, X, Y) && Tangent(r4, r5, X, Y))
		result = NewDisc(X, Y, -A1);
	else
		result = NULL;
	return result;
} /*S3*/



TDisc* S2_P1(Point *r0, Point *r1, Point *r2, Point *r3, Point *r4)
/* r0,r1 и r2,r3 - сегменты, r4 точка все против ЧС! */
{
	TDisc* result;
	double e1x = 0.0, e1y = 0.0, e2x = 0.0, e2y = 0.0, ex = 0.0, ey = 0.0, sx = 0.0, sy = 0.0, t1 = 0.0, t2 = 0.0, l = 0.0, r1x = 0.0, r1y = 0.0, r2x = 0.0, r2y = 0.0, z1x = 0.0, z1y = 0.0, z2x = 0.0, z2y = 0.0, XC = 0.0, YC = 0.0, Rad = 0.0;
	e1x = (r1->X - r0->X);
	e1y = (r1->Y - r0->Y);
	e2x = (r2->X - r3->X);
	e2y = (r2->Y - r3->Y);
	if ((r0 == r3) && (r0 == r4))
	{
		if (e1x * e2y - e1y * e2x >= 0.0)
			result = NewDisc((r0->X), (r0->Y), 0.0);
		else
			result = NULL;
		return result;
	}
	t1 = sqrt(e1x * e1x + e1y * e1y);
	e1x = e1x / t1;
	e1y = e1y / t1;
	t2 = sqrt(e2x * e2x + e2y * e2y);
	e2x = e2x / t2;
	e2y = e2y / t2;
	ex = e1x + e2x;
	ey = e1y + e2y;
	t1 = sqrt(ex * ex + ey * ey);
	if (abs(t1) < Mnull)  /* r0,r1 и r2,r3 - коллинеарны */
	{
		result = NULL;
		return result;
	}
	if (abs(t1) < 0.001) /* r0,r1 и r2,r3 - почти коллинеарны*/
	{
		ex = -e1y + e2y;
		ey = e1x - e2x;
		t1 = sqrt(ex * ex + ey * ey);
	}
	sx = ey / t1;
	sy = -ex / t1;
	t2 = double(((r2->X - r4->X) * e2y - (r2->Y - r4->Y) * e2x)) / (sx * e2y - sy * e2x);
	if (abs(t2) < Mnull)
	{
		t2 = 0.0;
		if (!(r4 == r3) && !(r4 == r2))
		{
			result = NULL;
			return result;
		}
	}
	t1 = double(((r1->X - r4->X) * e1y - (r1->Y - r4->Y) * e1x)) / (sx * e1y - sy * e1x);
	if (abs(t1) < Mnull)
	{
		t1 = 0.0;
		if (!(r4 == r1) && !(r4 == r0))
		{
			result = NULL;
			return result;
		}
	}
	if (t1 * t2 > 0.0)
	{
		result = NULL;
		return result;
	}
	z1x = (r4->X) + sx * t1;
	z1y = (r4->Y) + sy * t1;
	z2x = (r4->X) + sx * t2;
	z2y = (r4->Y) + sy * t2;
	if (t1 * t2 == 0.0)
	{
		if (t1 == 0.0)
			l = 0.5 * (t1 + t2) * (sx * e1x + sy * e1y) / (ex * e1x + ey * e1y);
		else
			l = 0.5 * (t1 + t2) * (sx * e2x + sy * e2y) / (ex * e2x + ey * e2y);
		XC = 0.5 * (z1x + z2x) - ex * l;
		YC = 0.5 * (z1y + z2y) - ey * l;
		r1x = XC - z1x;
		r1y = YC - z1y;
		if (r1x * e1y - r1y * e1x > 0.0)
		{
			result = NULL;
			return result;
		}
		r2x = XC - z2x;
		r2y = YC - z2y;
		if (r2x * e2y - r2y * e2x < 0.0)
		{
			result = NULL;
			return result;
		}
		Rad = sqrt(r1x * r1x + r1y * r1y);
		if ((t2 == 0.0) && !Tangent(r0, r1, XC, YC))
		{
			result = NULL;
			return result;
		}
		if ((t1 == 0.0) && !Tangent(r2, r3, XC, YC))
		{
			result = NULL;
			return result;
		}
	}
	else
	{
		l = sqrt(abs(t1 * t2));
		r1x = z1x + l * e1x;
		r1y = z1y + l * e1y;
		if (!Tangent(r0, r1, r1x, r1y))
		{
			result = NULL;
			return result;
		}
		r2x = z2x + l * e2x;
		r2y = z2y + l * e2y;
		if (!Tangent(r2, r3, r2x, r2y))
		{
			result = NULL;
			return result;
		}
		Circle((r4->X), (r4->Y), r1x, r1y, r2x, r2y, &XC, &YC, &Rad);
	}
	if (Rad > Mnull)
		result = NewDisc(XC, YC, Rad);
	else
		result = NULL;
	return result;
} /*S2_P1*/


TDisc* S1_P2(Point *r0, Point *r1, Point *r2, Point *r3)
{
	TDisc* result;
	double ux = 0.0, uy = 0.0, vx = 0.0, vy = 0.0, wx = 0.0, wy = 0.0, al = 0.0, uv = 0.0, uw = 0.0, vw = 0.0, RRx = 0.0, RRy = 0.0, r1x = 0.0, r1y = 0.0, D = 0.0, t = 0.0, zx = 0.0, zy = 0.0, rx = 0.0, ry = 0.0, XC = 0.0, YC = 0.0, Rad = 0.0;
	ux = (r1->X - r0->X);
	uy = (r1->Y - r0->Y);
	vx = (r2->X - r0->X);
	vy = (r2->Y - r0->Y);
	wx = (r3->X - r2->X);
	wy = (r3->Y - r2->Y);
	uv = ux * vy - uy * vx;
	uw = ux * wy - uy * wx;
	if (uv * (uv + uw) < 0.0)   /* r2 и r3 лежат по разные стороны от r0,r1 */
	{
		result = NULL;
		return result;
	}
	if (abs(uw) > Mnull)
	{                   /* r2,r3 не перпендикулярен r0,r1 */
		vw = (vx * wy - vy * wx);
		al = vw / uw;
		rx = (r0->X) + ux * al;
		ry = (r0->Y) + uy * al;
		RRx = (rx - (r2->X));
		RRy = (ry - (r2->Y));
		r1x = (rx - (r3->X));
		r1y = (ry - (r3->Y));
		t = sqrt((RRx * RRx + RRy * RRy) * (r1x * r1x + r1y * r1y));
		t = sqrt(t);
		D = sqrt(ux * ux + uy * uy);
		if (uw > 0.0)
		{
			zx = rx - ux * t / D;
			zy = ry - uy * t / D;
		}
		else
		{
			zx = rx + ux * t / D;
			zy = ry + uy * t / D;
		}
	}
	else
	{                   /* r2,r3 перпендикулярен r0,r1 */
		t = -(vx * uy - vy * ux) / (wy * uy + wx * ux);
		zx = (r2->X) + wx / 2.0 + wy * t;
		zy = (r2->Y) + wy / 2.0 - wx * t;
	}
	if (Tangent(r0, r1, zx, zy))
	{
		Circle((r2->X), (r2->Y), (r3->X), (r3->Y), zx, zy, &XC, &YC, &Rad);
		result = NewDisc(XC, YC, Rad);
	}
	else
		result = NULL;
	return result;
} /*S1_P2*/


bool InterPoint(Point *r, TDisc *C)
{
	bool result = false;
	double u = 0.0, V = 0.0, r1 = 0.0;
	u = r->X - C->X;
	V = r->Y - C->Y;
	r1 = u * u + V * V;
	if (C->Rad == 0.0)
		result = (r1 < Mnull);
	else
		result = (r1 <= C->Rad * C->Rad - Mnull);
	return result;
} /*InterPoint*/


bool InterEdge(Point *r1, Point *r2, TDisc *C)
{
	bool result = false;
	double u1 = 0.0, V1 = 0.0, u2 = 0.0, V2 = 0.0, u3 = 0.0, v3 = 0.0, r = 0.0, S = 0.0;
	u1 = (r1->X) - C->X;
	V1 = (r1->Y) - C->Y;
	u2 = (r2->X) - C->X;
	V2 = (r2->Y) - C->Y;
	if (Tangent(r1, r2, C->X, C->Y))
	{
		r = u1 * V2 - u2 * V1;
		r = r * r;
		u3 = u1 - u2;
		v3 = V1 - V2;
		S = (u3 * u3 + v3 * v3);
		if (S > Mnull)
			r = r / S;
	}
	else
	{
		r = (u1 * u1 + V1 * V1);
		u3 = (u2 * u2 + V2 * V2);
		if (u3 < r)
			r = u3;
	}
	if (C->Rad == 0.0)
		result = (r < Mnull);
	else
		result = (r <= C->Rad * C->Rad - Mnull);
	return result;
} /*InterEdge*/


double DistPoint(Point *r, Point *rp)
{
	double result = 0.0;
	double u = 0.0, V = 0.0;
	u = (r->X - rp->X);
	V = (r->Y - rp->Y);
	result = sqrt(u * u + V * V);
	return result;
} /*DistPoint*/


double DistEdge(Point *r1, Point *r2, Point *rp)
{
	double result = 0.0;
	double u1 = 0.0, V1 = 0.0, u2 = 0.0, V2 = 0.0, u3 = 0.0, v3 = 0.0, r = 0.0, X = 0.0, Y = 0.0;
	u1 = (r1->X - rp->X);
	V1 = (r1->Y - rp->Y);
	u2 = (r2->X - rp->X);
	V2 = (r2->Y - rp->Y);
	X = (rp->X);
	Y = (rp->Y);
	if (Tangent(r1, r2, X, Y))
	{
		r = u1 * V2 - u2 * V1;
		r = r * r;
		u3 = u1 - u2;
		v3 = V1 - V2;
		r = r / (u3 * u3 + v3 * v3);
		r = sqrt(r);
	}
	else
		r = MAXDOUBLE;
	result = r;
	return result;
} /*DistEdge*/


bool SameCircle(TDisc *CircP, TDisc *CircQ)
{
	bool result = false;
	double E = 0.0;
	E = 0.01 * (CircP->Rad + CircQ->Rad);
	if (E > 0.1)
		E = 0.1;
	result = ((abs(CircP->X - CircQ->X) < E) && (abs(CircP->Y - CircQ->Y) < E) && (abs(CircP->Rad - CircQ->Rad) < E));
	return result;
} /*SameCircle*/


bool InterTriangle(Point *r0, Point *r1, Point *r2, Point *r)
{
	bool result = false;
	double S1 = 0.0, S2 = 0.0, S3 = 0.0;  /* Попадание точки r внутрь треугольника r0,r1,r2 */
	S1 = (r0->X - r->X) * (r1->Y - r->Y) - (r0->Y - r->Y) * (r1->X - r->X);
	S2 = (r1->X - r->X) * (r2->Y - r->Y) - (r1->Y - r->Y) * (r2->X - r->X);
	S3 = (r2->X - r->X) * (r0->Y - r->Y) - (r2->Y - r->Y) * (r0->X - r->X);
    result = ((S1 > 0) && (S2 > 0) && (S3 > 0)) || ((S1 < 0) && (S2 < 0) && (S3 < 0));
	return result;
}


bool InterPointEdge(Point *r, Point *r0, Point *r1, Point *r2, TDisc *Circ1, TDisc *Circ2)
{
	bool result = false;
	Point *P1, *P2;
	double X1 = 0.0, Y1 = 0.0, X2 = 0.0, Y2 = 0.0, m1 = 0.0, a = 0.0; /* Попадание точки r внутрь треугольника с вершинами: r0 и точки
		   касания сегмента r1,r2 и кругов Circ1,Circ2 */
	X1 = r2->X - r1->X;
	Y1 = r2->Y - r1->Y;
	m1 = X1 * X1 + Y1 * Y1;  /*квадрат длины сегмента r1,r2 */
	X2 = Circ1->X - r1->X;
	Y2 = Circ1->Y - r1->Y;
	a = (X1 * X2 + Y1 * Y2) / m1;
	P1 = new Point(r1->X + X1 * a, r1->Y + Y1 * a); /* касание r1,r2 и Circ1 */
	X2 = Circ2->X - r1->X;
	Y2 = Circ2->Y - r1->Y;
	a = (X1 * X2 + Y1 * Y2) / m1;
	P2 = new Point(r1->X + X1 * a, r1->Y + Y1 * a);  /* касание r1,r2 и Circ2 */
	result = InterTriangle(r0, P1, P2, r);
	delete P1;
	delete P2;
	return result;
}



bool InterTwoEdges(Point *r, Point *r1, Point *r2, Point *r3, Point *r4, TDisc *Circ1, TDisc *Circ2)
{
	bool result = false;
	double X1 = 0.0, X2 = 0.0, Y1 = 0.0, Y2 = 0.0, a = 0.0, B = 0.0, C = 0.0, m1 = 0.0;
	TDisc *C1, *C2;  /* Попадание точки r внутрь трапеции с вершинами в точках
		 касания кругов Circ1,Circ2 с сегментами r1,r2 и r3,r4 */
	result = false;
	/*Последовательность кругов*/
	X1 = r2->X - r1->X;
	Y1 = r2->Y - r1->Y;
	X2 = Circ2->X - Circ1->X;
	Y2 = Circ2->Y - Circ1->Y;
	a = X1 * X2 + Y1 * Y2;
	if (a >= 0)
	{
		C1 = Circ1;
		C2 = Circ2;
	}
	else
	{
		C1 = Circ2;
		C2 = Circ1;
	}
	/*Положение точки относительно линии центров*/
	X1 = C2->X - C1->X;
	Y1 = C2->Y - C1->Y;
	X2 = r->X - C1->X;
	Y2 = r->Y - C1->Y;
	if (X1 * Y2 - X2 * Y1 <= 0)
	{
		/* первый сегмент r1,r2 */
		X1 = r2->X - r1->X;
		Y1 = r2->Y - r1->Y;
		m1 = sqrt(X1 * X1 + Y1 * Y1);  /* длина сегмента r1,r2 */
		X2 = r->X - r1->X;
		Y2 = r->Y - r1->Y;
		if (X1 * Y2 - X2 * Y1 <= 0)
			return result; /* точка r справа от сегмента r1,r2 */
		B = (X1 * X2 + Y1 * Y2) / m1;    /* проекция точки r на r1,r2 */
		X2 = C1->X - r1->X;
		Y2 = C1->Y - r1->Y;
		a = (X1 * X2 + Y1 * Y2) / m1;    /* проекция центра C1 на r1,r2 */
		X2 = C2->X - r1->X;
		Y2 = C2->Y - r1->Y;
		C = (X1 * X2 + Y1 * Y2) / m1;    /* проекция центра Circ2 на r1,r2 */
		if ((a > B) || (C < B))
			return result;
	}
	else
	{
		/* второй сегмент r3,r4 */
		X1 = r4->X - r3->X;
		Y1 = r4->Y - r3->Y;
		m1 = sqrt(X1 * X1 + Y1 * Y1);  /* длина сегмента r3,r4 */
		X2 = r->X - r3->X;
		Y2 = r->Y - r3->Y;
		if (X1 * Y2 - X2 * Y1 <= 0)
			return result; /* точка r справа от сегмента r3,r4 */
		B = (X1 * X2 + Y1 * Y2) / m1;    /* проекция точки r на r3,r4 */
		X2 = C1->X - r3->X;
		Y2 = C1->Y - r3->Y;
		a = (X1 * X2 + Y1 * Y2) / m1;    /* проекция центра Circ1 на r3,r4 */
		X2 = C2->X - r3->X;
		Y2 = C2->Y - r3->Y;
		C = (X1 * X2 + Y1 * Y2) / m1;    /* проекция центра Circ2 на r3,r4 */
		if ((a < B) || (C > B))
			return result;
	}
	result = true;
	return result;
} /*Geometry*/
