#pragma once

/* Модуль содержит процедуры непосредственных геометрических рассчетов */

#include "StructureTD.h"
#include <cmath>
using namespace std;

#define Sqr(x) ((x)*(x))

TDisc* P3(Point *r0, Point *r1, Point *r2);
/* Построение круга, касающегося трех точек */
TDisc* S1_P1(Point *r0, Point *r1, Point *r2, Point *r3);
/* Построение круга, касающегося сегмента в его конце и еще одной точки */
TDisc* S3(Point *r0, Point *r1, Point *r2, Point *r3, Point *r4, Point *r5);
/* Построение круга, касающегося трех сегментов */
TDisc* S2_P1(Point *r0, Point *r1, Point *r2, Point *r3, Point *r4);
/* Построение круга, касающегося двух сегментов и точки */
TDisc* S1_P2(Point *r0, Point *r1, Point *r2, Point *r3);
/* Построение круга, касающегося сегмента и двух точек */
bool InterPoint(Point *r, TDisc *C);
/* Установление инцидентности точки и окружности */
bool InterEdge(Point *r1, Point *r2, TDisc *C);
/* Установление инцидентности сегмента и окружности */
TDisc* NewDisc(double X, double Y, double Rad);
/* Порождение нового круга */
double DistPoint(Point *r, Point *rp);
/* Расстояние между двумя точками */
double DistEdge(Point *r1, Point *r2, Point *rp);
/* Расстояние между сегментом и точкой */
bool SameCircle(TDisc *CircP, TDisc *CircQ);
/* Установление совпадения окружностей */
bool InterTriangle(Point *r0, Point *r1, Point *r2, Point *r);
/* Попадание точки r внутрь треугольника r0,r1,r2 */
bool InterPointEdge(Point *r, Point *r0, Point *r1, Point *r2, TDisc *Circ1, TDisc *Circ2);
/* Попадание точки r внутрь треугольника с вершинами: r0 и точки
   касания сегмента r1,r2 и кругов Circ1,Circ2 */
bool InterTwoEdges(Point *r, Point *r1, Point *r2, Point *r3, Point *r4, TDisc *Circ1, TDisc *Circ2);
/* Попадание точки r внутрь трапеции с вершинами в точках
   касания кругов Circ1,Circ2 с сегментами r1,r2 и r3,r4 */
bool Collinear(Point *r0, Point *r1, Point *r2);
/*коллинеарность трех собственных точек*/
bool Codirect(Point *r0, Point *r1, Point *r2);
/*монотонность трех собственных точек*/
PlanePosition Classify(Point *P1, Point *P2, Point *q);
/*Положение точки q относительно вектора p1,p2 */
