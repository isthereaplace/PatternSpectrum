#pragma once

/*Модуль содержит процедуры построения графа смежности сайтов
  многоугольной фигуры в виде обобщённой триангуляции Делоне,
  состоящей из вершин, рёбер и граней.*/


#include "LinkedList.h"
#include "StructureTD.h"
#include "Geometry.h"
#include <QDebug>
#include <vector>
using std::vector;

void ProduceElements(Domain* Collect);
/* Построение элементов триангуляции */
void CreateTriangulation(Domain* Collect);
/* Построение трангуляции и скелета */
