#pragma once

#include "StructureSkel.h"
#include "SpanTree.h"
#include "BitRaster.h"
#include "StructureTD.h"
#include "TriDel.h"
#include "LinkedList.h"

/*
Подпрограмма гранично-скелетного преобразования осуществляет построение
граничных контуров и скелета. Контура, ограничивающие площадь менее
AreaIgnore, уничтожаются (игнорируются), полученый скелет регуляризируется
с помощью операции стрижки с параметром PruningSize
*/
void BondSkeletTrans(BitRaster* InputImg  /*указатель на исходный образ*/
                    ,int PruningSize  /*размер стрижки скелета*/
                    ,int AreaIgnore  /*площадь игнорируемых контуров*/
                    ,TPolFigure*& Figure  /*гранично-скелетное представление фигуры*/);
