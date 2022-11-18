#include "BSTrans.h"

void BondSkeletTrans(BitRaster* InputImg  /*указатель на исходный образ*/
                    ,int PruningSize  /*размер стрижки скелета*/
                    ,int AreaIgnore  /*площадь игнорируемых контуров*/
                    ,TPolFigure*& Figure  /*гранично-скелетное представление фигуры*/)
{
    Figure = new TPolFigure(InputImg, AreaIgnore);
    Figure->MakeTriangDel();
    Figure->CutSkeleton(PruningSize);
}
