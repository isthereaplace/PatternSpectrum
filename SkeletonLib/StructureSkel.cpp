#include "StructureSkel.h"
#include "delaunator.hpp"
using namespace delaunator;

#define _USE_MATH_DEFINES
#include <math.h>

typedef unsigned int unsignedint;


int NodeCount = 0, BoneCount = 0, CompCount = 0;


unsignedint TimeInMilSecond()
/*Текущее время в миллисекундах*/
{
	return 1000 * (clock() / float(CLOCKS_PER_SEC));
}


double TNode::X()
{
	return Disc->X;	
}


double TNode::Y()
{
	return Disc->Y;
}


double TNode::r()
{
	return Disc->Rad;
}


NodeKind TNode::Kind()
{
	if (Bones[1] == NULL)
		return Out1_In0;
	if (Bones[2] == NULL)
		if (Bones[0]->dest == this)		// Есть входящее ребро
			return Out0_In2;
		else if (Bones[1]->dest == this)
			return Out1_In1;
		else
			return Out2_In0;
	else if (Bones[0]->org == this)		// Есть исходящее ребро
		return Out1_In2;
	else
		return Out0_In3;
}


TConnected::TConnected()
	: Nodes(new LinkedListTail<TNode>),
	Bones(new LinkedListTail<TBone>)
{
	CompCount++;
}


TNode::TNode()
	: Number(0),
	Disc(NULL),
	Depth(0.0),
	Merge(NULL)
{
	for (int i = 1; i <= 3; i++) {
		Bones[i - 1] = NULL;
		Sites[i - 1] = NULL;
	}
	NodeCount++;
}


TNode::~TNode()
{
	NodeCount--;
	if (Disc != NULL)
		delete Disc;
	// todo check:  inherited::Destroy;
}


void TNode::DetachBone(TBone* Bone) /*отцепить кость*/
{
	if (Bones[1 - 1] == Bone)
	{
		Bones[1 - 1] = Bones[2 - 1];
		Bones[2 - 1] = Bones[3 - 1];
		Bones[3 - 1] = NULL;
	}
	else
		if (Bones[2 - 1] == Bone)
		{
		Bones[2 - 1] = Bones[3 - 1];
		Bones[3 - 1] = NULL;
		}
		else
			if (Bones[3 - 1] == Bone)
				Bones[3 - 1] = NULL;
}


TBone::TBone()
	: Com(NULL),
	org(NULL),
	dest(NULL),
	Met(false),
	Virt(NULL),
	type(Linear),
	t(0),
	s(0),
	p(0),
	x0(0),
	y0(0),
	dx(1),
	dy(0),
	dir(1),
	Lacuna(NULL)
{
	BoneCount++;
	Sites[0] = Sites[1];
}


TNode* TBone::GetNextNode(TNode* Source)
/*узел, противоположный Node на этой же кости*/
{
	TNode* result = NULL;
	if (org == Source)
		result = dest;
	else
		if (dest == Source)
			result = org;
		else
			result = NULL;
	return result;
}


void TBone::DetachNode(TNode* Node)
{     /*отцепиться от узла*/
	if (Node == org)
		org = NULL;
	else
		if (Node == dest)
			dest = NULL;
}


TBone::~TBone()
{
	BoneCount--;
	if (Virt != NULL) {
		delete Virt;
		Virt = NULL;
	}
}


void TBone::DestroyWithDetach()
/*отцепиться и уничтожить кость*/
{
	if (org != NULL)
		org->DetachBone(this);
	if (dest != NULL)
		dest->DetachBone(this);
	BoneCount--;
	if (Virt != NULL) {
		delete Virt;
		Virt = NULL;
	}
	// todo check:  inherited::Destroy;
	// TODO find usage
}


bool TBone::Fiction()
/*кость вдоль фиктивного контура */
{
	bool result = false;
	bool b1 = false, b2 = false;
    b1 = false;
	b2 = false;
	for (int i = 1; i <= 3; i++)
	{
		b1 = b1 || org->Sites[i - 1]->Cont->Fiction;
		b2 = b2 || dest->Sites[i - 1]->Cont->Fiction;
	}
	result = b1 && b2;
	return result;
}


void TBone::Bisector(TSite *E1, TSite *E2)
/*Вычисление параболического бисектора*/
{
	double xv = 0.0, yv = 0.0, Z = 0.0;
	TPoint *s0 = NULL, *S1 = NULL, *S2 = NULL, *St = NULL;
	TDisc *C1 = NULL, *C2 = NULL, *ct = NULL;
	double ex = 0.0, ey = 0.0, fx = 0.0, fy = 0.0;
	double d0 = 0.0, d1 = 0.0, d2 = 0.0, t0 = 0.0, t1 = 0.0, t2 = 0.0, a = 0.0, B = 0.0, M = 0.0, XT = 0.0, yt = 0.0;
	/*#
	St1:
	St2:
	St3:
	St4:
	St5:
	St6:
	St7:
	St8:
	*/
	/*определение и правильная ориентация всех точек*/
	if (E1->IsVertex())
	{
		s0 = ((Vertex*)E1)->p;
		S1 = ((TEdge*)E2)->org;
		S2 = ((TEdge*)E2)->dest;
	}
	else
	{
		s0 = ((Vertex*)E2)->p;
		S1 = ((TEdge*)E1)->org;
		S2 = ((TEdge*)E1)->dest;
	}
	Z = (S1->X - s0->X) * (S2->Y - s0->Y) - (S2->X - s0->X) * (S1->Y - s0->Y);
	if (Z == 0)
		return;
	if (Z < 0)
	{
		St = S1;
		S1 = S2;
		S2 = St;
	}
	C1 = org->Disc;
	C2 = dest->Disc;
	Z = (S2->X - S1->X) * (C2->X - C1->X) + (S2->Y - S1->Y) * (C2->Y - C1->Y);
	if (Z < 0)
	{
		ct = C1;
		C1 = C2;
		C2 = ct;
	}
	/*вычисление контрольной точки*/
//St1:
	ex = S2->X - S1->X;
	ey = S2->Y - S1->Y;
	fx = -ey;
	fy = ex;
//St2:
	M = Sqr(ex) + Sqr(ey);
	d1 = (ex * (C1->Y - S1->Y) - ey * (C1->X - S1->X)) / M;
	d2 = (ex * (C2->Y - S1->Y) - ey * (C2->X - S1->X)) / M;
	d0 = (ex * (s0->Y - S1->Y) - ey * (s0->X - S1->X)) / M;
//St3:
	t1 = (ex * (C1->X - S1->X) + ey * (C1->Y - S1->Y)) / M;
	t2 = (ex * (C2->X - S1->X) + ey * (C2->Y - S1->Y)) / M;
	t0 = (ex * (s0->X - S1->X) + ey * (s0->Y - S1->Y)) / M;
	if (abs(t1 - t2) < 0.0001)
		return;
//St4:
	B = d0 / 2;
//St5:
	if (abs(t1 - t0) > abs(t2 - t0))
		a = (d1 - B) / Sqr(t1 - t0);
	else
		a = (d2 - B) / Sqr(t2 - t0);
//St6:
	XT = ((t2 - t0) * t2 - (t1 - t0) * t1 - (d2 - d1) / (2 * a)) / (t2 - t1);
	yt = d1 + (XT - t1) * 2 * a * (t1 - t0);
//St7:
	xv = S1->X + XT * ex + yt * fx;
	yv = S1->Y + XT * ey + yt * fy;
	Virt = new TPoint(xv, yv);
}


void TBone::BezierPoints(double& X1, double& Y1, double& X2, double& Y2)
/*Вычисление контрольных точек для параболического ребра*/
{
	if (Virt == NULL)
	{
		X1 = org->X();
		Y1 = org->Y();
		X2 = dest->X();
		Y2 = dest->Y();
	}
	else
	{
		X1 = org->X() + double(2) * (Virt->X - org->X()) / 3;
		Y1 = org->Y() + double(2) * (Virt->Y - org->Y()) / 3;
		X2 = dest->X() + double(2) * (Virt->X - dest->X()) / 3;
		Y2 = dest->Y() + double(2) * (Virt->Y - dest->Y()) / 3;
		//       x1:=Virt.X; y1:=Virt.Y;
		//       x2:=Virt.X; y2:=Virt.Y;
	}
}


TConnected::~TConnected()
/*Уничтожение компоненты*/
{
	TBone* B = NULL;
	TNode* n = NULL;
	HoleList.clear();
	while (!Bones->isEmpty())
	{
		B = (TBone*)Bones->first();
		B->removeFromCurrentList();
		delete B;
	}
	while (!Nodes->isEmpty())
	{
		n = (TNode*)Nodes->first();
		n->removeFromCurrentList();
		delete n;
	}
	delete Bones;
	delete Nodes;
	CompCount--;
	// todo check:  inherited::Destroy;
}


void TConnected::CutSkeleton(double Eps)
/*стрижка-укорачивание на Eps - точность*/
{
	LinkedListTail<TNode>* Hairs;
	TNode* Node = NULL, *Node1 = NULL, *t = NULL;
	TBone* Bone = NULL;
	double l = 0.0, D = 0.0;
	/*поиск всех терминальных узлов*/
	Hairs = new LinkedListTail < TNode > ;
	Node = (TNode*)Nodes->first();
	while (Node != NULL)
	{
		Node1 = Node->getNext();
		Node->Depth = Node->r();
		if (Node->Kind() == Out1_In0)
			Node->moveIntoTail(Hairs);
		Node = Node1;
	}
	/*собственно стрижка*/
	while (!Hairs->isEmpty())
	{
		t = (TNode*)Hairs->first();
		Bone = t->Bones[1 - 1];
		if (Bone != NULL)
		{
			Node = Bone->GetNextNode(t);
			l = sqrt(Sqr(t->X() - Node->X()) + Sqr(t->Y() - Node->Y()));
			if ((l + t->Depth - Node->r() < Eps) && (Node->Kind() != Out1_In0))
			{
				D = l + t->Depth;
				if (D > Node->Depth)
					Node->Depth = D;
				t->DetachBone(Bone);
				Bone->DetachNode(t);
				delete t;
				Bone->DetachNode(Node);
				Node->DetachBone(Bone);
				Bone->DestroyWithDetach(); /// !!!
				delete Bone;
				if (Node->Kind() == Out1_In0)
					Node->moveIntoTail(Hairs);
			}
			else
				t->moveIntoTail(Nodes);
		}
		else
			delete t;
	}
	delete Hairs;
}


void TPolFigure::ClearAll()
{
	TConnected* Com = NULL; /*уничтожение компонент*/
	while (!Components->isEmpty())
	{
		Com = (TConnected*)Components->first();
		delete Com;
	}
}


TPolFigure::~TPolFigure()
{
	ClearAll();
	delete Components;
	Components = NULL;
	// todo check:  inherited::Destroy;
}


TPolFigure::TPolFigure(BitRaster* bitRaster, double Amin)
/*Построение компонент для матрицы с отбрасыванием малых контуров*/
{
	unsigned int TT, TT1;
	TT1 = TimeInMilSecond();
	Components = new LinkedListTail < TConnected >;		/*Порождение списка компонент*/
	
	ContourTracer* BinIm = new ContourTracer(bitRaster, Amin);
	BinIm->traceContours();

	CreateContours(BinIm);
	
	/*Упорядочение контуров АВЛ деревом*/
	if (!ElementsExist)
	{
		ProduceElements(this);
		ElementsExist = true;
	}
	TT = TimeInMilSecond();
	SpaningTree(this);
	MakeComponents();

	RTab.TimeTree = TimeInMilSecond() - TT;
	RTab.ConnectComp = CompCount;
	delete BinIm;
	RTab.TimeTrace = TimeInMilSecond() - TT1 - RTab.TimeTree;
	RTab.Total = RTab.TimeTrace + RTab.TimeTree + RTab.TimeSkelet + RTab.TimePrun;
}


TPolFigure::TPolFigure(QString path)
/*Построение компонент, записанных в текстовом файле*/
{
	unsigned int TT, TT1;
	TT1 = TimeInMilSecond();
	Components = new LinkedListTail < TConnected >;		/*Порождение списка компонент*/

	TContour* S;       /* контур для скелетизации (модуль Structure) */
	int i = 0;
	RTab.Points = 0;
	RTab.Polygons = 0;

	ifstream src;
	src.open(path.toStdString());
	int nComp, nCont, nPoint;
	double x, y;
	src >> nComp;
	for (int iComp = 0; iComp < nComp; iComp++) {
		src >> nCont;
		for (int iCont = 0; iCont < nCont; iCont++) {
			S = AddContour();
			S->Internal = (iCont > 0);
			src >> nPoint;
			for (int iPoint = 0; iPoint < nPoint; iPoint++) {
				src >> x >> y;
				S->AddPoint(x, -y);
			}
			if (nPoint < 3)
				delete S;
			else
			{
				S->ShiftHead();
				if (S->Internal == S->ConterClockWise())
					S->Invert();
				S->ClosedPath = true;
				RTab.Polygons++;
				RTab.Points = RTab.Points + S->ListPoints->cardinal();
			}
		}
	}
	src.close();

	/*Упорядочение контуров АВЛ деревом*/
	if (!ElementsExist)
	{
		ProduceElements(this);
		ElementsExist = true;
	}
	TT = TimeInMilSecond();
	SpaningTree(this);
	MakeComponents();

	RTab.TimeTree = TimeInMilSecond() - TT;
	RTab.ConnectComp = CompCount;
	RTab.TimeTrace = TimeInMilSecond() - TT1 - RTab.TimeTree;
	RTab.Total = RTab.TimeTrace + RTab.TimeTree + RTab.TimeSkelet + RTab.TimePrun;
}


void TPolFigure::CreateContours(ContourTracer *BinIm)
/*Построение контуров для бинарного образа BinIm*/
{
    ClosedPath* CP; /* контур из растрового образа (модуль Kontur) */
	RasterPoint* r;  /* точка контура (модуль Kontur) */
	TContour* S;       /* контур для скелетизации (модуль Structure) */
	int i = 0;
	RTab.Points = 0;
	RTab.Polygons = 0;
	CP = BinIm->initialContour();
	while (CP != NULL)
	{
		S = AddContour();
		S->Internal = CP->Internal;
		r = CP->initialPoint();
		i = 0;
		while (r != NULL)
		{
			S->AddPoint(r->x, r->y);
			i++;
			r = r->getNext();
		}
		if (i < 3)
			delete S;
		else
		{
			S->ShiftHead();
			if (S->Internal == S->ConterClockWise())
				S->Invert();
			S->ClosedPath = true;
			RTab.Polygons++;
			RTab.Points = RTab.Points + S->ListPoints->cardinal();
		}
		CP = CP->getNext();
	}
}

void TPolFigure::MakeComponents()
/*Формирование компонент из контуров*/
{
	TContour *S = NULL, *S1 = NULL, *S2 = NULL;
	TConnected* Com = NULL;
	/*Временное размещение внутренних контуров во внешние*/
	S = (TContour*)Boundary->first();
	while (S != NULL)
	{
		S1 = S->getNext();
		if (S->Internal)
		{
			S2 = S->Container;
			if (S2->MySons == NULL)
				S2->MySons = new LinkedListTail < TContour > ;
			S->moveIntoTail(S2->MySons);
		}
		S = S1;
	}
	/*Из каждого внешнего создается компонента*/
	S = (TContour*)Boundary->first();
	while (S != NULL)
	{
		if (!S->Internal)
		{
			Com = new TConnected;
			Com->Border = S;
			Com->moveIntoTail(Components);
			while (S->MySons != NULL)
			{
				while (!S->MySons->isEmpty())
				{
					S1 = (TContour*)S->MySons->first();
					Com->HoleList.push_back(S1);
					S1->moveIntoTail(Boundary);
				}
				delete S->MySons;
				S->MySons = NULL;
			}
		}
		S = S->getNext();
	}
}


void TPolFigure::Invert()/*Инверсия фигуры*/
{
	TContour* S = NULL;
	Point *p = NULL;
	double xmin = 0.0, Xmax = 0.0, ymin = 0.0, Ymax = 0.0;
	int W = 0;
	W = 100;
	ClearAll();
	xmin = 10000;
	Xmax = -10000;
	ymin = 10000;
	Ymax = -10000;
	S = (TContour*)Boundary->first();
	while (S != NULL)
	{
		S->Internal = !S->Internal;
		p = (Point*)S->ListPoints->first();
		while (p != NULL)
		{
			if (p->X < xmin)
				xmin = p->X;
			if (p->X > Xmax)
				Xmax = p->X;
			if (p->Y < ymin)
				ymin = p->Y;
			if (p->Y > Ymax)
				Ymax = p->Y;
			p = p->getNext();
		}
		S->Invert();
		S->Container = NULL;
		S->ClosestSite = NULL;
		if (S->MySons != NULL)
			delete S->MySons;
		S->MySons = NULL;
		S = S->getNext();
	}
	S = new TContour;
	S->Internal = false;
	p = new Point(xmin - W, ymin - W);
	p->moveIntoTail(S->ListPoints);
	p = new Point(Xmax + W, ymin - W);
	p->moveIntoTail(S->ListPoints);
	p = new Point(Xmax + W, Ymax + W);
	p->moveIntoTail(S->ListPoints);
	p = new Point(xmin - W, Ymax + W);
	p->moveIntoTail(S->ListPoints);
	S->moveIntoTail(Boundary);
	S->Fiction = true;
	S->ShiftHead();
	if (S->Internal == S->ConterClockWise())
		S->Invert();
	ElementsExist = false;
	ProduceElements(this);
	ElementsExist = true;
}


void TPolFigure::RestoreInversion()
/*Восстановление фигуры после инверсии*/
{
	TContour *S = NULL, *s0 = NULL;
    /*Удаление частей скелета, инцидентных фиктивному контуру*/
	ClearAll();
	SkelExist = false;
	/*  Com:=Components.first AS TConnected;
	  WHILE Com<>NIL DO
	  BEGIN
	  Bone:=Com.Bones.first AS TBone;
	  WHILE Bone<>NIL DO
	  BEGIN
	  Node:=Bone.Org;
	  Node1:=Bone.Dest;
	  Bone1:=Bone.suc AS TBone;
	  IF Bone.Fiction THEN (*ребро инцидентное фиктивному контуру*)
	  BEGIN
	  Bone.DestroyWithDetach;
	  IF Node.Kind=Isolated THEN Node.Destroy;
	  IF Node1.Kind=Isolated THEN Node1.Destroy;
	  END;
	  Bone:=Bone1;
	  END;
	  Com:=Com.suc AS TConnected;
	  END; */
	/*Инвертирование направления контуров*/
	S = (TContour*)Boundary->first();
	while (S != NULL)
	{
		S->Internal = !S->Internal;
		if (S->Fiction)
			s0 = S;
		else
		{
			S->Invert();
			S->Container = NULL;
			S->ClosestSite = NULL;
		}
		S = S->getNext();
	}
	delete s0;
}

void printNode(TNode *n)
{
    qDebug() << n->Disc->X << " " << n->Disc->Y << " " << n->Disc->Rad << endl;
}

void TPolFigure::MakeNodeBoneRepresentation()
/* Формирование структуры типа TPolFigure, описывающей
границу и внутренний скелет бинарной области из триангуляции.
Сама триангуляция при этом разрушается*/
{
	vector<Triplet*> ListOld;
	vector<TNode*> ListNew; /*Рабочие списки*/
    int M = 0;
	TNode* Node = NULL, *Node1 = NULL;
	TSite *S1 = NULL, *S2 = NULL;
	Triplet *Tr = NULL, *Tr1 = NULL;
	TBone* Bone = NULL;
	TConnected* Com = NULL;
	if ((Boundary == NULL) || Boundary->isEmpty())
		return;
	Com = (TConnected*)Components->first();
	while (Com != NULL)
	{
		Tr = (Triplet*)Com->Border->Map->MapTriplet->first();
		M = 0;
		while (Tr != NULL)
		{
			Node = new TNode;
			Node->Sites[1 - 1] = Tr->E1;
			Node->Sites[2 - 1] = Tr->E2;
			Node->Sites[3 - 1] = Tr->e3;
			Node->Disc = Tr->Circ;
			Tr->Circ = NULL;
			M++;
			Tr->Numb = M;
			ListOld.push_back(Tr);
			ListNew.push_back(Node);
			Node->moveIntoTail(Com->Nodes);  /* В список узлов */
			Tr = Tr->getNext();
            //printNode(Node);
		}
		/* M - число узлов*/

		/*Образование костей*/
		for (int stop = M - 1, i = 0; i <= stop; i++)
		{
			Tr = ListOld[i];
			Node = ListNew[i];
			for (int stop = 3, j = 1; j <= stop; j++)
			{
				switch (j)
				{
				case 1:
				{
					Tr1 = Tr->t1;
					S1 = Tr->e3;
					S2 = Tr->E1;
				}
					break;
				case 2:
				{
					Tr1 = Tr->t2;
					S1 = Tr->E1;
					S2 = Tr->E2;
				}
					break;
				case 3:
				{
					Tr1 = Tr->t3;
					S1 = Tr->E2;
					S2 = Tr->e3;
				}
					break;
				default:
					Tr1 = NULL;
				}
				if ((Tr1 != NULL) && (Tr->Numb < Tr1->Numb))
				{
					Node1 = ListNew[Tr1->Numb - 1];
					Bone = new TBone;
					Bone->org = Node;
					Bone->dest = Node1;
					if (S1->IsVertex() == !S2->IsVertex())
						Bone->Bisector(S1, S2); /*параболическое ребро*/
					Bone->moveIntoTail(Com->Bones);   /*В список костей*/
					Node->AddBone(Bone);
					Node1->AddBone(Bone);
				}
			}
		}
		ListOld.clear();
		ListNew.clear();
		delete Com->Border->Map;
		Com->Border->Map = NULL;
		Com = Com->getNext();
	}
	SkelExist = true;
	ListOld.clear();
	ListNew.clear();
	RTab.Vertex = NodeCount;
	RTab.Edges = BoneCount;
}


void TPolFigure::MakeTriangDel()
{
	unsigned int TT = 0;
	TT = TimeInMilSecond();
	if (!ElementsExist)
	{
		ProduceElements(this);
		ElementsExist = true;
	}
	CreateTriangulation(this);
	MakeNodeBoneRepresentation();
	MapExist = true;
	RTab.TimeSkelet = TimeInMilSecond() - TT;
	RTab.Total = RTab.TimeTrace + RTab.TimeTree + RTab.TimeSkelet + RTab.TimePrun;
}


void TPolFigure::CutSkeleton(double Eps)
/*Стрижка скелета*/
{
	TConnected* Com = NULL;
	unsigned int TT = 0;
	TT = TimeInMilSecond();
	Com = (TConnected*)Components->first();
	while (Com != NULL)
	{
		Com->CutSkeleton(Eps);
		Com = Com->getNext();
	}
	RTab.TimePrun = TimeInMilSecond() - TT;
	RTab.Vertex = NodeCount;
	RTab.Edges = BoneCount;
}


void TNode::AddBone(TBone* Bone) /*добавить кость*/
{
	if (Bones[1 - 1] == NULL)
		Bones[1 - 1] = Bone;
	else
		if (Bones[2 - 1] == NULL)
			Bones[2 - 1] = Bone;
		else
			if (Bones[3 - 1] == NULL)
				Bones[3 - 1] = Bone;
}


bool TNode::IntersectInside(TBone* BoneA, TBone* BoneB) {
	QPointF ExtrA = BoneA->GetExtremePoint(r());
	QPointF ExtrB = BoneB->GetExtremePoint(r());
	double d = sqrt(pow(ExtrA.x() - ExtrB.x(), 2) + pow(ExtrA.y() - ExtrB.y(), 2));
	if (d < 2*r()) {
		double alpha = atan2(ExtrB.y() - ExtrA.y(), ExtrB.x() - ExtrA.x());
		double beta = acos(d / (2 * r()));
		double x1 = ExtrA.x() + r()*cos(alpha - beta);
		double y1 = ExtrA.y() + r()*sin(alpha - beta);
		double x2 = ExtrA.x() + r()*cos(alpha + beta);
		double y2 = ExtrA.y() + r()*sin(alpha + beta);
		if (sqrt(pow(x1 - X(), 2) + pow(y1 - Y(), 2)) <= r() + 1e-6 && 
			sqrt(pow(x2 - X(), 2) + pow(y2 - Y(), 2)) <= r() + 1e-6)
			return true;
	}
	return false;
}


void TNode::RebuildNeighborhood() {
	list<TBone*> Arcs;
	switch (Kind())
	{
	case Out1_In0:
		return;
	case Out2_In0:
		Bones[0]->Adjacent.push_front(AdjElem(Bones[1], std::list<AdjElem>::iterator()));
		Bones[1]->Adjacent.push_front(AdjElem(Bones[0], std::list<AdjElem>::iterator()));
		Bones[0]->Adjacent.front().Link = Bones[1]->Adjacent.begin();
		Bones[1]->Adjacent.front().Link = Bones[0]->Adjacent.begin();
		return;
	case Out1_In1:
	case Out1_In2:
		for (int i = 1; i <= (Kind() == Out1_In1 ? 1 : 2); i++){
			for (auto Elem = Bones[i]->Adjacent.begin(); Elem != Bones[i]->Adjacent.end(); Elem++) {
				if (IntersectInside(Bones[0], Elem->Bone)) {
					Bones[0]->Adjacent.push_front(*Elem);
					Elem->Link->Bone = Bones[0];
					Elem->Link->Link = Bones[0]->Adjacent.begin();
				}
				else
					Elem->Bone->Adjacent.erase(Elem->Link);
			}
			Bones[i]->Adjacent.clear();
		}
		return;
	case Out0_In2:
	case Out0_In3:
		for (int i = 0; i <= (Kind() == Out0_In2 ? 1 : 2); i++) {
			for (auto Elem = Bones[i]->Adjacent.begin(); Elem != Bones[i]->Adjacent.end(); Elem++) {
				Arcs.push_front(Elem->Bone);
				Elem->Bone->Adjacent.erase(Elem->Link);
			}
			Bones[i]->Adjacent.clear();
		}
		for (auto iArc = Arcs.begin(); iArc != Arcs.end(); iArc++) {
			auto jArc = iArc; jArc++;
			while (jArc != Arcs.end()) {
				if (IntersectInside(*iArc, *jArc)) {
					(*iArc)->Adjacent.push_front(AdjElem(*jArc, std::list<AdjElem>::iterator()));
					(*jArc)->Adjacent.push_front(AdjElem(*iArc, std::list<AdjElem>::iterator()));
					(*iArc)->Adjacent.front().Link = (*jArc)->Adjacent.begin();
					(*jArc)->Adjacent.front().Link = (*iArc)->Adjacent.begin();
				}
				jArc++;
			}
		}
		return;
	}
}


Lacuna::Lacuna()
	: color(rand() % 256, rand() % 256, rand() % 256), Bones(0) {}


void Lacuna::Absorb(Lacuna* Inflow) {
	if (Inflow != NULL && Inflow != this) {
		if (this->Bones.size() > Inflow->Bones.size()) {
			for (auto iBone = Inflow->Bones.begin(); iBone != Inflow->Bones.end(); iBone++) {
				(*iBone)->Lacuna = this;
			}
			this->Bones.splice(this->Bones.end(), Inflow->Bones);
		}
		else {
			for (auto iBone = this->Bones.begin(); iBone != this->Bones.end(); iBone++) {
				(*iBone)->Lacuna = Inflow;
			}
			Inflow->Bones.splice(Inflow->Bones.end(), this->Bones);
		}
	}
}


Lacuna* TNode::UpdateLacunas() {
	int NumIn;
	Lacuna* L = NULL;
	switch (Kind())
	{
	case Out1_In0:
		L = new Lacuna();
		L->Bones.push_back(Bones[0]);
		Bones[0]->Lacuna = L;
		break;
	case Out2_In0:
		L = new Lacuna();
		L->Bones.push_back(Bones[0]);
		L->Bones.push_back(Bones[1]);
		Bones[0]->Lacuna = L;
		Bones[1]->Lacuna = L;
		break;
	case Out1_In1:
		if (Bones[1]->Lacuna != NULL) {
			Bones[1]->Lacuna->Bones.push_back(Bones[0]);
			Bones[0]->Lacuna = Bones[1]->Lacuna;
		}
		break;
	case Out1_In2:
		if (Bones[1]->Lacuna != NULL) {
			Bones[1]->Lacuna->Bones.push_back(Bones[0]);
			Bones[0]->Lacuna = Bones[1]->Lacuna;
			if (Bones[2]->Lacuna != NULL) {
				Bones[1]->Lacuna->Absorb(Bones[2]->Lacuna);
			}
		}
		else if (Bones[2]->Lacuna != NULL) {
			Bones[2]->Lacuna->Bones.push_back(Bones[0]);
			Bones[0]->Lacuna = Bones[2]->Lacuna;
		}
		break;
	case Out0_In2:
	case Out0_In3:
		NumIn = Kind() == Out0_In2 ? 2 :3;
		for (int i = 0; i < NumIn; i++) {
			for (int j = i + 1; j < NumIn; j++) {
				if (Bones[i]->Lacuna != NULL && Bones[j]->Lacuna != NULL) {
					Bones[i]->Lacuna->Absorb(Bones[j]->Lacuna);
				}
			}
		}
		break;
	}
	return L;
}


void TPolFigure::MonotonicSubdivision(){
	// Сделать подразбиение на монотонные рёбра
	for (TConnected* iComp = Components->first(); iComp != NULL; iComp = iComp->getNext()){
		LinkedListTail<TBone>* Update = new LinkedListTail<TBone>;
		list<TBone*> Skewed;

		for (TBone* iBone = iComp->Bones->first(); iBone != NULL; iBone = iBone->getNext()){
			iBone->DetermineType();
			if (iBone->type == Parabolic || iBone->type == Hyperbolic){
				QPair<QPointF, double> MinPoint = iBone->GetMinimum();
				if (MinPoint.second >= 0){
					// Бицикл содержит точку минимума
					if (MinPoint.second > iBone->org->r())
						MinPoint.second = iBone->org->r();
					if (MinPoint.second > iBone->dest->r())
						MinPoint.second = iBone->dest->r();

					TBone* nBone = new TBone();			// Создаём новое ребро

					TNode* nNode = new TNode();			// Создаём новый узел в точке минимума
					nNode->Disc = new TDisc(MinPoint.first.x(), MinPoint.first.y(), MinPoint.second);
					nNode->Sites[0] = iBone->Sites[0];
					nNode->Sites[1] = iBone->Sites[1];
					nNode->Bones[0] = iBone;
					nNode->Bones[1] = nBone;
					nNode->moveIntoTail(iComp->Nodes);

					nBone->Com = iBone->Com;
					nBone->org = nNode;
					nBone->dest = iBone->dest;
					nBone->Sites[0] = iBone->Sites[0];
					nBone->Sites[1] = iBone->Sites[1];
					nBone->type = iBone->type;
					nBone->moveIntoTail(Update);

					iBone->dest = iBone->org;
					iBone->org = nNode;
					bool stop = false;
					for (int i = 0; i < 3 && !stop; i++){
						if (nBone->dest->Bones[i] == iBone){
							nBone->dest->Bones[i] = nBone;
							stop = true;
						}
					}
					nBone->FillInfo();

					AllBones.push_back(nBone);
				}
			}
			iBone->FillInfo();
			AllBones.push_back(iBone);

			if (iBone->org->X() == iBone->dest->X() &&
				iBone->org->Y() == iBone->dest->Y() &&
				iBone->org->r() != iBone->dest->r())
			{
				iBone->org->Disc->Rad = iBone->dest->r();
				Skewed.push_back(iBone);
			}
		}
		bool repeat = true;
		while (repeat){
			repeat = false;
			for (auto iBone = Skewed.begin(); iBone != Skewed.end(); iBone++){
				if ((*iBone)->org->X() == (*iBone)->dest->X() &&
					(*iBone)->org->Y() == (*iBone)->dest->Y() &&
					(*iBone)->org->r() != (*iBone)->dest->r())
				{
					(*iBone)->org->Disc->Rad = (*iBone)->dest->r();
					repeat = true;
				}
			}
		}

		Update->moveAllElementsToAnotherTailEnd(iComp->Bones);
		delete Update;

		for (TNode* iNode = iComp->Nodes->first(); iNode != NULL; iNode = iNode->getNext()) {
			for (TBone** iBone = iNode->Bones; iBone - iNode->Bones < 3 && *iBone != NULL; iBone++) {
				for (TBone** jBone = iBone + 1; jBone - iNode->Bones < 3 && *jBone != NULL; jBone++) {
					if ((*iBone)->dest == iNode && (*jBone)->org == iNode) {
						TBone* Temp = *iBone;
						*iBone = *jBone;
						*jBone = Temp;
					}
				}
			}
			NodeKind a = iNode->Kind();
			AllNodes.push_back(iNode);
		}
	}
}


void TBone::DetermineType(){
	int found = 0;
	for (int i = 0; (i < 3) && (found < 2); i++){
		for (int j = 0; (j < 3) && (found < 2); j++){
			if (org->Sites[i] == dest->Sites[j]){
				Sites[found++] = org->Sites[i];
			}
		}
	}
	if (!Sites[0]->IsVertex() && !Sites[1]->IsVertex())
		type = Linear;
	else if (Sites[0]->IsVertex() && Sites[1]->IsVertex())
		type = Hyperbolic;
	else{
		type = Parabolic;
		if (!Sites[0]->IsVertex()){
			// Сайт-точка будет всегда идти первым
			TSite* temp = Sites[0];
			Sites[0] = Sites[1];
			Sites[1] = temp;
		}
	}
}


QPair<QPointF, double> TBone::GetMinimum(){
	if (type == Linear){
		return (QPair<QPointF, double>(QPointF(0, 0), -1));
		// Точка минимума линейного ребра не является внутренней
	}
	else{
		double r = org->r();
		double R = dest->r();
		if (r > R){
			double temp = r;
			r = R;
			R = temp;
		}
		double l = sqrt((dest->X() - org->X()) * (dest->X() - org->X()) +
			(dest->Y() - org->Y()) * (dest->Y() - org->Y()));
		if (type == Parabolic){
			double x1 = ((TVertex*)Sites[0])->p->X;
			double y1 = ((TVertex*)Sites[0])->p->Y;
			double x2 = ((TEdge*)Sites[1])->org->X;
			double y2 = ((TEdge*)Sites[1])->org->Y;
			double x3 = ((TEdge*)Sites[1])->dest->X;
			double y3 = ((TEdge*)Sites[1])->dest->Y;
			double rate = ((x1 - x2)*(x3 - x2) + (y1 - y2)*(y3 - y2)) /
				((x3 - x2)*(x3 - x2) + (y3 - y2)*(y3 - y2));
			// Отношение длины проекции к длине сайта-сегмента
			double x4 = x2 + rate * (x3 - x2);
			double y4 = y2 + rate * (y3 - y2);

			double cOrg = (org->X() - x1)*(y4 - y1) - (org->Y() - y1)*(x4 - x1);
			double cDest = (dest->X() - x1)*(y4 - y1) - (dest->Y() - y1)*(x4 - x1);
			if (cOrg * cDest < 0){
				// Точки по разные стороны оси симметрии параболы
				double x = (x1 + x4) / 2;
				double y = (y1 + y4) / 2;
				double r = sqrt((x1 - x4)*(x1 - x4) + (y1 - y4)*(y1 - y4)) / 2;

				if (x == org->X() && y == org->Y() || x == dest->X() && y == dest->Y())
					return QPair<QPointF, double>(QPointF(0, 0), -1);

				return (QPair<QPointF, double>(QPointF(x, y), r));
			}
			else{
				return (QPair<QPointF, double>(QPointF(0, 0), -1));
			}
		}
		else{
			double x1 = ((TVertex*)Sites[0])->p->X;
			double y1 = ((TVertex*)Sites[0])->p->Y;
			double x2 = ((TVertex*)Sites[1])->p->X;
			double y2 = ((TVertex*)Sites[1])->p->Y;

			double cOrg = (org->X() - x1)*(y2 - y1) - (org->Y() - y1)*(x2 - x1);
			double cDest = (dest->X() - x1)*(y2 - y1) - (dest->Y() - y1)*(x2 - x1);
			if (cOrg * cDest < 0){
				double x = (x1 + x2) / 2;
				double y = (y1 + y2) / 2;
				double r = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2)) / 2;
				return (QPair<QPointF, double>(QPointF(x, y), r));
			}
			else{
				return (QPair<QPointF, double>(QPointF(0, 0), -1));
			}
		}
	}
}


#define xMin	org->X()
#define yMin	org->Y()
#define rMin	org->r()
#define xMax	dest->X()
#define yMax	dest->Y()
#define rMax	dest->r()


void TBone::FillInfo(){

	if (org->r() > dest->r()){
		TNode* temp = org;
		org = dest;
		dest = temp;
	}
	double l = sqrt((xMax - xMin)*(xMax - xMin) + (yMax - yMin)*(yMax - yMin));

	if (type == Linear){
		t = sqrt(max(l*l - (rMax - rMin)*(rMax - rMin), 0.0));
		p = (l > 0) ? acos(min((rMax - rMin) / l, 1.0)) : M_PI_2;
		s = t * (rMin + rMax);
		if (l > 0){
			dx = (xMax - xMin) / l;
			dy = (yMax - yMin) / l;
		}
		else{
			double dx1 = ((TEdge*)Sites[0])->dest->X - ((TEdge*)Sites[0])->org->X;
			double dy1 = ((TEdge*)Sites[0])->dest->Y - ((TEdge*)Sites[0])->org->Y;
			double dl = sqrt(pow(dx1, 2) + pow(dy1, 2));
			dx1 /= dl;
			dy1 /= dl;
			double dx2 = ((TEdge*)Sites[1])->dest->X - ((TEdge*)Sites[1])->org->X;
			double dy2 = ((TEdge*)Sites[1])->dest->Y - ((TEdge*)Sites[1])->org->Y;
			double l2 = sqrt(pow(dx2, 2) + pow(dy2, 2));
			dx2 /= dl;
			dy2 /= dl;
			if (dx1*dx2 + dy1*dy2 > 0){
				dx = dx1 + dx2;
				dy = dy1 + dy2;
			}
			else{
				dx = dx1 - dx2;
				dy = dy1 - dy2;
			}
			dl = sqrt(pow(dx, 2) + pow(dy, 2));
			dx /= dl;
			dy /= dl;
		}
	}

	if (type == Parabolic){
		if (!Sites[0]->IsVertex()){
			// Сайт-точка будет всегда идти первым
			TSite* temp = Sites[0];
			Sites[0] = Sites[1];
			Sites[1] = temp;
		}
		double x1 = ((TVertex*)Sites[0])->p->X;
		double y1 = ((TVertex*)Sites[0])->p->Y;
		double x2 = ((TEdge*)Sites[1])->org->X;
		double y2 = ((TEdge*)Sites[1])->org->Y;
		double x3 = ((TEdge*)Sites[1])->dest->X;
		double y3 = ((TEdge*)Sites[1])->dest->Y;
		double a = ((x1 - x2)*(x3 - x2) + (y1 - y2)*(y3 - y2)) /
			((x3 - x2)*(x3 - x2) + (y3 - y2)*(y3 - y2));
		// Отношение длины проекции к длине сайта-сегмента
		double x4 = x2 + a * (x3 - x2);
		double y4 = y2 + a * (y3 - y2);
		x0 = (x1 + x4) / 2;
		y0 = (y1 + y4) / 2;
		p = sqrt((x1 - x4)*(x1 - x4) + (y1 - y4)*(y1 - y4));
		t = sqrt(max(2 * p*(rMax - p / 2), 0.0));
		s = (rMax + p) / 2 * t -
			(rMin + p) / 2 * sqrt(2 * p * max(rMin - p / 2, 0.0));
		dx = (y1 - y4) / p;
		dy = -(x1 - x4) / p;
		dir = 1;
		// Направляющий вектор оси OY
		if ((xMax - x0) * dx + (yMax - y0) * dy < 0){
			// Поворот от OX к OY производится по часовой стрелке
			dx = -dx;
			dy = -dy;
			dir = -1;
		}
	}

	if (type == Hyperbolic){
		double x1 = ((TVertex*)Sites[0])->p->X;
		double y1 = ((TVertex*)Sites[0])->p->Y;
		double x2 = ((TVertex*)Sites[1])->p->X;
		double y2 = ((TVertex*)Sites[1])->p->Y;
		x0 = (x1 + x2) / 2;
		y0 = (y1 + y2) / 2;
		p = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
		t = sqrt(max(rMax*rMax - p*p / 4, 0.0));
		s = p / 2 * t -
			p / 2 * sqrt(max(rMin*rMin - p*p / 4, 0.0));
		if (t > 0){
			dx = (xMax - x0) / t;
			dy = (yMax - y0) / t;
		}
		else{
			dx = (y1 - y2) / p;
			dy = (x2 - x1) / p;
		}
	}

	if (_isnan(s))
		throw("Square is not a number");
	return;

}


double TBone::Square(double rad){
	if (rad <= rMin){
		return s;
	}
	if (rad > rMax){
		return 0;
	}

	if (type == Linear){
		double q = (t > 0) ? t * (rMax - rad) / (rMax - rMin) : 0;
		return q * (rad + rMax);
	}
	if (type == Parabolic){
		double q = sqrt(max(2 * p*(rad - p / 2), 0.0));
		return (rMax + p) / 2 * t - (rad + p) / 2 * q;
	}

	// Бицикл - гиперболический
	return p/2 * (t - sqrt(max(rad*rad - p*p/4, 0.0)));

}


QPointF TBone::GetExtremePoint(double r)
{
	if (r <= rMin){
		return QPointF(xMin, yMin);
	}
	else if (r >= rMax){
		return QPointF(xMax, yMax);
	}
	double a, x, y;
	
	switch (type){
	case Linear:
		a = (r - rMin) / (rMax - rMin);
		x = xMin + a * (xMax - xMin);
		y = yMin + a * (yMax - yMin);
		break;
	case Parabolic:
		a = sqrt(max(2*p*(r-p/2), 0.0));
		x = x0 + a*dx - dir*(r-p/2)*dy;
		y = y0 + a*dy + dir*(r-p/2)*dx;
		break;
	case Hyperbolic:
		a = sqrt(max(r*r - p*p/4, 0.0));
		x = x0 + a*dx;
		y = y0 + a*dy;
		break;
	}
	if (_isnan(x) || _isnan(y))
		throw("Coordinate is not a number");
	return QPointF(x, y);
}


double TBone::SectorArea(double rad){
	double angle;
	switch (type){
	case Linear:
		angle = p;
		break;
	case Parabolic:
		angle = asin(min((p - rad) / rad, 1.0));
		angle = (angle + M_PI_2) / 2;
		break;
	case Hyperbolic:
		angle = asin(min((p/2) / rad, 1.0));
		break;
	}
	return angle * rad*rad;
}


std::list<QPair<int, int>> GeneratePairs(list<int> A, list<int> B) {
	std::list<QPair<int, int>> res;
	for (auto i = A.begin(); i != A.end(); i++) {
		for (auto j = B.begin(); j != B.end(); j++) {
			res.push_back(QPair<int, int>(std::min(*i, *j), std::max(*i, *j)));

		}
	}
	return res;
}


void Lacuna::TriangulateForViz(double rad, list<QPair<int, int>>& Edges, vector<double>& Vert) {
	Edges.clear();
	Vert.clear();
	set<QPair<double, double>> Ends;
	for (auto iBone : Bones) {
		if (iBone->dest->r() >= rad) {
			QPointF extr = iBone->GetExtremePoint(rad);
			Ends.insert({ extr.x(), extr.y() });
		}
	}
	for (auto Pt: Ends) {
		Vert.push_back(Pt.first);
		Vert.push_back(Pt.second);
	}

	if (Vert.size() / 2 > 1) {
		if (Vert.size() / 2 > 3) {
			double x_min = Vert[0];
			double x_max = Vert[0];
			double y_min = Vert[1];
			double y_max = Vert[1];
			for (int i = 1; i < Vert.size() / 2; i++) {
				x_min = std::min(x_min, Vert[2 * i]);
				x_max = std::max(x_max, Vert[2 * i]);
				y_min = std::min(y_min, Vert[2 * i + 1]);
				y_max = std::max(y_max, Vert[2 * i + 1]);
			}
			if (x_min == x_max || y_min == y_max) {
				int shift = (x_min == x_max) ? 1 : 0;
				vector<pair<double, int>> temp;
				for (int i = 0; i < Vert.size() / 2; i++) {
					temp.push_back(make_pair(Vert[2 * i + shift], i));
				}
				sort(temp.begin(), temp.end());
				for (int i = 0; i < temp.size() - 1; i++) {
					Edges.push_back({ temp[i].second, temp[i + 1].second });
				}
			}
			else {
				auto DT = delaunator::Delaunator(Vert);
				for (int i = 0; i < DT.triangles.size() / 3; i++) {
					vector<int> Idx = { (int)DT.triangles[3 * i + 0], (int)DT.triangles[3 * i + 1], (int)DT.triangles[3 * i + 2] };
					Edges.push_back({ Idx[0], Idx[1] });
					Edges.push_back({ Idx[0], Idx[2] });
					Edges.push_back({ Idx[1], Idx[2] });
				}
			}
		}
		else if (Vert.size() == 3) {
			Edges.push_back({ 0, 1 });
			Edges.push_back({ 0, 2 });
			Edges.push_back({ 1, 2 });
		}
		else {
			Edges.push_back({ 0, 1 });
		}
	}
}


double Lacuna::Triangulate(double rad) {
	list<QPair<int, int>> Edges;
	vector<QPointF> Ends;
	vector<double> Vert;
	vector<list<int>> Conn;
	vector<TBone*> Active;

	for (auto iBone = Bones.begin(); iBone != Bones.end(); iBone++) {
		(*iBone)->Lacuna = NULL;
		if ((*iBone)->dest->r() >= rad)
			Active.push_back(*iBone);	
	}
	Bones.clear();

	for (auto iBone = Active.begin(); iBone != Active.end(); iBone++) {
		QPointF Pt = (*iBone)->GetExtremePoint(rad);
		Ends.push_back(Pt);
		int idx = Vert.size() / 2;
		for (int i = 0; i < Vert.size() / 2 && idx == Vert.size() / 2; i++) {
			if (Pt.x() == Vert[2 * i] && Pt.y() == Vert[2 * i + 1]) {
				idx = i;
			}
		}
		if (idx == Vert.size() / 2) {
			Vert.push_back(Pt.x());
			Vert.push_back(Pt.y());
			Conn.push_back(list<int>());
		}
		Conn[idx].push_back(Ends.size() - 1);
	}

	list<QPair<int, int>> Add;
	if (Vert.size() / 2 > 1) {
		if (Vert.size() / 2 > 3) {
			double x_min = Vert[0];
			double x_max = Vert[0];
			double y_min = Vert[1];
			double y_max = Vert[1];
			for (int i = 1; i < Vert.size() / 2; i++) {
				x_min = std::min(x_min, Vert[2 * i]);
				x_max = std::max(x_max, Vert[2 * i]);
				y_min = std::min(y_min, Vert[2 * i + 1]);
				y_max = std::max(y_max, Vert[2 * i + 1]);
			}
			if (x_min == x_max || y_min == y_max) {
				int shift = (x_min == x_max) ? 1 : 0;
				vector<pair<double, int>> temp;
				for (int i = 0; i < Vert.size() / 2; i++) {
					temp.push_back(make_pair(Vert[2 * i + shift], i));
				}
				sort(temp.begin(), temp.end());
				for (int i = 0; i < temp.size() - 1; i++) {
					Add = GeneratePairs(Conn[temp[i].second], Conn[temp[i + 1].second]);
					Edges.splice(Edges.end(), Add);
				}
			}
			else {
				auto DT = delaunator::Delaunator(Vert);
				for (int i = 0; i < DT.triangles.size() / 3; i++) {
					vector<uint> Idx = { (uint)DT.triangles[3 * i + 0], (uint)DT.triangles[3 * i + 1], (uint)DT.triangles[3 * i + 2] };
					Add = GeneratePairs(Conn[Idx[0]], Conn[Idx[1]]); Edges.splice(Edges.end(), Add);
					Add = GeneratePairs(Conn[Idx[0]], Conn[Idx[2]]); Edges.splice(Edges.end(), Add);
					Add = GeneratePairs(Conn[Idx[1]], Conn[Idx[2]]); Edges.splice(Edges.end(), Add);
				}
			}
		}
		else if (Vert.size() / 2 == 3) {
			Add = GeneratePairs(Conn[0], Conn[1]); Edges.splice(Edges.end(), Add);
			Add = GeneratePairs(Conn[0], Conn[2]); Edges.splice(Edges.end(), Add);
			Add = GeneratePairs(Conn[1], Conn[2]); Edges.splice(Edges.end(), Add);
		}
		else {
			Add = GeneratePairs(Conn[0], Conn[1]); Edges.splice(Edges.end(), Add);
		}
	}

	auto iEdge = Edges.begin();
	while (iEdge != Edges.end())
	{
		QPointF P = QPointF((Ends[iEdge->first].x() + Ends[iEdge->second].x()) / 2,
							(Ends[iEdge->first].y() + Ends[iEdge->second].y()) / 2);
		QPointF C = Ends[iEdge->first];
		QPointF D = QPointF(Active[iEdge->first]->org->X(), Active[iEdge->first]->org->Y());
		TSite* SiteA = Active[iEdge->first]->Sites[0];
		TSite* SiteB = Active[iEdge->first]->Sites[1];
		QPointF A = SiteA->IsVertex() ? QPointF(((TVertex*)SiteA)->p->X, ((TVertex*)SiteA)->p->Y) : ProjectPoint(C, (TEdge*)SiteA);
		QPointF B = SiteB->IsVertex() ? QPointF(((TVertex*)SiteB)->p->X, ((TVertex*)SiteB)->p->Y) : ProjectPoint(C, (TEdge*)SiteB);

		if (CheckHalfPlane(P, D, QLineF(C.x(), C.y(), A.x(), A.y())) &&
			CheckHalfPlane(P, D, QLineF(C.x(), C.y(), B.x(), B.y())))
		{
			C = Ends[iEdge->second];
			D = QPointF(Active[iEdge->second]->org->X(), Active[iEdge->second]->org->Y());
			SiteA = Active[iEdge->second]->Sites[0];
			SiteB = Active[iEdge->second]->Sites[1];
			A = SiteA->IsVertex() ? QPointF(((TVertex*)SiteA)->p->X, ((TVertex*)SiteA)->p->Y) : ProjectPoint(C, (TEdge*)SiteA);
			B = SiteB->IsVertex() ? QPointF(((TVertex*)SiteB)->p->X, ((TVertex*)SiteB)->p->Y) : ProjectPoint(C, (TEdge*)SiteB);
			if (CheckHalfPlane(P, D, QLineF(C.x(), C.y(), A.x(), A.y())) &&
				CheckHalfPlane(P, D, QLineF(C.x(), C.y(), B.x(), B.y())))
			{
				iEdge++;
			}
			else {
				iEdge = Edges.erase(iEdge);
			}
		}
		else {
			iEdge = Edges.erase(iEdge);
		}
	}
	double res = 0;
	if (Edges.size() > 0) {
		double add = 0;
		vector<bool> Used(Active.size());
		for (iEdge = Edges.begin(); iEdge != Edges.end(); iEdge++) {
			add = PairSquare(Ends[iEdge->first], Ends[iEdge->second], rad);
			if (add > 0) {
				Used[iEdge->first] = Used[iEdge->second] = true;
				res += add;
			}
		}
		for (int i = 0; i < Used.size(); i++) {
			if (Used[i]) {
				Bones.push_back(Active[i]);
				Active[i]->Lacuna = this;
			}
		}
	}
	return res;
}


double PairSquare(QPointF A, QPointF B, double rad) {
	double xc = (A.x() + B.x()) / 2;
	double yc = (A.y() + B.y()) / 2;

	double l = sqrt(pow(B.x() - A.x(), 2) + pow(B.y() - A.y(), 2));
	if (l > 2 * rad) {
		return 0;
	}
	double h = sqrt(rad*rad - l * l / 4);
	double angle = asin(h / rad);
	double s = 2 * (angle*rad*rad - (l / 2)*h);
	if (_isnan(s))
		throw("Square is not a number");
	return s;
}


bool CheckHalfPlane(QPointF A, QPointF B, QLineF E) {
	double a = E.dy();
	double b = -E.dx();
	double c = -a*E.x1() - b*E.y1();
	double p = a * A.x() + b * A.y() + c;
	double q = a * B.x() + b * B.y() + c;
	return p * q >= 0;
}


QPointF ProjectPoint(QPointF A, TEdge* E) {
	double dx = E->dest->X - E->org->X;
	double dy = E->dest->Y - E->org->Y;
	double t = (dx*(A.x() - E->org->X) + dy * (A.y() - E->org->Y)) /
		(dx*dx + dy * dy);
	return QPointF(E->org->X + t * dx, E->org->Y + t * dy);
}