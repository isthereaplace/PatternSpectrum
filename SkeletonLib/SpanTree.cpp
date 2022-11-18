#include "SpanTree.h"

const int NMem = 20;   /*Константа для отладочной визуализации*/


enum EventKind {
	EastCape,
	WestCape,
	PassSite,
	AddDisc,
	RemoveDisc
};
/*типы событий*/
typedef Element TSite;
typedef Point TPoint;

class TEventPoint : public LinkedListElement < TEventPoint > /* событие-изменение статуса */ {
public:
	TSite *Site;
	TSite *St1, *St2;
	/* В зависимости от типа события Kind сайты St1,St2 - это:
	   если Kind=EastCape или WestCape, то St1 - нижний, St2 - верхний;
	   если Kind=PassSite, то St1 - правее точки Pt, St2 - левее Pt*/
	TPoint* Pt; /*точка, с которой связано событие*/
	EventKind Kind;   /*тип события*/
	TEventPoint(TSite *S, TSite *SPred, TSite* Ssuc, bool Bubble); // Create
	TEventPoint(TSite *S, TDisc *d); // CreateByDisc
	virtual ~TEventPoint();
	int  ComparePos(TEventPoint * Ev);
	/* положение относительно Ev */
	TEventPoint() {}
};


class TPriorLine {
public:
	LinkedListTail<TEventPoint> * EvList;    /*список статических событий "концы отрезков"*/
	TAVL * EvTree;    /*дерево динамических событий*/
	int Memb[NMem + 1];  /*массив событий для отладки*/
	double XPos;                       /*абсциса заметающей прямой*/
	void  Insert(TEventPoint * Ev);  /*вставить событие в очередь*/
	void  Remove(TEventPoint * Ev);  /*удалить событие из очереди*/
	TEventPoint *  MinEvent();    /*ближайшее событие*/
	TPriorLine(Domain* Collect);/*создать очередь*/
	virtual ~TPriorLine();
	bool  empty();            /* очередь пустая? */
	void  View();                     /*визуализация*/
	TPriorLine() {}
};


class TSiteRecord {
public:
	TSite* Site; /* сайт, пересекающийся с заметающей прямой */
	TSiteRecord(TSite *St);
	TSiteRecord() {}
	virtual ~TSiteRecord() {};
};


class TDiscRecord : public TSiteRecord {
public:
	TDisc* Disc; /* пузырь */
	TSite* BornSite; /*образующий сайт для диска*/
	TEventPoint * EventOwn; /*запланированное событие*/
	TDiscRecord(TSite* St); // CreateTop
	void  Correction(TSite *St); /*коррекция пузыря */
	virtual ~TDiscRecord();
	TDiscRecord() {}
};


class TStatus {
public:
	TAVL * StTree;  /* упорядоченное множество отрезков в статусе*/
	int Memb[NMem + 1];
	/*упорядоченный список сайтов для отладки*/
	void  Insert(TSiteRecord * SRec);  /*вставка записи в статус*/
	void  Remove(TSiteRecord * SRec);  /*удаление записи из статуса*/
	void  Replace2(TSiteRecord * SRec, TSiteRecord * SBelow, TSiteRecord * SAbove);
	/*замена записи двумя другими*/
	void  Replace1(TSiteRecord * SRec, TSiteRecord * SAdj);
	/*замена одной записи другой*/
	void  Neighbours(TSiteRecord * SRec, TSiteRecord * SAdj);
	/*обработка записей, ставших соседними*/
	TSiteRecord *  Find(TSite* St);
	TSiteRecord *  Above(TSiteRecord * SRec);
	/* сайт выше */
	TSiteRecord *  Below(TSiteRecord * SRec);
	/* сайт ниже */
	TStatus(TPriorLine * PrLine); // Create
	virtual ~TStatus();
	void  WorkEvent(TEventPoint * EvP);
	/* обработка очередного события */
	void  TestUp(TSiteRecord * SRec);
	/*коррекция смежных кругов, лежащих выше*/
	void  TestDown(TSiteRecord * SRec);
	/*коррекция смежных кругов, лежащих ниже*/
	void  View();             /*визуализация*/
	TStatus() {}
};



TPriorLine * PriorLine;  /*очередь событий*/
TStatus * Status;          /*статус заметающей прямой*/

/* Описание методов и подпрограмм */


TEventPoint::TEventPoint(TSite *S, TSite *SPred, TSite *Ssuc, bool Bubble)
{
	Site = S;
	Kind = EastCape;  // default value, ok?
	bool PP = false, SP = false;
	double Prod = 0.0;
	Vertex* V = (Vertex*)S;
	Pt = V->p;
	Edge *EdP = (Edge*)SPred;
	Edge *EdS = (Edge*)Ssuc;
	PP = EdP->WestDirect();
	SP = EdS->WestDirect();
	if (PP == SP)
	{ /*проходной узел*/
		Kind = PassSite;
		if (EdP->org->X > EdS->dest->X)
		{
			St1 = EdP;
			St2 = EdS;
		}
		else
		{
			St1 = EdS;
			St2 = EdP;
		}
	}
	else
	{
		Prod = (EdP->dest->X - EdP->org->X) * (EdS->dest->Y - EdS->org->Y) - (EdP->dest->Y - EdP->org->Y) * (EdS->dest->X - EdS->org->X);
		if (PP && (!SP))
		{ /*западный мыс*/
			if (Bubble)
				Kind = AddDisc;
			else
				Kind = WestCape;
			if (Prod < 0)
			{
				St1 = EdP;
				St2 = EdS;
			}
			else
			{
				St1 = EdS;
				St2 = EdP;
			}
		}
		else
		{ /*западный мыс*/
			Kind = EastCape;
			if (Prod > 0)
			{
				St1 = EdP;
				St2 = EdS;
			}
			else
			{
				St1 = EdS;
				St2 = EdP;
			}
		}
	}
}


TEventPoint::TEventPoint(TSite *S, TDisc *D)
	: Site(S),
	Pt(new TPoint(D->X - D->Rad, D->Y)),
	Kind(RemoveDisc)
{
}


TEventPoint::~TEventPoint()
{
	if (Kind == RemoveDisc)
		delete Pt;
	// todo check:  inherited::Destroy;
}


int TEventPoint::ComparePos(TEventPoint* Ev)
{
	// 0 => equal, 1 => Item1 > Item2, -1 => Item1 < Item2
	int result = 1;
	if (Pt->X < Ev->Pt->X)
		result = -1; /*левый конец левее*/
	else
		if (Pt->X == Ev->Pt->X)
		{
		if (Pt->Y < Ev->Pt->Y)
			result = -1;/*абсциссы совпадают*/
		else
			if (Pt->Y == Ev->Pt->Y)
				result = 0;
		}
	return result;
}


int CompareEvents(void* Item1, void* Item2)
{
	TEventPoint *Ev1 = (TEventPoint*)Item1;
	TEventPoint *Ev2 = (TEventPoint*)Item2;
	int result = -Ev1->ComparePos(Ev2); // descending
	return result;
}

bool CompareForSTL(void* Item1, void* Item2) {
	int result = CompareEvents(Item1, Item2);

	// convert to C++ "less" function
	if (result == -1) {
		return true;
	}
	else {
		return false;
	}
}


TPriorLine::TPriorLine(Domain* Collect)
{
	EvList = new LinkedListTail<TEventPoint>(),
		EvTree = new TAVL(CompareEvents),
		XPos = 0.0;
	vector<TEventPoint*> MyList;
	int i = 0;
	TSite *E, *E1, *E2;
	TEventPoint* EvP = NULL;
	// inherited::Create();
	/*Составление общего списка событий*/
	TContour* S = Collect->Boundary->first();
	while (S != NULL)
	{
		/*событие AddDisc для начального сайта-точки*/
		E1 = S->Elements[S->NumbElem - 1];
		i = 0;
		while (i < S->NumbElem)
		{
			E = S->Elements[i];
			E2 = S->Elements[i + 1];
			EvP = new TEventPoint(E, E1, E2, (i == 0) && S->Internal);
			MyList.push_back(EvP);
			i = i + 2;
			E1 = E2;
		}
		S = S->getNext();
	}
	/*Упорядочение списка сайтов справа налево*/
	sort(MyList.begin(), MyList.end(), CompareForSTL);
	/*Формирование очереди событий*/
	for (int i = 0; i <= int(MyList.size()) - 1; i++)
	{
		EvP = MyList[i];
		EvP->moveIntoTail(EvList);
	}
}


void TopologySort(Domain* Collect)
{
	TContour *St, *S1;
	/*упорядочить контура по касанию пузырей*/
	TContour* S = Collect->Boundary->first();
	while (S != NULL)
	{
		S1 = S->getNext();
		if (S->ClosestSite != NULL)
		{
			St = S->ClosestSite->Cont;
			S->Container = St;
			if (St->MySons == NULL)
				St->MySons = new LinkedListTail < TContour > ;
			S->moveIntoTail(St->MySons);
		}
		S = S1;
	}
	/*установить вложенность контуров*/
	St = Collect->Boundary->first();
	while (St != NULL)
	{
		if (St->MySons != NULL)
		{
			while (!St->MySons->isEmpty())
			{
				S = St->MySons->first();
				if (St->Internal == S->Internal)
					S->Container = St->Container;
				else
					S->Container = St;
				S->moveIntoTail(Collect->Boundary);
			}
			delete St->MySons;
			St->MySons = NULL;
		}
		St = St->getNext();
	}
}


void SpaningTree(Domain* Collect)
{
	TEventPoint* EvP = NULL;
	/*построить очередь событий*/
	PriorLine = new TPriorLine(Collect);
	/*построить статус заметающей прямой*/
	Status = new TStatus(PriorLine);
	/*провести прямую справа налево*/
	while (!PriorLine->empty())
	{
		//      PriorLine.View;
		EvP = PriorLine->MinEvent();
		PriorLine->XPos = EvP->Pt->X;
		/*обработка очередного события*/

		// debug
		//d_Output().clear(); d_Output().s() << "EvP " << PriorLine->XPos << " " << EvP->Pt->X << " " << EvP->Pt->Y << endl; d_Output().show();

		try
		{
			Status->WorkEvent(EvP);
		}
		catch (...)
		{
            qDebug() << "Ошибка на сайте " << EvP->Site->NEl;
			break;
		}
		delete EvP;
	}
	delete Status;
	delete PriorLine;
	TopologySort(Collect);
}


void TPriorLine::Insert(TEventPoint* Ev)
{
	EvTree->Insert(Ev);
	//    View;
}


void TPriorLine::Remove(TEventPoint* Ev)
{
	EvTree->Remove(Ev);
	//    View;
}


TEventPoint* TPriorLine::MinEvent()
{
	TEventPoint* result = NULL;
	TEventPoint* EvT = NULL, *EvL = NULL;
	double XT = 0.0, XL = 0.0;
	TKnot* Kn;
	XT = -10000000;
	XL = -10000000;
	Kn = EvTree->Ruler->first();
	if (Kn != NULL)
	{
		EvT = (TEventPoint*)Kn->key;
		XT = EvT->Pt->X;
	}
	EvL = EvList->first();
	if (EvL != NULL)
		XL = EvL->Pt->X;
	if (XT >= XL)
		result = (TEventPoint*)EvTree->MinKey();
	else
	{
		EvL->removeFromCurrentList();
		result = EvL;
	}
	//    View;
	return result;
}


TPriorLine::~TPriorLine()
{
	while (!EvList->isEmpty())
		delete EvList->first();
	delete EvList;
	delete EvTree;
}


bool TPriorLine::empty()
{
	bool result = false;
	result = (EvList->isEmpty()) && (EvTree->Ruler->isEmpty());
	return result;
}


void TStatus::WorkEvent(TEventPoint* EvP)
{
	TSiteRecord* SR = NULL, *SRA = NULL, *SRB = NULL, *SAb = NULL, *SBe = NULL, *SRab = NULL;
	TDiscRecord* SD = NULL;
	switch (EvP->Kind)
	{
	case EastCape:
	{ /*добавление восточного мыса в сатус*/
		SR = new TSiteRecord(EvP->Site);
		Insert(SR);
		TestUp(SR);
		TestDown(SR);
		SAb = new TSiteRecord(EvP->St2);
		SBe = new TSiteRecord(EvP->St1);
		Replace2(SR, SBe, SAb);
		TestUp(SAb);
		TestDown(SBe);
		TestDown(SAb);
		TestUp(SBe);
		delete SR;
	}
		break;
	case WestCape:
	{ /*удаление западного мыса из статуса*/
		SR = Find(EvP->St1);
		SRB = Above(SR);
		Remove(SRB);
		SRA = new TSiteRecord(EvP->Site);
		Replace1(SR, SRA);
		delete SR;
		delete SRB;
		SR = SRA;
		TestUp(SR);
		TestDown(SR);
		SRA = Above(SR);
		SRB = Below(SR);
		Remove(SR);
		delete SR;
		TestUp(SRB);
		TestDown(SRA);
	}
		break;
	case PassSite:
	{ /*проходной промежуточный сайт*/
		SR = Find(EvP->St1);
		SRab = new TSiteRecord(EvP->Site);
		Replace1(SR, SRab);
		delete SR;
		SR = SRab;
		TestUp(SR);
		TestDown(SR);
		SRab = new TSiteRecord(EvP->St2);
		Replace1(SR, SRab);
		delete SR;
		SR = SRab;
		TestUp(SR);
		TestDown(SR);
	}
		break;
	case AddDisc:
	{ /* западный мыс - с пузырем */
		SR = Find(EvP->St1);
		SRB = Above(SR);
		if (SRB != NULL)
		{
			Remove(SRB);
			delete SRB;
		}
		SRab = new TSiteRecord(EvP->Site);
		Replace1(SR, SRab);
		delete SR;
		SR = SRab;
		TestUp(SR);
		TestDown(SR);
		SD = /*CreateTop*/ new TDiscRecord(EvP->Site);
		Replace1(SR, SD);
		delete SR;
		PriorLine->Insert(SD->EventOwn);
		SRA = Above(SD);
		SRB = Below(SD);
		if (SRA != NULL)
			SD->Correction(SRA->Site);
		if (SRB != NULL)
			SD->Correction(SRB->Site);
	}
		break;
	case RemoveDisc:
	{
		SR = Find(EvP->Site);
		SD = (TDiscRecord*)SR;
		SD->Site->Cont->ClosestSite = SD->BornSite;
		if (SD->BornSite != NULL)
			SD->Site->Cont->Container = SD->BornSite->Cont;
		SRA = Above(SR);
		SRB = Below(SR);
		Remove(SR);
		delete SD;
		Neighbours(SRA, SRB);
	}
		break;
	}
}


TDiscRecord::TDiscRecord(TSite* St)
	: TSiteRecord(St),
	Disc(new TDisc(-1000 - ((Vertex*)(St))->p->Y, 0, ((Vertex*)(St))->p->X, true)),
	BornSite(NULL),
	EventOwn(/*CreateByDisc*/ new TEventPoint(St, Disc))
{
}


TDiscRecord::~TDiscRecord()
{
	delete Disc;
	// todo check:  inherited::Destroy;
}


void TDiscRecord::Correction(TSite* St)
{
	TPoint *r, *r0, *r1;
	double ax = 0.0, ay = 0.0, bx = 0.0, by = 0.0, S = 0.0, M = 0.0, t = 0.0, lam = 0.0; /*коррекция пузыря */
	t = -1;
	if (St->IsVertex())
	{
		r0 = ((Vertex*)St)->p;
		r = ((Vertex*)Site)->p;
		bx = r->X - r0->X;
		by = r->Y - r0->Y;
		if (bx > 0)
			t = 0.5 * (bx + double(Sqr(by)) / bx);
		else
			if ((bx == 0) && (by > 0))
				t = -2;
	}
	else
	{
		r = ((Vertex*)Site)->p;
		r0 = ((Edge*)St)->org;
		r1 = ((Edge*)St)->dest;
		ax = r1->X - r0->X;
		ay = r1->Y - r0->Y;
		bx = r->X - r0->X;
		by = r->Y - r0->Y;
		S = ax * by - ay * bx;
		M = Sqr(ax) + Sqr(ay);
		if (ax == 0)
			t = -S / (2 * ay);
		else
			t = (S * ay + abs(S) * sqrt(M)) / Sqr(ax);
		lam = (ax * (bx - t) + ay * by) / M;
		if ((lam <= 0) || (lam >= 1))
			t = -1;
		//      IF (lam<0) OR (lam>1) THEN t:=-1;
	}
	if ((t > 0) && (Disc->HalfPlane || (t < Disc->Rad)))
	{
		Disc->X = r->X - t;
		Disc->Y = r->Y;
		Disc->Rad = t;
		Disc->HalfPlane = false;
		BornSite = St;
		if (EventOwn != NULL)
		{ /*перепланирование события*/
			PriorLine->Remove(EventOwn);
			EventOwn->Pt->X = Disc->X - t;
			EventOwn->Pt->Y = Disc->Y;
			PriorLine->Insert(EventOwn);
		}
	}
	else
		if ((t == -2) && Disc->HalfPlane && (BornSite == NULL))
			BornSite = St;
}


bool Lower(TSite* S1, TSite* S2)
/*TRUE если S1 ниже S2 в точке x*/
{
	bool result = false;
	Vertex *V1, *V2;
	Edge *E1, *E2;
	TPoint *P1;
	PlanePosition Pos;
	double Y1 = 0.0, Y2 = 0.0, X = 0.0;
	if (S1->IsVertex())
	{
		V1 = (Vertex*)S1;
		if (S2->IsVertex())
		{ /*оба сайта - точки*/
			V2 = (Vertex*)S2;
			result = (V1->p->Y < V2->p->Y) || ((V1->p->Y == V2->p->Y) && (S1->NEl < S2->NEl));
		}
		else /* первый сайт - точка, а второй - ребро */
		{
			E2 = (Edge*)S2;
			X = PriorLine->XPos;
			Y1 = V1->p->Y;
			if (abs(E2->org->X - E2->dest->X) < 0.1)
				Y2 = double((E2->org->Y + E2->dest->Y)) / 2;
			else
				Y2 = E2->org->Y + (E2->dest->Y - E2->org->Y) * (X - E2->org->X) / (E2->dest->X - E2->org->X);
			result = (Y1 < Y2) || ((Y1 == Y2) && (S1->NEl < S2->NEl));
		}
	}
	else
	{
		E1 = (Edge*)S1;
		if (S2->IsVertex())
		{
			V2 = (Vertex*)S2;
			X = PriorLine->XPos;
			Y2 = V2->p->Y;
			if (abs(E1->org->X - E1->dest->X) < 0.1)
				Y1 = double((E1->org->Y + E1->dest->Y)) / 2;
			else
				Y1 = E1->org->Y + (E1->dest->Y - E1->org->Y) * (X - E1->org->X) / (E1->dest->X - E1->org->X);
			result = (Y1 < Y2) || ((Y1 == Y2) && (S1->NEl < S2->NEl));
		}
		else
		{
			E2 = (Edge*)S2;
			X = PriorLine->XPos;
			if (abs(E1->org->X - E1->dest->X) < 0.1)
				Y1 = double((E1->org->Y + E1->dest->Y)) / 2;
			else
				Y1 = E1->org->Y + (E1->dest->Y - E1->org->Y) * (X - E1->org->X) / (E1->dest->X - E1->org->X);
			if (abs(E2->org->X - E2->dest->X) < 0.1)
				Y2 = double((E2->org->Y + E2->dest->Y)) / 2;
			else
				Y2 = E2->org->Y + (E2->dest->Y - E2->org->Y) * (X - E2->org->X) / (E2->dest->X - E2->org->X);
			if (Y1 != Y2)
				result = (Y1 < Y2) || ((Y1 == Y2) && (S1->NEl < S2->NEl));
			else
			{
				if ((E1->org == E2->org) || (E1->dest == E2->org))
					P1 = E2->dest;
				else
					P1 = E2->org;
				Pos = Classify((E1->org), (E1->dest), P1);
				result = (E1->WestDirect() != (LEFT_POS == Pos));
			}
		}
	}
	return result;
}


TSiteRecord::TSiteRecord(TSite* St)
{
	Site = St;
}


void TStatus::Insert(TSiteRecord *SRec)
{
	StTree->Insert(SRec);
	//    View;
}


void TStatus::Remove(TSiteRecord *SRec)
{
	if (SRec != NULL)
		StTree->Remove(SRec);
	//    View;
}


void TStatus::Replace2(TSiteRecord *SRec, TSiteRecord *SBelow, TSiteRecord *SAbove)
/*замена записи двумя другими*/
{
	StTree->Remove(SRec);
	//    View;
	StTree->Insert(SBelow);
	//    View;
	StTree->Insert(SAbove);
	//    View;
}


void TStatus::Replace1(TSiteRecord *SRec, TSiteRecord *SAdj)
/*замена одной записи другой*/
{
	StTree->Remove(SRec);
	StTree->Insert(SAdj);
	//    View;
}


void  TStatus::Neighbours(TSiteRecord * SRec, TSiteRecord * SAdj)
{
	TDiscRecord * SD; /*обработка записей, ставших соседними*/
	if (((SRec == NULL) || (SAdj == NULL)) || (dynamic_cast<TDiscRecord *>((SRec)) && dynamic_cast<TDiscRecord *>((SAdj))))
		return;
	else
		if (dynamic_cast<TDiscRecord *>(SRec))
		{
		SD = (TDiscRecord *)SRec;
		SD->Correction((SAdj->Site));
		}
		else
			if (dynamic_cast<TDiscRecord *>(SAdj))
			{
		SD = (TDiscRecord *)SAdj;
		SD->Correction((SRec->Site));
			}
}


void  TStatus::TestUp(TSiteRecord * SRec)
/*Коррекция смежных кругов, лежащих выше*/
{
	TDiscRecord * SD;
	TSiteRecord * SAdj;
	if (SRec == NULL)
		return;
	SAdj = Above(SRec);
	while ((SAdj != NULL) && dynamic_cast<TDiscRecord *>((SAdj)))
	{
		SD = (TDiscRecord *)SAdj;
		SD->Correction((SRec->Site));
		if (SD->EventOwn->Pt->X < PriorLine->XPos)
			SAdj = Above(SAdj);
		else
			SAdj = NULL;
	}
}


void  TStatus::TestDown(TSiteRecord * SRec)
/*Корркция смежных кругов, лежащих ниже*/
{
	TDiscRecord * SD;
	TSiteRecord * SAdj;
	if (SRec == NULL)
		return;
	SAdj = Below(SRec);
	while ((SAdj != NULL) && dynamic_cast<TDiscRecord *>((SAdj)))
	{
		SD = (TDiscRecord *)SAdj;
		SD->Correction((SRec->Site));
		if (SD->EventOwn->Pt->X < PriorLine->XPos)
			SAdj = Below(SAdj);
		else
			SAdj = NULL;
	}
}


TSiteRecord *  TStatus::Find(TSite* St)
{
	TSiteRecord * result;
	TSiteRecord * Fict;
	TKnot * Kn;
	Fict = new TSiteRecord(St);
	Kn = StTree->FindKey(Fict);
	if (Kn == NULL) {
        qDebug() << "Не найден сайт " << St->NEl;
	}
	if (Kn == NULL)
		result = NULL;
	else
		result = (TSiteRecord *)Kn->key;
	delete Fict;
	return result;
}


TSiteRecord *  TStatus::Above(TSiteRecord * SRec)
{
	TSiteRecord * result;
	TKnot * Kn;
	TSiteRecord * SR;
	Kn = StTree->FindKey(SRec);
	Kn = StTree->AboveKnot(Kn);
	if (Kn != NULL)
		SR = (TSiteRecord *)Kn->key;
	else
		SR = NULL;
	result = SR;
	return result;
}


TSiteRecord *  TStatus::Below(TSiteRecord * SRec)
{
	TSiteRecord * result;
	TKnot * Kn;
	Kn = StTree->FindKey(SRec);
	Kn = StTree->BelowKnot(Kn);
	if (Kn != NULL)
		result = (TSiteRecord *)Kn->key;
	else
		result = NULL;
	return result;
}


int  CompareRecords(void* p1, void* p2)
{
	int result;
	TSiteRecord * RecP1, *RecP2;
	RecP1 = (TSiteRecord *)p1;
	RecP2 = (TSiteRecord *)p2;
	if (RecP1->Site == RecP2->Site)
		result = 0;
	else
		if (Lower(RecP1->Site, RecP2->Site))
			result = -1;
		else
			result = 1;
	return result;
}


TStatus::TStatus(TPriorLine *)
{
	StTree = new TAVL(CompareRecords);
}


TStatus::~TStatus()
{
	TSiteRecord * SR;
	TDiscRecord * SD;
	TKnot * Kn;
	Kn = StTree->Ruler->first();
	while (Kn != NULL)
	{
		SR = (TSiteRecord *)Kn->key;
		if (dynamic_cast<TDiscRecord*>(SR))
		{
			SD = (TDiscRecord *)SR;
			delete SD;
		}
		else
			delete SR;
		Kn = Kn->getNext();
	}
	delete StTree;
}


void  TStatus::View() /*визуализация*/
{
	TKnot * Kn;
	TSiteRecord * SR;
	int i;
	Kn = StTree->Ruler->first();
	for (i = 0; i != NMem; i++)
		Memb[i] = 0;
	i = 0;
	while ((Kn != NULL) && (i <= NMem))
	{
		SR = (TSiteRecord *)Kn->key;
		Memb[i] = SR->Site->NEl;
		Kn = Kn->getNext();
		i++;
	}
}


void  TPriorLine::View() /*визуализация*/
{
	TKnot * Kn;
	TEventPoint * EvP;
	int i;
	Kn = EvTree->Ruler->first();
	for (i = 0; i != NMem; i++)
		Memb[i] = 0;
	i = 0;
	while ((Kn != NULL) && (i <= NMem))
	{
		EvP = (TEventPoint*)Kn->key;
		Memb[i] = EvP->Site->NEl;
		Kn = Kn->getNext();
		i++;
	}
}
