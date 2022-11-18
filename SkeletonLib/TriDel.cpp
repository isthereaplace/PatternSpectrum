#include "TriDel.h"

typedef unsigned char unsignedchar;

void d_Vertex(Vertex* v)
{
    qDebug() << "Vertex(" << v->p->X << ", " << v->p->Y << ")";
}

void d_Edge(Edge* v)
{
    qDebug() << "Edge(" << v->org->X << ", " << v->org->Y << " -> " << v->dest->X << ", " << v->dest->Y << ")";
}

void d_Element(Element* e)
{
    if (e->isVertex) {
        d_Vertex((Vertex*)e);
    }
    else {
        d_Edge((Edge*)e);
    }
}

void d_Triplet(Triplet* t)
{
    if (t != NULL) {
        d_Element(t->E1);
        d_Element(t->E2);
        d_Element(t->e3);
    } else {
        qDebug() << "NULL";
    }
}

void d_MapTriplet(LinkedListTail<Triplet>* v)
{
    qDebug() << "Triplets:" << endl;
    Triplet* c = v->first();
    while (c != NULL)
    {
        d_Triplet(c);
        qDebug() << "----";
        c = c->getNext();
    }
}

unsignedchar KindTriangle(Element *E1, Element *E2, Element *e3)
{
	unsignedchar result = '\0';
	unsigned char n = '\0';
	if (E1->IsVertex())
		n = 100;
	else
		n = 200;
	if (E2->IsVertex())
		n = n + 10;
	else
		n = n + 20;
	if (e3->IsVertex())
		n = n + 1;
	else
		n = n + 2;
	result = n;
	return result;
} /*KindTriangle*/


TDisc* CreateCircl(Element *E1, Element *E2, Element *e3)
{
	TDisc *result = NULL;
	Point *r0 = NULL, *r1 = NULL, *r2 = NULL, *r3 = NULL, *r4 = NULL, *r5 = NULL;
	TDisc *Cr = NULL;
	Edge *El1 = NULL, *El2 = NULL, *El3 = NULL;
	Vertex *t1 = NULL, *t2 = NULL, *t3 = NULL;
	unsigned char n = '\0';
	bool B = false;
	if ((E1 == NULL) || (E2 == NULL) || (e3 == NULL))
		result = NULL;
	else
	{
		Cr = NULL;
		El1 = NULL;
		El2 = NULL;
		El3 = NULL;
		t1 = NULL;
		t2 = NULL;
		t3 = NULL;
		n = KindTriangle(E1, E2, e3);
		switch (n)
		{
		case 222:
		{
			El1 = (Edge*)E1;
			El2 = (Edge*)E2;
			El3 = (Edge*)e3;
		}
			break;
		case 221:
		{
			El1 = (Edge*)E1;
			El2 = (Edge*)E2;
			t1 = (Vertex*)e3;
		}
			break;
		case 212:
		{
			El1 = (Edge*)e3;
			El2 = (Edge*)E1;
			t1 = (Vertex*)E2;
		}
			break;
		case 122:
		{
			El1 = (Edge*)E2;
			El2 = (Edge*)e3;
			t1 = (Vertex*)E1;
		}
			break;
		case 211:
		{
			El1 = (Edge*)E1;
			t1 = (Vertex*)E2;
			t2 = (Vertex*)e3;
		}
			break;
		case 121:
		{
			El1 = (Edge*)E2;
			t1 = (Vertex*)e3;
			t2 = (Vertex*)E1;
		}
			break;
		case 112:
		{
			El1 = (Edge*)e3;
			t1 = (Vertex*)E1;
			t2 = (Vertex*)E2;
		}
			break;
		case 111:
		{
			t1 = (Vertex*)E1;
			t2 = (Vertex*)E2;
			t3 = (Vertex*)e3;
		}
			break;
		}
		switch (n)
		{
		case 222:
		{
			r0 = El1->org;
			r1 = El1->dest;
			r2 = El2->org;
			r3 = El2->dest;
			r4 = El3->org;
			r5 = El3->dest;
			Cr = S3(r0, r1, r2, r3, r4, r5);
			if (Cr == NULL)
				Cr = S3(r2, r3, r4, r5, r0, r1);
			if (Cr == NULL)
				Cr = S3(r4, r5, r0, r1, r2, r3);
		}
			break;
		case 221: case 212: case 122:
		{
			r0 = El1->org;
			r1 = El1->dest;
			r2 = El2->org;
			r3 = El2->dest;
			r4 = t1->p;
			Cr = S2_P1(r0, r1, r2, r3, r4);
		}
			break;
		case 211: case 121: case 112:
		{
			r0 = El1->org;
			r1 = El1->dest;
			r2 = t1->p;
			r3 = t2->p;
			if ((r3 == r0) || (r3 == r1) || (r2 == r0) || (r2 == r1))
				Cr = S1_P1(r0, r1, r2, r3);
			else
				Cr = S1_P2(r0, r1, r2, r3);
		}
			break;
		case 111:
		{
			r0 = t1->p;
			r1 = t2->p;
			r2 = t3->p;
			Cr = P3(r0, r1, r2);
		}
			break;
		}
		if (Cr != NULL)
		{
			switch (n)
			{
			case 221: case 212: case 122:
				B = t1->ConformWithDisc(Cr);
				break;
			case 211: case 121: case 112:
				B = t1->ConformWithDisc(Cr) && t2->ConformWithDisc(Cr);
				break;
			case 111:
				B = t1->ConformWithDisc(Cr) && t2->ConformWithDisc(Cr) && t3->ConformWithDisc(Cr);
				break;
			default: /*222:*/
				B = true;
			}
			if (!B)
			{
				delete Cr;
				Cr = NULL;
			}
		}
		result = Cr;
	}
	return result;
}


void HullCorrection(Couple *&r)
{
	Triplet *Z = NULL;
	Element *E1 = NULL, *E2 = NULL, *e3 = NULL;
	Couple *p = NULL, *q = NULL, *S = NULL;
	Z = r->right;
	E1 = r->Frst;
	e3 = r->Scnd;
	E2 = Z->FollowingElement(E1);
	p = new Couple(E1, E2);
	q = new Couple(E2, e3);
	p->right = Z->AdjacentForElement(E1);
	q->right = Z->AdjacentForElement(E2);
	p->moveAsNextFor(r);
	q->moveAsNextFor(p);
	delete r;
	r = NULL;
	if (Z->ExtCouple != NULL)
	{ /*Подготовка к склеиванию рёбер*/
		if (p->EquCoples(Z->ExtCouple))
			S = p;
		else
			if (q->EquCoples(Z->ExtCouple))
				S = q;
			else
				S = NULL;
		Z->ExtCouple->PasteCouple = S;
		if (S != NULL)
			S->PasteCouple = Z->ExtCouple;
	}
	delete Z;
	Z = NULL;
	if ((p->right != NULL) && (p->right->Depth > 0))
		r = p;
	else
		if ((q->right != NULL) && (q->right->Depth > 0))
			r = q;
		else
			r = NULL;
}


bool InterElement(Element *e, TDisc *Circ)
/* Проверка попадания элемента E в круг Circ */
{
	bool result = false;
	Vertex * T = NULL;
	Edge * Ed = NULL;
	if (Circ == NULL)
		result = false;
	else
	{
		if (e->IsVertex())
		{
			T = (Vertex *)e;
			result = InterPoint(T->p, Circ);
		}
		else
		{
			Ed = (Edge *)e;
			result = InterEdge(Ed->org, Ed->dest, Circ);
		}
	}
	return result;
}


Triplet* InterCouples(Element *E, Triplet *t)
/* Проверка попадания элемента E в смежную пару симплексов */
{
	Triplet *result = NULL;
	Element *E1 = NULL, *E2 = NULL;
    Triplet *TS = NULL;
	Vertex *V = NULL;
	Edge *Ed1 = NULL, *Ed2 = NULL;
	Point *r = NULL;
	bool B = false;
	result = NULL;
	if (!E->IsVertex())
		return result;
	E1 = t->E1;
	for (int i = 1; i <= 3; i++)
	{
		E2 = t->FollowingElement(E1);
		if ((!E1->IsVertex()) || (!E2->IsVertex()))
		{
			TS = t->AdjacentForElement(E1);
			r = ((Vertex*)E)->p;
			if ((TS != NULL) && (t->Circ->Rad <= TS->Circ->Rad))
			{
				if (E1->IsVertex())
				{
					V = (Vertex*)E1;
					Ed1 = (Edge*)E2;
					B = InterPointEdge(r, V->p, Ed1->org, Ed1->dest, t->Circ, TS->Circ);
				}
				else
					if (E2->IsVertex())
					{
					V = (Vertex*)E2;
					Ed1 = (Edge*)E1;
					B = InterPointEdge(r, V->p, Ed1->org, Ed1->dest, t->Circ, TS->Circ);
					}
					else
					{
						Ed1 = (Edge*)E1;
						Ed2 = (Edge*)E2;
						B = InterTwoEdges(r, Ed1->org, Ed1->dest, Ed2->org, Ed2->dest, t->Circ, TS->Circ);
					}
				if (B)
				{
					result = TS;
					break;
				}
			}
		}
		E1 = E2;
	}
	return result;
}


bool TestOnIntersection(Couple *&r, Element *E)
/* Случай попадания элемента E в пустой круг симплекса, смежного с R */
{
	bool result = false;
	result = (r != NULL) && (r->right != NULL) && InterElement(E, r->right->Circ);
	if (result)
		HullCorrection(r);
	return result;
}


bool FullTestOnIntersection(Couple *&r, Element *E)
/* Случай попадания элемента E в пустой круг и бициклы симплекса,
								смежного с R */
{
	bool result = false;
	result = (r != NULL) && (r->right != NULL) && (E != r->right->E1) && (E != r->right->E2) && (E != r->right->e3) && (InterElement(E, r->right->Circ));
	if (result)
		HullCorrection(r);
	return result;
}

bool BestNewTriangle(LinkedListTail<Triplet> *ListTriplet, Couple *&r, Couple *&p, Couple *&q, Triplet *&Z1, Triplet *&Z2) /* Построение симплекса Делоне по цепочке пар P,R,Q по ЧС */ /* Z1,Z2 - построенные симплексы, z2#NIL если построено два симплекса */
{
	bool result = false;
	TDisc *CircP = NULL, *CircQ = NULL, *d0 = NULL;
	Point *S1 = NULL, *s0 = NULL;
	Triplet *Z = NULL;
	bool Pq = false;
	double A0 = 0.0, AP = 0.0, AQ = 0.0, dx = 0.0, dy = 0.0;
	Z1 = NULL;
	Z2 = NULL;
	Pq = false;
	result = false;
	CircP = NULL;
	if (p != NULL)
		CircP = CreateCircl(r->Frst, r->Scnd, p->Frst);
	CircQ = NULL;
	if (q != NULL)
		CircQ = CreateCircl(r->Frst, r->Scnd, q->Scnd);
	if ((CircP != NULL) && (CircQ != NULL))
	{
		if ((q != NULL) && InterElement(q->Scnd, CircP))
		{
			delete CircP;
			CircP = NULL;
		}
		else
			if ((p != NULL) && InterElement(p->Frst, CircQ))
			{
			delete CircQ;
			CircQ = NULL;
			}
	}
	if ((CircP == NULL) && (CircQ == NULL))
		return result;/* нельзя построить новый симплекс */
	if ((CircP != NULL) && (CircQ != NULL))
	{
		if (r->right == NULL)
			Pq = true; /*первый шаг -выбор любой*/
		else
		{
			d0 = r->right->Circ;
			if (!r->Frst->IsVertex())
			{  /*слева сайт-сегмент*/
				s0 = ((Edge*)r->Frst)->dest;
				S1 = ((Edge*)r->Frst)->org;
				dx = S1->X - s0->X;
				dy = S1->Y - s0->Y;
			}
			else
				if (!r->Scnd->IsVertex())
				{  /*справа сайт-сегмент*/
				s0 = ((Edge*)r->Scnd)->org;
				S1 = ((Edge*)r->Scnd)->dest;
				dx = S1->X - s0->X;
				dy = S1->Y - s0->Y;
				}
				else
				{ /*оба сайта - точки*/
					s0 = ((Vertex*)r->Frst)->p;
					S1 = ((Vertex*)r->Scnd)->p;
					dx = -(S1->Y - s0->Y);
					dy = (S1->X - s0->X);
				}
			A0 = (d0->X - s0->X) * dx + (d0->Y - s0->Y) * dy;
			AP = (CircP->X - s0->X) * dx + (CircP->Y - s0->Y) * dy;
			AQ = (CircQ->X - s0->X) * dx + (CircQ->Y - s0->Y) * dy;
			if ((AP < A0) && (AQ >= A0))
				Pq = false;
			else
				if ((AQ < A0) && (AP >= A0))
					Pq = true;
				else
					if (AP <= AQ)
						Pq = true;
					else
						if (AQ < AP)
							Pq = false;
		}
	}
	if ((CircP != NULL) && (CircQ == NULL))
		Pq = true;
	if ((CircQ != NULL) && (CircP == NULL))
		Pq = false;
	if (Pq) /* P главнее Q или круги совпадают */
	{
		Z1 = new Triplet(p->Frst, p->Scnd, r->Scnd, CircP);
		Z1->t2 = p->right;
		Z1->t3 = r->right;
		if (p->right != NULL)
			p->right->PasteWithTriplet(p->Scnd, Z1);
		if (r->right != NULL)
			r->right->PasteWithTriplet(r->Scnd, Z1);
		r->Frst = p->Frst;
		r->right = Z1;
		if (p->PasteCouple != NULL)
			p->PasteCouple->right = Z1;
		delete p;
		p = NULL;
		Z1->moveIntoTail(ListTriplet);
		result = true;
	}
	if (!Pq) /* Q главнее P или совпадают */
	{  /* Для случая Same=TRUE создаются два симплекса */
		Z = new Triplet(r->Frst, r->Scnd, q->Scnd, CircQ);
		Z->t2 = r->right;
		Z->t3 = q->right;
		if (r->right != NULL)
			r->right->PasteWithTriplet(r->Scnd, Z);
		if (q->right != NULL)
			q->right->PasteWithTriplet(q->Scnd, Z);
		r->Scnd = q->Scnd;
		r->right = Z;
		if (q->PasteCouple != NULL)
			q->PasteCouple->right = Z;
		delete q;
		q = NULL;
		Z->moveIntoTail(ListTriplet);
		if (Z1 == NULL)
			Z1 = Z;
		else
			Z2 = Z;
		result = true;
	}
  {  /* Уничтожить лишние круги */
	  if (!Pq && (CircP != NULL))
		  delete CircP;
	  if (Pq && (CircQ != NULL))
		  delete CircQ;
  }
	return result;
}


void Sew(LinkedListTail<Triplet> *ListTriplet, bool Direct, Element *&E, LinkedListTail<Couple> *&HullLeft, LinkedListTail<Couple> *&HullRight, Triplet *&TFirst, Triplet *&TLast, Triplet *&Z1, Triplet *&Z2, Couple *&p, Couple *&q, Couple *&r)
{  /*Сшивка - образование новых треугольников*/
	bool B, C;
	C = true;
	while (C)
	{
		B = true;
		while (B)
		{
			E = r->Scnd;
			if (Direct)
			{  /*Шов в прямом направлении*/
				p = HullLeft->last();
				q = HullRight->first();
			}
			else
			{  /*Шов в обратном направлении*/
				p = HullRight->last();
				q = HullLeft->first();
			}
			B = TestOnIntersection(p, E);
			E = r->Frst;
			B = B | TestOnIntersection(q, E);
		}
		C = BestNewTriangle(ListTriplet, r, p, q, Z1, Z2);

        //qDebug() << "$$$$ START";
        //d_MapTriplet(ListTriplet);
//        d_Triplet(Z1);
//        d_Triplet(Z2);
//        qDebug() << "$$$$ END";

		if ((TFirst == NULL) && (Z1 != NULL))
			TFirst = Z1;
		if (Z1 != NULL)
		{
			TLast = Z1;
			Z1->Depth = 1;
		}
		if (Z2 != NULL)
		{
			TLast = Z2;
			Z2->Depth = 1;
		}
	}
	if (r->right != NULL)
		r->right->Depth = 0;
}


void Merge(LinkedListTail<Triplet> *ListTriplet, LinkedListTail<Couple>* HullLeft, LinkedListTail<Couple>* HullRight, LinkedListTail<Couple>* ResHull) /* Слияние двух триангуляций */
{
	/*d_Output().clear(); d_Output().s() << "Left" << endl; d_Output().show();
	d_PrintHull(HullLeft);
	d_Output().clear(); d_Output().s() << "Right" << endl; d_Output().show();
	d_PrintHull(HullRight);*/

	Couple *t = NULL, *p = NULL, *q = NULL, *r = NULL, *r1 = NULL, *r2 = NULL;
	Element *E = NULL;
    bool B = false;
	Triplet *TFirst = NULL, *TLast = NULL, *Z1 = NULL, *Z2 = NULL;
	Element *ELeft = NULL, *ERight;
	Element *ELeftOrig = NULL, *ERightFin = NULL;
	TFirst = NULL;
	TLast = NULL;
	p = HullLeft->last();
	q = HullRight->first();
	ELeft = p->Scnd;
	ERight = q->Frst;
	p = HullLeft->first();
	q = HullRight->last();
	ELeftOrig = p->Frst;
	ERightFin = q->Scnd;
	r = new Couple(ELeft, ERight); /* Начальный мостик */
	if (ELeft->Cont == ERight->Cont)
	{ /*Сшиваются части одного контура*/
		Sew(ListTriplet, true, E, HullLeft, HullRight, TFirst, TLast, Z1, Z2, p, q, r);
		/*Объединение краёв в общую цепь ResHull*/
		while (!HullLeft->isEmpty())
		{
			t = HullLeft->first();
			t->moveIntoTail(ResHull);
			if (t->right != NULL)
				t->right->Depth = 0;
		}
		r->moveIntoTail(ResHull);
		while (!HullRight->isEmpty())
		{
			t = HullRight->first();
			t->moveIntoTail(ResHull);
			if (t->right != NULL)
				t->right->Depth = 0;
		}
		/*Контрольная расчистка при "загибающемся крае"*/
		E = ELeft;
		while (E != NULL)
		{
			B = true;
			while (B && (r != NULL) && (r->right != NULL))
				B = FullTestOnIntersection(r, E);
			if (E == ELeftOrig)
				E = NULL;
			else
				E = E->getPrev();
		}
		E = ERight;
		while (E != NULL)
		{
			B = true;
			while (B && (r != NULL) && (r->right != NULL))
				B = FullTestOnIntersection(r, E);
			if (E == ERightFin)
				E = NULL;
			else
				E = E->getNext();
		}
		if (r == NULL)
			TLast = NULL;
		else
			TLast = r->right;
	}
	else
	{ /*Сшиваются разные контура*/
		r2 = new Couple(r->Scnd, r->Frst);
		Sew(ListTriplet, true, E, HullLeft, HullRight, TFirst, TLast, Z1, Z2, p, q, r);
		r2->right = TFirst;
		r1 = r;
		if ((TLast != NULL) && TLast->ConsistsElements(ERight, ELeft))
		{ /* Склейка кольцевой триангулированной полосы */
			TFirst->PasteWithTriplet(ELeft, TLast);
			TLast->PasteWithTriplet(ERight, TFirst);
			delete r2;
			delete r1;
			return;
		}
		else
		{ /*Шов не кольцевой*/
		}
		/* IF R2.Frst.Cont.Internal<>R2.Scnd.Cont.Internal THEN
		 BEGIN R2.Destroy; R2:=NIL END; */
		/*Объединение краёв в общую цепь ResHull*/
		while (!HullLeft->isEmpty())
		{
			t = HullLeft->first();
			t->moveIntoTail(ResHull);
			if (t->right != NULL)
				t->right->Depth = 0;
		}
		r1->moveIntoTail(ResHull);
		while (!HullRight->isEmpty())
		{
			t = HullRight->first();
			t->moveIntoTail(ResHull);
			if (t->right != NULL)
				t->right->Depth = 0;
		}
		if (r2 != NULL)
			r2->moveIntoTail(ResHull);
	}
}

void TriangulateIteration(AdArrElem *Elements, TContour *Cont) /* Триангуляция элементов контура итерационным путем*/
{
	TriangleMap *Map = NULL, *Map1 = NULL, *Map2 = NULL;
	int i = 0;
	LinkedListTail<TriangleMap> *MapsList = NULL;
	Couple *p = NULL;
	/*Формирование начального списка карт*/
	MapsList = new LinkedListTail < TriangleMap > ;
	i = 0;
	while (i < Cont->NumbElem)
	{
		Map = new  TriangleMap();
		Map->NSite = 0;
		Map->moveIntoTail(MapsList);
		(new Couple(Elements[i], Elements[i + 1]))->moveIntoTail(Map->MapHull);
		i = i + 2;
		Map->NSite += 2;
		while ((i < Cont->NumbElem) && ((Vertex*)Cont->Elements[i])->Reflex)
		{
			(new Couple(Elements[i - 1], Elements[i]))->moveIntoTail(Map->MapHull);
			(new Couple(Elements[i], Elements[i + 1]))->moveIntoTail(Map->MapHull);
			i = i + 2;
			Map->NSite += 2;
		}
	}
	/*Итерации слияния*/
    while (MapsList->first() != MapsList->last())
	{
        //qDebug() << "la";
		Map = MapsList->first();
		while ((Map != NULL) && (Map->getNext() != NULL))
		{
			Map1 = Map;
			Map2 = Map->getNext();
			if (Map2 != NULL)
			{
				/* Слить Map1 и Map2 в Map */
				Map = new TriangleMap();
				LinkedListTail<Triplet> *ListTriplet = Map->MapTriplet;
				//d_PrintHull(Map1->MapHull);
				//d_PrintHull(Map2->MapHull);

				Merge(ListTriplet, Map1->MapHull, Map2->MapHull, Map->MapHull);

//                qDebug() << "#################### Map1";
//                d_MapTriplet(Map1->MapTriplet);
//                qDebug() << "#################### Map2";
//                d_MapTriplet(Map2->MapTriplet);

                Map1->MapTriplet->moveAllElementsToAnotherTailEnd(ListTriplet);
				Map2->MapTriplet->moveAllElementsToAnotherTailEnd(ListTriplet);

//                qDebug() << "#################### ListTriplet";
//                d_MapTriplet(ListTriplet);

				Map->NSite = Map1->NSite + Map2->NSite;
				Map->moveAsPrevFor(Map1);
				//d_PrintHull(Map->MapHull);
				Map = Map2->getNext();
				delete Map1;
                delete Map2;
			}
		}
	}
	Map = MapsList->first();

    //d_MapTriplet(Map->MapTriplet);

	if (Cont->Internal)
	{
		Map1 = Map;
		Map2 = Map;
		Map = new TriangleMap();
		LinkedListTail<Triplet> *ListTriplet = Map->MapTriplet;
		Merge(ListTriplet, Map1->MapHull, Map2->MapHull, Map->MapHull);
		Map1->MapTriplet->moveAllElementsToAnotherTailEnd(ListTriplet);
		delete Map1;
	}
	else
	{
		p = Map->MapHull->first();
		delete p;
	}
	delete Cont->Map;
	Cont->Map = Map;
	Map->removeFromCurrentList();
	//d_PrintHull(Map->MapHull);
	delete MapsList;
}

void TriangulateSegment(AdArrElem *Elements, LinkedListTail<Triplet> *ListTriplet, bool Internal, LinkedListTail<Couple> *ResHull) /* Триангуляция контура с двумя вершинами */
{
	TDisc *C;
	Triplet *t1, *t2;
	Vertex *V;
	if (Internal)
	{
		(new Couple(Elements[0], Elements[1]))->moveIntoTail(ResHull);
		(new Couple(Elements[1], Elements[2]))->moveIntoTail(ResHull);
		(new Couple(Elements[2], Elements[3]))->moveIntoTail(ResHull);
		(new Couple(Elements[3], Elements[0]))->moveIntoTail(ResHull);
	}
	else
	{
		V = (Vertex*)Elements[0];
		C = new TDisc(V->p->X, V->p->Y, 0);
		t1 = new Triplet(Elements[0], Elements[1], Elements[3], C);
		V = (Vertex*)Elements[2];
		C = new TDisc(V->p->X, V->p->Y, 0);
		t2 = new Triplet(Elements[2], Elements[3], Elements[1], C);
		t1->PasteWithTriplet(Elements[1], t2);
		t2->PasteWithTriplet(Elements[3], t2);
		t1->moveIntoTail(ListTriplet);
		t2->moveIntoTail(ListTriplet);
	}
}


void TriangulateContour(TContour *Cont)
{ /* Триангуляция отдельного контура */
	if (Cont->Map != NULL)
		delete Cont->Map;
	Cont->Map = new TriangleMap();
	AdArrElem *Elements = Cont->Elements;
	LinkedListTail<Triplet> *ListTriplet = Cont->Map->MapTriplet;
	if (Cont->NumbElem == 1)
		TriangulateSegment(Elements, ListTriplet, Cont->Internal, Cont->Map->MapHull);
	else
		TriangulateIteration(Elements, Cont);
}


void ProduceElements(Domain *Collect)
{
	TContour *S;
	int i = 0;
	i = 0;
	S = Collect->Boundary->first();
	while (S != NULL)
	{
		S->CreateElements(i);
		S = S->getNext();
	}
}


void RotationHull(TContour *Father, TContour *Son)
/* Смещение начала внешнего контура FatherHull таким образом, чтобы
  образовать начальный мостик с внутренним контуром SonHull */
{
	Couple *r, *q, *RMin;
	TDisc *Cr;
	Point *t1, *t2;
	Edge *E;
	double RadMin = 0.0;
	LinkedListTail<Couple> *FatherHull, *SonHull;
	RadMin = 100000;
	RMin = NULL;
	FatherHull = Father->Map->MapHull;
	SonHull = Son->Map->MapHull;
	q = SonHull->first();
	while (q->Frst != Son->WestElement)
	{
		q->moveIntoTail(SonHull);
		q = SonHull->first();
	}
	t1 = ((Vertex*)q->Frst)->p;
	t2 = new Point(t1->X, t1->Y + 1);
	E = new Edge(t1, t2);
	r = FatherHull->first();
	while (r != NULL)
	{
		Cr = CreateCircl(E, r->Frst, q->Frst);
		if ((Cr != NULL) && (Cr->Rad < RadMin))
		{
			RadMin = Cr->Rad;
			RMin = r;
		}
		if (Cr != NULL)
			delete Cr;
		r = r->getNext();
	}
	if (RMin != NULL)
	{
		r = FatherHull->first();
		while (RMin != r)
		{
			r->moveIntoTail(FatherHull);
			r = FatherHull->first();
		}
	}
	delete E;
	delete t2;
	//E.Destroy;
	//t2.Destroy;
}


void CuttingIn(TContour *Father, TContour *Son)
/* Врезка внутреннего контура SonHull во внешний FatherHull */
{
	Triplet *t, *TS;
	Couple *p, *q, *r, *S;
	Element *E1, *E2;
	Element *E;
	bool InCircle = false /*Попадание в круг или в бицикл*/, InBycicle = false;
	LinkedListTail<Triplet> *ListTriplet;
	LinkedListTail<Couple> *FatherHull;
	ListTriplet = Father->Map->MapTriplet;
	FatherHull = Father->Map->MapHull;
	E = Son->WestElement;
	t = ListTriplet->first();
	while ((t != NULL) && !InterElement(E, t->Circ))
		t = t->getNext();
	InCircle = (t != NULL);
	if (!InCircle) /* E не попадает в пустые круги симплекса */
	{
		/* Поиск бицикла (T,TS), содержащего E */
		TS = NULL;
		t = ListTriplet->first();
		while (t != NULL)
		{
			TS = InterCouples(E, t);
			if (TS != NULL)
				break;
			t = t->getNext();
		}
		InBycicle = (TS != NULL);
	}
	if (InCircle) /* Пробный элемент попадает в пустой круг симплекса */
	{
		p = new Couple(t->E1, t->E2);
		p->right = t->t2;
		p->moveIntoTail(FatherHull);
		q = new Couple(t->E2, t->e3);
		q->right = t->t3;
		q->moveIntoTail(FatherHull);
		r = new Couple(t->e3, t->E1);
		r->right = t->t1;
		r->moveIntoTail(FatherHull);
		if (t->ExtCouple != NULL)
		{ /*Подготовка к склеиванию рёбер*/
			if (p->EquCoples(t->ExtCouple))
				S = p;
			else
				if (q->EquCoples(t->ExtCouple))
					S = q;
				else
					if (r->EquCoples(t->ExtCouple))
						S = r;
					else
						S = NULL;
			t->ExtCouple->PasteCouple = S;
			if (S != NULL)
				S->PasteCouple = t->ExtCouple;
		}
		delete t;
		t = NULL;
		p = FatherHull->first();
		while (p != NULL)
		{
			t = p->right;
			if ((t != NULL) && InterElement(E, t->Circ))
			{
				q = p->getPrev();
				HullCorrection(p);
				if (q == NULL)
					p = FatherHull->first();
				else
					p = q->getNext();
			}
			else
				p = p->getNext();
		}
	}
	else
		if (InBycicle) /* Пробный элемент попадает в бицикл */
		{
		E1 = t->AdjacentForTriplet(TS);
		E2 = TS->FollowingElement(E1);
		t->BreakAdjacence(TS);
		TS->BreakAdjacence(t);
		p = new Couple(E1, E2);
		p->right = t;
		p->moveIntoTail(FatherHull);
		q = new Couple(E2, E1);
		q->right = TS;
		q->moveIntoTail(FatherHull);
		}
		else /*Никуда не попадает - это контакт двух внутренних контуров*/
		{
		}
}


void NextTopInHull(Couple *&PCur, Vertex *&VCur)
{
	while (!((PCur->Frst->IsVertex()) && (PCur->Frst != VCur)))
		PCur = PCur->getPrevLooped();
	VCur = (Vertex*)PCur->Frst;
}


bool InBack(Couple *r, Vertex *V)
{
	bool result = false;
	Point *P1, *P2, *q;
	PlanePosition Cl;
	P1 = ((Vertex*)r->Frst)->p;
	P2 = ((Vertex*)r->Scnd)->p;
	q = V->p;
	Cl = Classify(P1, P2, q);
	result = (Cl == RIGHT_POS) || (Cl == Between);
	return result;
}


bool CommonTangent(TContour *Father, TContour *Son)
/*Существование общего касательного ребра к оболочкам двух триангуляций.
  Если такое есть, то точки концы ребра выводятся в начало оболочек*/
{
	bool result = false;
	LinkedListTail<Couple> *FatherHull, *SonHull;
	Couple *p, *PF, *PS, *r;
	Vertex *VF, *VS;
	bool B = false, Sh = false;
	result = false;
	FatherHull = Father->Map->MapHull;
	SonHull = Son->Map->MapHull;
	Sh = false;  /*Признак первого сдвига нач. точки*/
	/*Начальное положение связующего ребра PF-PS*/
	PF = FatherHull->first();
	while (PF->Frst != Father->WestElement)
		PF = PF->getNext();
	PS = SonHull->first();
	while (PS->Frst != Son->WestElement)
		PS = PS->getNext();
	/*Проход вдоль контуров*/
	VF = (Vertex*)PF->Frst;
	VS = (Vertex*)PS->Frst;
	r = new Couple(VF, VS);
	NextTopInHull(PS, VS);
	B = true;
	while (B)
	{
		while (InBack(r, VS))
		{
			r->Scnd = VS;
			NextTopInHull(PS, VS);
		}
		NextTopInHull(PF, VF);
		if (Sh && (r->Frst == Father->WestElement))
		{
			delete r;
			return result;
		}
		if (InBack(r, VF))
		{
			Sh = true;
			r->Frst = VF;
		}
		else
			B = false;
	}
	/*Вращение на начальные точки*/
	if (r != NULL)
	{
		p = FatherHull->first();
		while (p->Frst != r->Frst)
		{
			p->moveIntoTail(FatherHull);
			p = FatherHull->first();
		}
		p = SonHull->first();
		while (p->Frst != r->Scnd)
		{
			p->moveIntoTail(SonHull);
			p = SonHull->first();
		}
		delete r;
	}
	result = true;
	return result;
}


void PasteTriangulations(TContour *Father, TContour *Son)
/* Вклейка триангуляции внутреннего контура в триангуляцию внешнего */
{
	LinkedListTail<Couple> *FatherHull, *SonHull, *ResHull, *Temp;
	Couple *p;

	////////////////  kk,ss: INTEGER;
	ResHull = new LinkedListTail < Couple > ;
	Temp = new LinkedListTail < Couple > ;
	LinkedListTail<Triplet> *ListTriplet = Father->Map->MapTriplet;
	FatherHull = Father->Map->MapHull;
	SonHull = Son->Map->MapHull;
	//////////////  MessageDlg('перед CommonTangent: '+IntToStr(AllocMemCount),mtWarning,[mbOk], 0);
	//////////////////  kk:=AllocMemCount;
	if (Father->Internal && CommonTangent(Father, Son))  /*Склейка двух внутренних контуров*/
	{
		//////////////    ss:=AllocMemCount;
		///////////////   MessageDlg('CommonTangent: '+IntToStr(kk)+'  '+IntToStr(ss),mtWarning,[mbOk], 0);
		Merge(ListTriplet, FatherHull, SonHull, ResHull);
		while (!ResHull->isEmpty())
		{
			p = ResHull->first();
			p->moveIntoTail(FatherHull);
		}
		////////////    ResHull.Destroy;
	}
	else
	{ /*Вклейка внутреннего во внешний либо в оболочку внутреннего*/
		while (!FatherHull->isEmpty())
		{
			p = FatherHull->first();
			if (p->right != NULL)
				p->right->ExtCouple = p;
			p->moveIntoTail(Temp);
		}
		CuttingIn(Father, Son);
		RotationHull(Father, Son);
		Merge(ListTriplet, FatherHull, SonHull, ResHull);
		if (ResHull->cardinal() == 1)
		{
			p = ResHull->first();
			delete p;
			///////////////      ResHull.Destroy;
		}
		while (!Temp->isEmpty())
		{
			p = Temp->first();
			if (p->right != NULL)
				p->right->ExtCouple = NULL;
			p->PasteCouple = NULL;
			p->moveIntoTail(FatherHull);
		}
		/////////////////    Temp.Destroy;
	}
	////////////  ss:=AllocMemCount;
	////////////  MessageDlg('PasteTriangulations: '+IntToStr(kk)+'  '
	////////////                +IntToStr(ss),mtWarning,[mbOk], 0);
	delete Temp;
	delete ResHull;
	while (!Son->Map->MapTriplet->isEmpty())
		Son->Map->MapTriplet->first()->moveIntoTail(ListTriplet);
	//////////////Son.Map.Destroy;
	////////////Son.Map:=NIL;
}


void Enumeration(Domain *Collect)
{
	TContour *S;
	Triplet *t;
	int n = 0;
	S = Collect->Boundary->first();
	while (S != NULL)
	{
		if (S->Map != NULL)
		{
			t = S->Map->MapTriplet->first();
			while (t != NULL)
			{
				t->Numb = n;
				n++;
				t = t->getNext();
			}
		}
		S = S->getNext();
	}
}


void CreateTriangulation(Domain *Collect)
{
	TContour *C;
	if ((Collect == NULL) || (Collect->Boundary == NULL))
		return;

	/* Триангуляция контуров */
	vector<TContour*> contours;
	C = Collect->Boundary->first();
	while (C != NULL)
	{
		contours.push_back(C);
		C = C->getNext();
	}

    //# pragma omp parallel for
    for (ulong i = 0; i < contours.size(); ++i) {
		TriangulateContour(contours[i]);
	}
	
	/* Склейка триангуляций контуров */
	C = Collect->Boundary->last();
	while (C != NULL)
	{
		if (C->Internal)
		{
			if (C->ClosestSite != NULL)
				PasteTriangulations(C->ClosestSite->Cont, C);
			else
				PasteTriangulations(C->Container, C);
			//////////////////      C.Map.Destroy;
			/////////////////      C.Map:=NIL;
		}
		C = C->getPrev();
	}
	Enumeration(Collect);
}
