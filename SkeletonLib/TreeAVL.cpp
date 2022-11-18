#include "TreeAVL.h"


TKey::TKey()
	: Val(0)
{
}


TKnot::TKnot()
	: key(0),
	left(0),
	right(0)
{
	bal = 0;
}


TKnot::~TKnot()
{
}


TAVL::TAVL(KeyCompare CompareFunc)
	: Root(0),
	Ruler(new LinkedListTail<TKnot>),
	Compare(CompareFunc),
	Count(0)
{
}


TAVL::~TAVL()
{
	TKnot* Knot = 0;
	Root = 0;
	while (!Ruler->isEmpty())
	{
		Knot = Ruler->first();
		delete Knot;
	}
	delete Ruler;
}


void TAVL::Search(void* X, TKnot*& p, bool& h)
/* Установка элемента в АВЛ-дерево*/
{
	TKnot* P1 = 0, *P2 = 0;
	bool NewKnot = false;
	if (p == 0) /* Включить новый узел */
	{
		p = new TKnot;
		h = true;
		p->key = X;
		//       p.count:=1;
		if (Ruler->isEmpty())
			p->moveIntoTail(Ruler);
	}
	else
		if (Compare(p->key, X) == 1) /*x встраивается в левое поддерево*/
		{
		NewKnot = (p->left == 0);
		Search(X, p->left, h);
		if (NewKnot && (p->left != 0))
			p->left->moveAsPrevFor(p);
		if (h) /* Выросла левая ветвь */
		{
			switch (p->bal)
			{
			case 1:
			{
				p->bal = 0;
				h = false;
			}
				break;
			case 0:
				p->bal = -1;
				break;
			case -1: /* Балансировка */
			{
				P1 = p->left;
				if (P1->bal == -1) /* Однократный LL-поворот */
				{
					p->left = P1->right;
					P1->right = p;
					p->bal = 0;
					p = P1;
				}
				else /* Двойной LR-поворот */
				{
					P2 = P1->right;
					P1->right = P2->left;
					P2->left = P1;
					p->left = P2->right;
					P2->right = p;
					if (P2->bal == -1)
						p->bal = 1;
					else
						p->bal = 0;
					if (P2->bal == 1)
						P1->bal = -1;
					else
						P1->bal = 0;
					p = P2;
				}
				p->bal = 0;
				h = false;
			}
				break;
			}
		}
		}
		else
			if (Compare(p->key, X) == -1) /*x встраивается в правое поддерево*/
			{
		NewKnot = (p->right == 0);
		Search(X, p->right, h);
		if (NewKnot && (p->right != 0))
			p->right->moveAsNextFor(p);
		if (h) /* Выросла правая ветвь */
		{
			switch (p->bal)
			{
			case -1:
			{
				p->bal = 0;
				h = false;
			}
				break;
			case 0:
				p->bal = 1;
				break;
			case 1: /* Балансировка */
			{
				P1 = p->right;
				if (P1->bal == 1) /* Однократный RR-поворот */
				{
					p->right = P1->left;
					P1->left = p;
					p->bal = 0;
					p = P1;
				}
				else /* Двойной RL-поворот */
				{
					P2 = P1->left;
					P1->left = P2->right;
					P2->right = P1;
					p->right = P2->left;
					P2->left = p;
					if (P2->bal == 1)
						p->bal = -1;
					else
						p->bal = 0;
					if (P2->bal == -1)
						P1->bal = 1;
					else
						P1->bal = 0;
					p = P2;
				}
				p->bal = 0;
				h = false;
			}
				break;
			}
		}
			}
			else; /*p.Key=x */ /*x совпал с ключем узла p*/
			//       p.count:=p.count+1;
}


void balanceL(TKnot*& p, bool& h)
/* Левая ветвь стала короче */
{
	TKnot* P1 = 0, *P2 = 0;
	Balance b1, b2;
	switch (p->bal)
	{
	case -1:
		p->bal = 0;
		break;
	case 0:
	{
		p->bal = 1;
		h = false;
	}
		break;
	case 1: /* Балансировка */
	{
		P1 = p->right;
		b1 = P1->bal;
		if (b1 >= 0) /* Однократный RR-поворот */
		{
			p->right = P1->left;
			P1->left = p;
			if (b1 == 0)
			{
				p->bal = 1;
				P1->bal = -1;
				h = false;
			}
			else
			{
				p->bal = 0;
				P1->bal = 0;
			}
			p = P1;
		}
		else /* Двойной LR-поворот */
		{
			P2 = P1->left;
			b2 = P2->bal;
			P1->left = P2->right;
			P2->right = P1;
			p->right = P2->left;
			P2->left = p;
			if (b2 == 1)
				p->bal = -1;
			else
				p->bal = 0;
			if (b2 == -1)
				P1->bal = 1;
			else
				P1->bal = 0;
			p = P2;
			P2->bal = 0;
		}
	}
		break;
	}
}


void balanceR(TKnot*& p, bool& h)
/* Правая ветвь стала короче */
{
	TKnot* P1 = 0, *P2 = 0;
	Balance b1, b2;
	switch (p->bal)
	{
	case 1:
		p->bal = 0;
		break;
	case 0:
	{
		p->bal = -1;
		h = false;
	}
		break;
	case -1: /* Балансировка */
	{
		P1 = p->left;
		b1 = P1->bal;
		if (b1 <= 0) /* Однократный LL-поворот */
		{
			p->left = P1->right;
			P1->right = p;
			if (b1 == 0)
			{
				p->bal = -1;
				P1->bal = 1;
				h = false;
			}
			else
			{
				p->bal = 0;
				P1->bal = 0;
			}
			p = P1;
		}
		else /* Двойной LR-поворот */
		{
			P2 = P1->right;
			b2 = P2->bal;
			P1->right = P2->left;
			P2->left = P1;
			p->left = P2->right;
			P2->right = p;
			if (b2 == -1)
				p->bal = 1;
			else
				p->bal = 0;
			if (b2 == 1)
				P1->bal = -1;
			else
				P1->bal = 0;
			p = P2;
			P2->bal = 0;
		}
	}
		break;
	}
}


void TAVL::Del(TKnot*& r, bool& h, TKnot*& q)
{
	if (r->right != 0)
	{
		Del(r->right, h, q);
		if (h)
			balanceR(r, h);
	}
	else
	{
		q->key = r->key;
		//         q.count:=r.count;
		q = r;
		r = r->left;
		h = true;
	}
}


void TAVL::Delete(void* X, TKnot*& p, bool& h)
/* Удаление элемента из АВЛ-дерева */
{
	TKnot* q = 0;
	if (p == 0) {} /* Ключа в дереве нет */
	else
		if (Compare(p->key, X) == 1)
		{
		Delete(X, p->left, h);
		if (h)
			balanceL(p, h);
		}
		else
			if (Compare(p->key, X) == -1)
			{
		Delete(X, p->right, h);
		if (h)
			balanceR(p, h);
			}
			else /* Исключение p */
			{
				q = p;
				if (q->right == 0)
				{
					p = q->left;
					h = true;
				}
				else
					if (q->left == 0)
					{
					p = q->right;
					h = true;
					}
					else
					{
						Del(q->left, h, q);
						if (h)
							balanceL(p, h);
					}
				delete q;
			}
}


TKnot* TAVL::FindKey(void* key)
/*поиск элемента в АВЛ-дереве*/
{
	TKnot* result = 0;
	result = SearchKey(key, Root);
	return result;
}


TKnot* TAVL::SearchKey(void* X, TKnot* p)
/*Поиск узла с ключем x в поддереве с корнем p*/
{
	TKnot* result = 0;
	int Comp = 0;
	if (p == 0) /* такого ключа нет */
		result = 0;
	else
	{
		Comp = Compare(p->key, X);
		if (Comp == 1)
			/*x находится в левом поддереве*/
			result = SearchKey(X, p->left);
		else
			if (Comp == -1)
				/*x находится в правом поддереве*/
				result = SearchKey(X, p->right);
			else /*x совпал с ключем p*/
				result = p;
	}
	return result;
}


void TAVL::InsertKnot(TKnot* Knot)
{
	bool h = false;  /* Вставка узла в АВЛ-дерево*/
	AddKnot(Knot, Root, h);
}

void TAVL::AddKnot(TKnot* Knot, TKnot*& p, bool& h)
/*Включение узла Knot в поддерево с корнем p*/
/*h=TRUE - это поддерево выросло по высоте*/
{
	TKnot* P1 = 0, *P2 = 0;
	bool NewKnot = false;
	void* X = 0;
	X = Knot->key;
	if (p == 0) /* Включить новый узел */
	{
		p = Knot;
		h = true;
		if (Ruler->isEmpty())
			p->moveIntoTail(Ruler);
	}
	else
		if (Compare(p->key, X) == 1) /*x встраивается в левое поддерево*/
		{
		NewKnot = (p->left == 0);
		Search(X, p->left, h);
		if (NewKnot && (p->left != 0))
			p->left->moveAsPrevFor(p);
		if (h) /* Выросла левая ветвь */
		{
			switch (p->bal)
			{
			case 1:
			{
				p->bal = 0;
				h = false;
			}
				break;
			case 0:
				p->bal = -1;
				break;
			case -1: /* Балансировка */
			{
				P1 = p->left;
				if (P1->bal == -1) /* Однократный LL-поворот */
				{
					p->left = P1->right;
					P1->right = p;
					p->bal = 0;
					p = P1;
				}
				else /* Двойной LR-поворот */
				{
					P2 = P1->right;
					P1->right = P2->left;
					P2->left = P1;
					p->left = P2->right;
					P2->right = p;
					if (P2->bal == -1)
						p->bal = 1;
					else
						p->bal = 0;
					if (P2->bal == 1)
						P1->bal = -1;
					else
						P1->bal = 0;
					p = P2;
				}
				p->bal = 0;
				h = false;
			}
				break;
			}
		}
		}
		else
			if (Compare(p->key, X) == -1) /*x встраивается в правое поддерево*/
			{
		NewKnot = (p->right == 0);
		Search(X, p->right, h);
		if (NewKnot && (p->right != 0))
			p->right->moveAsNextFor(p);
		if (h) /* Выросла правая ветвь */
		{
			switch (p->bal)
			{
			case -1:
			{
				p->bal = 0;
				h = false;
			}
				break;
			case 0:
				p->bal = 1;
				break;
			case 1: /* Балансировка */
			{
				P1 = p->right;
				if (P1->bal == 1) /* Однократный RR-поворот */
				{
					p->right = P1->left;
					P1->left = p;
					p->bal = 0;
					p = P1;
				}
				else /* Двойной RL-поворот */
				{
					P2 = P1->left;
					P1->left = P2->right;
					P2->right = P1;
					p->right = P2->left;
					P2->left = p;
					if (P2->bal == 1)
						p->bal = -1;
					else
						p->bal = 0;
					if (P2->bal == -1)
						P1->bal = 1;
					else
						P1->bal = 0;
					p = P2;
				}
				p->bal = 0;
				h = false;
			}
				break;
			}
		}
			}
			else; /*x совпал с ключем узла p*/
}


void TAVL::SubtracKnot(TKnot* Knot, TKnot*& p, bool& h)
/*Исключение узла Knot из поддерева с корнем p*/
{
	TKnot* q = 0;
	void* X = 0;
	X = Knot->key;
	if (p == 0) {} /* Ключа в дереве нет */
	else
		if (Compare(p->key, X) == 1)
		{
		SubtracKnot(Knot, p->left, h);
		if (h)
			balanceL(p, h);
		}
		else
			if (Compare(p->key, X) == -1)
			{
		SubtracKnot(Knot, p->right, h);
		if (h)
			balanceR(p, h);
			}
			else /* Исключение p */
			{
				q = p;
				if (q->right == 0)
				{
					p = q->left;
					h = true;
				}
				else
					if (q->left == 0)
					{
					p = q->right;
					h = true;
					}
					else
					{
						Del(q->left, h, q);
						if (h)
							balanceL(p, h);
					}
				q->removeFromCurrentList();
			}
}


void TAVL::RemoveKnot(TKnot* Knot)
{
	bool h = false;
	SubtracKnot(Knot, Root, h);
}


TKnot* TAVL::AboveKnot(TKnot* Knot)
{
	TKnot* result = 0;
	result = Knot->getNext();
	return result;
}


TKnot* TAVL::BelowKnot(TKnot* Knot)
{
	TKnot* result = 0;
	result = Knot->getPrev();
	return result;
}


void TAVL::ReplaceKey(TKnot* Knot1, TKnot* Knot2)
/*Замена ключа в Knot1 на ключ из Knot2*/
{
	Knot1->key = Knot2->key;
}


TKnot* TAVL::MinKnot()
{
	TKnot* result = 0;
	result = Ruler->first();
	return result;
}


void TAVL::Insert(void* X)
/*Вставка в дерево ключа, которого там заведомо нет*/
{
	bool h = false;
	Search(X, Root, h);
	Count++;
}


void TAVL::Remove(void* X)
/*Удаление из дерева ключа */
{
	bool h = false;
	Delete(X, Root, h);
	Count--;
}


void* TAVL::MinKey()
/*Минимальный ключ в дереве с удалением*/
{
	void* result = 0;
	TKnot* Kn = 0;
	bool h = false;
	Kn = Ruler->first();
	if (Kn != 0)
	{
		result = Kn->key;
		Delete(Kn->key, Root, h);
		Count--;
	}
	else
		result = 0;
	return result;
}


bool TAVL::empty()
{
	bool result = false;
	result = Ruler->isEmpty();
	return result;
}
