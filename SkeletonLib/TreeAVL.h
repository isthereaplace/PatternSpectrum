#pragma once

/*Модуль для работы с АВЛ-дервьями. Алгоритмы взяты из книги Н.Вирта
  "Алгоритмы+структуры данных-=программы"*/

#include "LinkedList.h"


class TAVL;
class TKey;
class TKnot;


typedef int Balance;
typedef int(*KeyCompare)(void*, void*);
/*Если P1>P2, то 1, если P1<P2, то 1, если P1=P2, то 0*/


class TKey {
	friend class TAVL;
	friend class TKnot;
public:
	int Val;
public:
	TKey();
};


class TKnot : public LinkedListElement < TKnot > {
	friend class TAVL;
	friend class TKey;
public:
	void* key;
	TKnot* left, *right;
	Balance bal;
	TKnot();
	virtual ~TKnot();
};


class TAVL {
	friend class TKey;
	friend class TKnot;
public:
	TKnot* Root;
	LinkedListTail<TKnot>* Ruler;
	KeyCompare Compare;
	int Count;
	TAVL(KeyCompare CompareFunc);
	virtual ~TAVL();
public:
	void Search(void* X, TKnot*& p, bool& h);
	/*Включение ключа x в поддерево с корнем p*/
	/*h=TRUE - это поддерево выросло по высоте*/
	void Delete(void* X, TKnot*& p, bool& h);
	/*Удаление ключа x из поддерева с корнем p*/
	void Insert(void* X);
	/*Вставка в дерево ключа, которого там заведомо нет*/
	void Remove(void* X);
	/*Удаление из дерева ключа */
	void* MinKey();
	/*Минимальный ключ в дереве*/
	TKnot* SearchKey(void* X, TKnot* p);
	/*Поиск узла с ключем x в поддереве с корнем p*/
	void AddKnot(TKnot* Knot, TKnot*& p, bool& h);
	/*Включение узла Knot в поддерево с корнем p*/
	/*h=TRUE - это поддерево выросло по высоте*/
	void SubtracKnot(TKnot* Knot, TKnot*& p, bool& h);
	/*Включение узла Knot из поддерева с корнем p*/
	TKnot* FindKey(void* key);
	void InsertKnot(TKnot* Knot);
	void RemoveKnot(TKnot* Knot);
	TKnot* AboveKnot(TKnot* Knot);
	TKnot* BelowKnot(TKnot* Knot);
	void ReplaceKey(TKnot* Knot1, TKnot* Knot2);
	TKnot* MinKnot();
	bool empty();
private:
	void Del(TKnot*& r, bool& h, TKnot*& q);
};
