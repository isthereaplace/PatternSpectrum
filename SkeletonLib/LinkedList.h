#pragma once

#ifndef NULL
    #define NULL 0
#endif

namespace LinkedListNS {
	template<class T> class LinkedListTail;

	template<class T> class LinkedListElement
	{
        template<class U> friend class LinkedListTail;
    protected:
		LinkedListElement<T>* prev;
		LinkedListElement<T>* next;
		bool isTail;
	public:
		LinkedListElement()
		{
			prev = next = NULL;
			isTail = false;
        }
		virtual ~LinkedListElement()
		{
			removeFromCurrentList();
        }
		T* getNext()
		{
			if (next && next->isTail) {
				return NULL;
			}
			return (T*)next;
        }
		T* getPrev()
		{
			if (prev && prev->isTail) {
				return NULL;
			}
			return (T*)prev;
        }
		T* getNextLooped()
		{
			if (next && next->isTail) {
				return (T*)(next->next);
			}
			return (T*)next;
        }
		T* getPrevLooped()
		{
			if (prev && prev->isTail) {
				return (T*)(prev->prev);
			}
			return (T*)prev;
        }
		void removeFromCurrentList()
		{
			if (next) {
				// if in list
				next->prev = prev;
				prev->next = next;
				next = prev = NULL;
			}
        }
		void moveAsNextFor(LinkedListElement<T>* L)
		{
			removeFromCurrentList();
			prev = L;
			next = L->next;
			L->next = this;
			next->prev = this;
        }
		void moveAsPrevFor(LinkedListElement<T>* L)
		{
			removeFromCurrentList();
			next = L;
			prev = L->prev;
			L->prev = this;
			prev->next = this;
        }
		// before: el1 -> el2 -> TAIL -> el1
		// after: el1 -> el2 -> this -> TAIL -> el1
		void moveIntoTail(LinkedListTail<T>* L)
		{
			removeFromCurrentList();
			next = L;
			prev = L->prev;
			L->prev = this;
			prev->next = this;
        }
	};

	// el1 -> el2 -> TAIL -> el1
	template<class T> class LinkedListTail : public LinkedListElement < T >
	{
		inline LinkedListElement<T>* _first()
		{
			if (isEmpty()) {
				return NULL;
			}
            return this->next;
        }
		inline LinkedListElement<T>* _last()
		{
			if (isEmpty()) {
				return NULL;
			}
            return this->prev;
        }
	public:
		LinkedListTail()
		{
            this->prev = this->next = this;
            this->isTail = true;
        }
		virtual ~LinkedListTail()
		{

        }
		bool isEmpty()
		{
            return this->next == this;
        }
		T* first()
		{
			return (T*)_first();
        }
		T* last()
		{
			return (T*)_last();
        }
		int cardinal()
		{
			int count = 0;
			LinkedListElement<T>* curr;
			if (!isEmpty()) {
				curr = _first();
				while (curr != this) {
					count += 1;
					curr = curr->next;
				}
			}
			return count;
        }
		void clear()
		{
			LinkedListElement<T>* curr = _first();
			while (curr != NULL) {
				curr->removeFromCurrentList();
				curr = _first();
			}
        }
        void moveAllElementsToAnotherTailEnd(LinkedListTail<T>* S)
		{
			if (!isEmpty()) {
				LinkedListElement<T>* f = _first();
				LinkedListElement<T>* l = _last();
				// append to S list
				f->prev = S->prev;
				f->prev->next = f;
				// move tail to new end
				l->next = S;
				S->prev = l;
				// clear current tail, it's empty now
                this->prev = this->next = this;
			}
        }
	};
}

using namespace LinkedListNS;
