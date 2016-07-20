/*
 * CHashmap.h
 *
 *  Created on: 2016/07/15
 *      Author: yasu7890v
 */

#ifndef CHASHMAP_H_
#define CHASHMAP_H_
#include "define.h"
#include <cassert>

template<int _N,typename T,int N=2*_N>
class CIndexHashmap {
private:
	template<typename U>
	struct CNode{
		int key;
		U data;
		CNode<U>* next;
		CNode<U>* prev;
		__device__ CNode() :key(-1), data(), next(0), prev(0){}
		__device__ void insert_next(CNode<U>* next_node){
			next->prev = next_node;
			next_node->next = next;
			next = next_node;
			next_node->prev = this;
		}

		__device__ void insert_prev(CNode<U>* prev_node){
			prev->next = prev_node;
			prev_node->prev = prev;
			prev = prev_node;
			prev_node->next = this;
		}

		__device__ void migrate(CNode<U>* dest){
			if (dest == this)return;
			dest->key = key;
			dest->data = data;
			next->prev = dest;
			prev->next = dest;
		}
		__device__ bool islast()const{
			return next == this;
		}

		__device__ bool isfirst()const{
			return prev == this;
		}
		__device__ void _fast_remove(CNode<U>* lptr){
			next->prev = prev;
			prev->next = next;
			lptr->migrate(this);
		}

		__device__ CNode<U>* search_forward(int _k){
			if (key == _k || next == this){
				return this;
			}
			else{
				return search_forward(_k);

			}
		}
	};

	CNode<T> bucket[N+1];
	CNode<T> first[N+1];
	CNode<T> last[N+1];
	int data_count;
public:
	__device__ int get_data_count()const{
		return data_count;
	}
	__device__ T& internal_data_itr(int idx){
		return bucket[idx].data;
	}
	__device__ CIndexHashmap():data_count(0){
		for(int i=0;i<N+1;i++){
			//bucket[i]=CDoubleLinkElem<T>();
			first[i].prev=&first[i];
			first[i].next=&last[i];
			last[i].next=&last[i];
			last[i].prev=&first[i];
		}

	}
	__device__ int _hash(int i)const{
		return i%(N+1);
	}
	__device__ void emplace(int key, T item){
		int didx=data_count++;
		bucket[didx].data=item;
		bucket[didx].key = key;
		int hash=_hash(key);
		last[hash].insert_prev(&bucket[didx]);
	}

	__device__ CNode<T>* _internal_emplace(int key){
		int didx = data_count++;
		bucket[didx].data = T();
		bucket[didx].key = key;
		int hash = _hash(key);
		last[hash].insert_prev(&bucket[didx]);
		return &bucket[didx];
	}
	__device__ T& operator[](int key){
		CNode<T>* fptr = first[_hash(key)].search_forward(key);
		if (fptr->islast()){
			return _internal_emplace(key)->data;
		}
		return fptr->data;
	}
	__device__ void remove(int key){
		CNode<T>* fptr = first[_hash(key)].search_forward(key);
		if (!fptr->islast()||data_count>0){
			
			fptr->_fast_remove(&bucket[--data_count]);
		}
	}

	__device__ void internal_remove(int idx){
		bucket[idx]._fast_remove(&bucket[--data_count]);
	}
	//virtual ~CHashmap();
};

#endif /* CHASHMAP_H_ */
