#include "define.h"
#include <cassert>


template<typename T,int N>
class ROHashmap{
	template<typename U>
	struct CNode{
		int key;
		U data;
		CNode<U>* next;
		__device__ CNode() :key(-1), data(), next(0){}
		__device__ void insert_next(CNode<U>* next_node){
			next = next_node;
		}

		__device__ bool islast()const{
			return next == this;
		}


		__device__ CNode<U>* search_forward(int _k){
			if (key == _k){
				return this;
			}
			else if (islast()){
				return 0;
			}
			else{
				return next->search_forward(_k);

			}
		}
	};
	CNode<T> pool[N];
	CNode<T>* nodes[N];
	int count;
public:
	__device__ ROHashmap() :count(0),nodes(){}
	__device__ int _hash(int key){
		return key%N;
	}
	__device__ void reset(){
		memset(nodes, 0, sizeof(CNode<T>*)*N);
		count = 0;
	}

	__device__ void add(int key, T value){
		assert(count < N);
		CNode<T>* node = &pool[count++];
		node->key = key;
		node->data = value;
		int hash = _hash(key);
		if (nodes[hash] == 0){
			nodes[hash] = node;
			node->next = node;
		}
		else{
			node->next = nodes[hash];
			nodes[hash] = node;
		}
	}
	__device__ int size()const{
		return count;
	}
	__device__ T* at(int key){
		CNode<T>* fptr = nodes[_hash(key)];
		if (fptr == 0){
			return 0;
		}
		else{
			CNode<T>* src = fptr->search_forward(key);
			if (src == 0){
				return 0;
			}
			else{
				return &src->data;
			}
		}
	}
};
