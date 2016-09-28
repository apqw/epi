#pragma once
#include "define.h"
#include <cassert>

template<typename U, typename UKey_t = int>
class Node{
public:
	UKey_t key;
	U value;
	Node* next = nullptr;
	__device__ __host__ Node* find(UKey_t fkey){
		if (key == fkey){
			return this;
		}
		else if (next == nullptr){
			return nullptr;
		}
		else{
			return next->find(fkey);
		}
	}
	__device__ __host__ bool find_assign(UKey_t fkey, U v){
		if (key == fkey){
			return false;
		}
		else if (next == nullptr){
			next = new Node<U, UKey_t>(fkey, v);
			return true;
		}
		else{
			return next->find_assign(fkey, v);
		}
	}

	__device__ __host__ bool find_assign(Node* ptr){
		if (key == ptr->key){ //?
			return false;
		}
		else if (next == nullptr){
			next = ptr;
			return true;
		}
		else{
			return next->find_assign(ptr);
		}
	}

	__device__ __host__ bool remove(UKey_t fkey){
		if (next == nullptr){
			return false;
		}
		else if (next->key == fkey){
			Node* tmp = next->next;
			delete next;
			next = tmp;
			return true;
		}
		else{
			return next->remove(fkey);
		}
	}
	template<typename F>
	__device__ __host__ void forward_foreach(F&& lmbd){

		/*
		consider order
		*/
		if (next != nullptr){
			next->forward_foreach(lmbd);
		}
		lmbd(this);
	}
	__device__ __host__ Node(UKey_t k, U v) :key(k), value(v){}
};

template<typename T,typename Key_t=int>
class IntegerHashmap{


	Node<T, Key_t>** bucket_slot;
	int current_size=16;
	int bucket_count = 0;
public:
	__device__ __host__ int _hash(Key_t key)const{
		return key%current_size;
	}

	__device__ __host__ bool should_rehash()const{
		return bucket_count * 4 > current_size * 3;
	}
	__device__ __host__ T& at(Key_t key){
		const int hash = _hash(key);
		Node<T, Key_t>* nptr = bucket_slot[hash];
		assert(nptr != nullptr);
		return nptr->find(key)->value; //assume key will found
	}

	__device__ __host__ void emplace(Key_t key, T value){
		rehash_check();
		const int hash = _hash(key);
		Node<T, Key_t>** pptr = &bucket_slot[hash];
		if (*pptr == nullptr){
			*pptr = new Node<T, Key_t>(key, value);
			bucket_count++;
		}
		else if ((*pptr)->find_assign(key, value)){
			bucket_count++;
		}
		dbgprint("emplaced count:%d\n", bucket_count);
	}

	__device__ __host__ void _emplace_ptr(Node<T, Key_t>* ptr){
		const int hash = _hash(ptr->key);
		Node<T, Key_t>** pptr = &bucket_slot[hash];
		if (*pptr == nullptr){
			*pptr = ptr;
		}
		else {
			(*pptr)->find_assign(ptr);
		}
	}

	__device__ __host__ void remove(Key_t key){
		const int hash = _hash(key);
		Node<T, Key_t>** pptr = &bucket_slot[hash];
		assert(*pptr != nullptr);
		if ((*pptr)->key == key){
			Node<T, Key_t>* tmp = (*pptr)->next;
			delete *pptr;
			*pptr = tmp;
			bucket_count--;
		} else if((*pptr)->remove(key)){
			bucket_count--;
		}
		dbgprint("removed count:%d\n", bucket_count);
	}
/*
	Node<T, Key_t>** alloc_bucket_slot(int size){
		void* tmp = malloc(sizeof(Node<T, Key_t>*)*size);
		memset(tmp, 0, sizeof(Node<T, Key_t>*)*size);//nullptr
		return (Node<T, Key_t>**)tmp;
	}
	*/
	__device__ Node<T, Key_t>** 	alloc_bucket_slot(int size){
		void* tmp = malloc(sizeof(Node<T, Key_t>*)*size);
		memset(tmp, 0, sizeof(Node<T, Key_t>*)*size);//nullptr
		return (Node<T, Key_t>**)tmp;
	}

	__device__ __host__ void rehash_check(){
		if (should_rehash()){
			dbgprint("rehash start...\n");
			current_size *= 2;
			Node<T, Key_t>** old_bucket_slot = bucket_slot;
			bucket_slot = alloc_bucket_slot(current_size);
			for (int i = 0; i < current_size/2; i++){
				if (old_bucket_slot[i] != nullptr){
					old_bucket_slot[i]->forward_foreach([&](Node<T, Key_t>* cn){
						cn->next = nullptr;
						_emplace_ptr(cn);
					});
				}
			}
			free(old_bucket_slot);
			dbgprint("rehash end.\n");
		}
	}

	/*
	__device__ __host__ IntegerHashmap(int size = 16) :current_size(size), bucket_slot(alloc_bucket_slot(size)), bucket_count(0){
	}
	*/
	
	__device__ void init(int size = 16){
		current_size = size;
		bucket_slot = alloc_bucket_slot(size);
		bucket_count = 0;
	}
	
	__device__ __host__ void destroy(){
		
		for (int i = 0; i < current_size ; i++){
			if (bucket_slot[i] != nullptr){
				bucket_slot[i]->forward_foreach([&](Node<T, Key_t>* cn){
					delete cn;
				});
			}
		}
		free(bucket_slot);
		
	}
	
};
