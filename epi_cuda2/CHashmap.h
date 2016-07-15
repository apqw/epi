/*
 * CHashmap.h
 *
 *  Created on: 2016/07/15
 *      Author: yasu7890v
 */

#ifndef CHASHMAP_H_
#define CHASHMAP_H_
/*
template<typename T>
struct CNode{
	T data;
	CNode<T>* next;
	CNode<T>* prev;
	CNode():data(),next(0),prev(0){}
	void insert_next(CNode<T>* next_node){
		next->prev=next_node;
		next_node->next=next;
		next=next_node;
		next_node->prev=this;
	}

	void insert_prev(CNode<T>* prev_node){
		prev->next=prev_node;
		prev_node->prev=prev;
		prev=prev_node;
		prev_node->next=this;
	}
	void remove(){
		next->prev=prev;
		prev->next=next;
	}
};

template<int N,typename T>
class CIndexHashmap {
private:
	CNode<T> bucket[N+1];
	CNode<T> first[N+1];
	CNode<T> last[N+1];
	int data_count;
public:
	CIndexHashmap():data_count(0),table{}{
		for(int i=0;i<N+1;i++){
			//bucket[i]=CDoubleLinkElem<T>();
			first[i].prev=&first[i];
			first[i].next=&last[i];
			last[i].next=&last[i];
			last[i].prev=&first[i];
		}

	}
	int _hash(int i)const{
		return i%(N+1);
	}
	void emplace(int key,T item){
		int didx=data_count++;
		bucket[didx].data=item;
		int hash=_hash(key);
		last[hash].insert_prev(&bucket[didx]);
	}
	T operator[](int idx);
	//virtual ~CHashmap();
};
*/
#endif /* CHASHMAP_H_ */
