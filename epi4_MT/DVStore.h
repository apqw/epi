#pragma once
#include <utility>
template<typename T>
class DVStore
{
private:
	T _internal1;
	T _internal2;

	T* _old;
	T* _new;
public:
	
	DVStore() {
		_old = &_internal1;
		_new = &_internal2;
	}
	DVStore(T&& init):_internal1(init),_internal2(init):DVStore() {
		
	}
	//~DVStore();

	const T& old() const {
		return *_old;
	}
	T& _raw_old() {
		return *_old;
	}
	T& operator()() {
		return *_new;
	}
	void swap() {
		std::swap(_old, _new);
	}
};

