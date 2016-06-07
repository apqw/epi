#ifndef ATOMICS_H
#define ATOMICS_H
#include <atomic>
#include <array>
#include <utility>
#include <cassert>
#include <algorithm>

struct atomic_double {
private:
	std::atomic<double> v;
public:
	atomic_double() :v(0) {}
	atomic_double(double _v) :v(_v) {}
	atomic_double(atomic_double&& _v) :v(_v.v.load()) {}
	atomic_double(const atomic_double& _v) :v(_v.v.load()) {}
	atomic_double& operator=(const double c) {
		double expected = v.load(std::memory_order_relaxed);
		while (
			!v.compare_exchange_weak(expected, c));
		return *this;
	}

	atomic_double& operator=(const atomic_double& c) {
		double expected = v.load(std::memory_order_relaxed);
		double cv = (double)c;
		while (
			!v.compare_exchange_weak(expected, cv));
		return *this;
	}


	operator double() const {
		return v.load(std::memory_order_relaxed);
	}
	/*
	T& operator()() {
		return v.load(std::memory_order_relaxed);
	}
	*/

	atomic_double& operator+=(double c) {
		auto expected = v.load(std::memory_order_relaxed);
		while (
			!v.compare_exchange_weak(expected, expected + c));
		return *this;
	}
	atomic_double& operator*=(double c) {
		auto expected = v.load(std::memory_order_relaxed);
		while (
			!v.compare_exchange_weak(expected, expected * c));
		return *this;
	}
	atomic_double& operator-=(double c) {
		auto expected = v.load(std::memory_order_relaxed);
		while (
			!v.compare_exchange_weak(expected, expected - c));
		return *this;
	}
	atomic_double& operator/=(double c) {
		auto expected = v.load(std::memory_order_relaxed);
		while (
			!v.compare_exchange_weak(expected, expected / c));
		return *this;
	}

};
template<typename T, unsigned N>
class Lockfree_push_stack {
protected:
	T _data[N] = {};
	std::atomic<size_t> _next = 0;
	/*
	template<typename ITR,size_t... I>
	Lockfree_push_stack(ITR it,std::index_sequence<I...>):_data{*(it+I)...}{}
	*/
public:
    Lockfree_push_stack():_next(0) {}

	void push_back(T&& item) {
		_data[_next++] = item;
	}

	void push_back(const T& item) {
		_data[_next++] = item;
	}

	size_t push_back_with_index(T&& item) {
		size_t tmp = _next++;
		_data[tmp] = item;
		return tmp;
	}

	size_t push_back_with_index(const T& item) {
		size_t tmp = _next++;
		_data[tmp] = item;
		return tmp;
	}

	void clear() {
		_next = 0;
	}

	size_t size() const {
		return _next;
	}

	T& operator[](size_t idx) {
		return _data[idx];
	}

	const T& operator[](size_t idx)const {
		return _data[idx];
	}

	void unordered_remove_non_concurrent(size_t idx) {
		assert(_next > 0);
		_data[idx] = _data[--_next];
	}

	template<class Fn>
	void foreach(const Fn& lmbd) {
		for (size_t i = 0; i < size(); ++i) {
            lmbd((*this)[i]);
		}
	}

	template<class Fn>
	void foreach(const Fn& lmbd)const {
		for (size_t i = 0; i < size(); ++i) {
			lmbd((*this)[i]);
		}
	}

	template<class Fn>
	void foreach_with_index(const Fn& lmbd) {
		for (size_t i = 0; i < size(); ++i) {
			lmbd((*this)[i],i);
		}
	}

	void force_set_count(size_t c) {
		_next = c;
	}

	bool exist(const T& k) const{
		return std::find(_data, _data+_next, k) != _data + _next;
	}
};


#endif // ATOMICS_H
