#ifndef ATOMICS_H
#define ATOMICS_H
#include <atomic>
#include <array>
#include <utility>
#include <cassert>
#include <algorithm>
#include <cstdlib>
#include <iostream>
template<typename FPTy>
struct atomic_fp {
private:
    std::atomic<FPTy> v;
public:
    atomic_fp() :v(0) {}
    atomic_fp(FPTy _v) :v(_v) {}
    atomic_fp(atomic_fp&& _v) :v(_v.v.load()) {}
    atomic_fp(const atomic_fp& _v) :v(_v.v.load()) {}
    atomic_fp& operator=(const FPTy c) {
        FPTy expected = v.load(std::memory_order_relaxed);
        while (
            !v.compare_exchange_weak(expected, c));
        return *this;
    }

    atomic_fp& operator=(const atomic_fp& c) {
        FPTy expected = v.load(std::memory_order_relaxed);
        FPTy cv = (FPTy)c;
        while (
            !v.compare_exchange_weak(expected, cv));
        return *this;
    }


    operator FPTy() const {
        return v.load(std::memory_order_relaxed);
    }

    atomic_fp& operator+=(FPTy c) {
        auto expected = v.load(std::memory_order_relaxed);
        while (
            !v.compare_exchange_weak(expected, expected + c));
        return *this;
    }
    atomic_fp& operator*=(FPTy c) {
        auto expected = v.load(std::memory_order_relaxed);
        while (
            !v.compare_exchange_weak(expected, expected * c));
        return *this;
    }
    atomic_fp& operator-=(FPTy c) {
        auto expected = v.load(std::memory_order_relaxed);
        while (
            !v.compare_exchange_weak(expected, expected - c));
        return *this;
    }
    atomic_fp& operator/=(FPTy c) {
        auto expected = v.load(std::memory_order_relaxed);
        while (
            !v.compare_exchange_weak(expected, expected / c));
        return *this;
    }

};

template<typename T>
class Lockfree_push_stack_dyn {
protected:
    std::vector<T> _data;
	std::atomic<size_t> _next = 0;
    //size_t allocated = 0;
	/*
	template<typename ITR,size_t... I>
	Lockfree_push_stack(ITR it,std::index_sequence<I...>):_data{*(it+I)...}{}
	*/
public:
    Lockfree_push_stack_dyn<T>(size_t N):_next(0),_data(N) {}

    Lockfree_push_stack_dyn<T>(Lockfree_push_stack_dyn<T>&& ot):_next(ot._next.load()),_data(std::move(ot._data)) {}
    Lockfree_push_stack_dyn<T>(const Lockfree_push_stack_dyn<T>& ot):_next(ot._next.load()),_data(ot._data) {}
    Lockfree_push_stack_dyn<T>& operator=(Lockfree_push_stack_dyn<T>&& ot){
    	_data=std::move(ot._data);
    	_next=ot._next.load();
    }

    Lockfree_push_stack_dyn<T>& operator=(const Lockfree_push_stack_dyn<T>& ot){
        	_data=ot._data;
        	_next=ot._next.load();
        }
    void test_realloc() {
        if (_data.size()==0||_next / (double)_data.size() > 0.8) {
            size_t origsize = _data.size();
            _data.resize(_data.size() == 0?2:_data.size() * 2);
            std::cout << "LFStack resized:"_s + std::to_string(origsize) + "->" + std::to_string(_data.size()) << std::endl;
        }
    }

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

    size_t max_size()const {
        return _data.size();
    }

    void resize(size_t sz) {
        if (sz <= _next)throw std::logic_error("Failed to resize the stack. Target size:"_s + std::to_string(sz) + " Current element num.:" + std::to_string(_next.load()));
        _data.resize(sz);
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
            lmbd(_data[i]);
		}
	}

	template<class Fn>
	void foreach(const Fn& lmbd)const {
		for (size_t i = 0; i < size(); ++i) {
			lmbd(_data[i]);
		}
	}

	template<class Fn>
	void foreach_with_index(const Fn& lmbd) {
		for (size_t i = 0; i < size(); ++i) {
			lmbd(_data[i],i);
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
