#pragma once
#include "define.h"
#include <initializer_list>
#include <cmath>
#include <cassert>
#include <array>
#include <algorithm>
#include <functional>
#include <numeric>
#include <atomic>

template<typename T,unsigned int N>
struct Vec {
	std::vector<T> value;

	Vec(const Vec<T,N>& v) :value(v.value) {}
	
	Vec(const std::vector<T>& v) :value(v) {}
	Vec(const std::array<T, N>& v){
		value.reserve(N);
		std::for_each(v.begin(), v.end(), [&](const T& a) {value.push_back(a); });
	}
	Vec(const T v[N]) {
		value.reserve(N);
		for (int i = 0; i < N; i++) {
			value.push_back(v[i]);
		}
	}
	
	Vec(const std::initializer_list<T>& v) :value(v){
	}

	explicit Vec(const T& _v) :value(std::vector<T>(N,_v)) {}

	template<typename U>
	Vec(const Vec<U,N>& _v) {
		value.reserve(N);
		for (int i = 0; i < N; i++) {
			value.push_back((T)(_v.value[i]));
		}
	}
	Vec():value(std::vector<T>(N)) {}
	template<typename L>
    void foreach(const L& lambda) {
		//std::for_each(value.begin(), value.end(), lambda); //gives const args
		for (auto it = value.begin(); it != value.end(); ++it) {
			lambda(*it);
		}
	}

	T& operator[](const std::size_t i) {
		return value[i];
	}
	
	Vec operator* (const Vec& v)const {
		std::vector<T> tmp; tmp.reserve(N);
		for (int i = 0; i < N; i++) {
			tmp.push_back(value[i] * v.value[i]);
		}
		return Vec<T, N>(tmp);
	}

	template<typename U>
	Vec operator* (const U& v)const {
		return *this*((Vec)v);
	}

	template<typename U>
	friend Vec operator* (const U& c, const Vec& v) {
		return ((Vec)c)*v;
	}

	Vec operator+ (const Vec& v)const {
		std::vector<T> tmp; tmp.reserve(N);
		for (int i = 0; i < N; i++) {
			tmp.push_back(value[i] + v.value[i]);
		}
		return Vec<T, N>(tmp);
	}

	template<typename U>
	Vec operator+ (const U& v)const {
		return *this+((Vec)v);
	}

	template<typename U>
	friend Vec operator+ (const U& c, const Vec& v) {
		return ((Vec)c)+v;
	}

	Vec operator-(const Vec& v) {
		std::vector<T> tmp; tmp.reserve(N);
		for (int i = 0; i < N; i++) {
			tmp.push_back(value[i] - v.value[i]);
		}
		return Vec<T, N>(tmp);
	}

	template<typename U>
	Vec operator-(const U& v)const {
		return *this - (Vec)v;
	}

	Vec operator/(const Vec& v)const {
		std::vector<T> tmp; tmp.reserve(N);
		for (int i = 0; i < N; i++) {
			tmp.push_back(value[i] / v.value[i]);
		}
		return Vec<T, N>(tmp);
	}

	template<typename U>
	Vec operator/(const U& v)const {
		return *this / (Vec)v;
	}

	/*
	template<typename U>
	friend bool operator== (const U& c, const Vec& v) {
		return ((Vec)c).value == v.value;
	}
	
	template<typename U>
	bool operator== (const U& v)const {
		return value == ((Vec)v).value;//std::array compares elements
	}
	*/
	bool operator== (const Vec& v)const {
		return value == v.value;//std::array compares elements
	}

	Vec operator-()const {
		std::vector<T> tmp; tmp.reserve(N);
		for (int i = 0; i < N; i++) {
			tmp.push_back(-value[i]);
		}
		return Vec<T, N>(tmp);
	}


	
	Vec<T, N>& operator=(const Vec<T, N>& ot) {
		//self check?
		value = ot.value;
		return *this;
	}
	
	template<typename U>
	Vec<T, N>& operator=(const Vec<U, N>& ot) {
		//self check?
		for (int i = 0; i < N; i++) {
			value[i] = (T)(ot.value[i]);
		}
		return *this;
	}

	Vec<T, N>& operator+=(const Vec<T, N>& ot) {
		for (int i = 0; i < N; i++) {
			value[i] += ot.value[i];
		}
		return *this;
	}

	template<typename U>
	Vec<T, N>& operator+=(const Vec<U, N>& ot) {
		for (int i = 0; i < N; i++) {
			value[i] += (T)(ot.value[i]);
		}
		return *this;
	}

	Vec<T, N>& operator-=(const Vec<T, N>& ot) {
		for (int i = 0; i < N; i++) {
			value[i] -= ot.value[i];
		}
		return *this;
	}

	template<typename U>
	Vec<T, N>& operator-=(const Vec<U, N>& ot) {
		for (int i = 0; i < N; i++) {
			value[i] -= (T)(ot.value[i]);
		}
		return *this;
	}

	Vec<T, N>& operator*=(const Vec<T, N>& ot) {
		for (int i = 0; i < N; i++) {
			value[i] *= ot.value[i];
		}
		return *this;
	}
	template<typename U>
	Vec<T, N>& operator*=(const Vec<U, N>& ot) {
		for (int i = 0; i < N; i++) {
			value[i] *= (T)(ot.value[i]);
		}
		return *this;
	}

	Vec<T, N>& operator/=(const Vec<T, N>& ot) {
		for (int i = 0; i < N; i++) {
			value[i] /= ot.value[i];
		}
		return *this;
	}
	template<typename U>
	Vec<T, N>& operator/=(const Vec<U, N>& ot) {
		for (int i = 0; i < N; i++) {
			value[i] /= (T)(ot.value[i]);
		}
		return *this;
	}

	T sum() {
		return std::accumulate(value.begin(), value.end(), 0.0,std::plus<T>());
	}

	T normSq() {
		auto tmp = (*this)*(*this);
		return tmp.sum();
	}

};

template<typename T>
struct Vec<T, 3> {
	std::vector<T> value;

	Vec(const Vec<T, 3>& v) :value(v.value) {}

	Vec(const std::vector<T>& v) :value(v) {}
	Vec(const std::array<T, 3>& v) :value{ v[0],v[1],v[2] }{}
	Vec(const T v[3]):value{v[0],v[1],v[2]}{}

	Vec(const std::initializer_list<T>& v) :value(v) {
	}

	explicit Vec(const T& _v) :value(std::vector<T>(3, _v)) {}

	template<typename U>
	Vec(const Vec<U, 3>& _v) {
		value.resize(3);
		value[0] = (T)(_v.value[0]);
		value[1] = (T)(_v.value[1]);
		value[2] = (T)(_v.value[2]);
	}
	Vec() :value(std::vector<T>(3)) {}
	template<typename L>
    void foreach(const L& lambda) {
		//std::for_each(value.begin(), value.end(), lambda); //gives const args
		lambda(value[0]);
		lambda(value[1]);
		lambda(value[2]);
	}

    const T& operator[](const std::size_t i)const {
		return value[i];
	}

    T& operator[](const std::size_t i) {
        return value[i];
    }


	Vec operator* (const Vec& v)const {
		return Vec<T, 3>({ value[0] * v.value[0], value[1] * v.value[1],value[2] * v.value[2] });
	}

	template<typename U>
	Vec operator* (const U& v)const {
		return Vec<T, 3>({ value[0] * (T)v, value[1] * (T)v,value[2] * (T)v });
	}

	template<typename U>
	friend Vec operator* (const U& c, const Vec& v) {
		return Vec<T, 3>({ v.value[0] * (T)c, v.value[1] * (T)c,v.value[2] * (T)c });
	}

	Vec operator+ (const Vec& v)const {
		return Vec<T, 3>({ value[0] + v.value[0], value[1] + v.value[1],value[2] + v.value[2] });
	}

	template<typename U>
	Vec operator+ (const U& v)const {
		return Vec<T, 3>({ value[0] + (T)v, value[1] +(T)v,value[2] + (T)v });
	}

	template<typename U>
	friend Vec operator+ (const U& c, const Vec& v) {
		return Vec<T, 3>({ v.value[0] + (T)c, v.value[1] + (T)c,v.value[2] + (T)c });
	}

	Vec operator-(const Vec& v) {
		return Vec<T, 3>({value[0]-v.value[0],value[1]-v.value[1],value[2]-v.value[2]});
	}

	template<typename U>
	Vec operator-(const U& v)const {
		return Vec<T, 3>({ value[0] - (T)v,value[1] - (T)v,value[2] - (T)v });
	}

	Vec operator/(const Vec& v) {
		return Vec<T, 3>({ value[0] / v.value[0],value[1] / v.value[1],value[2] / v.value[2] });
	}

	template<typename U>
	Vec operator/(const U& v)const {
		return Vec<T, 3>({ value[0] / (T)v,value[1] / (T)v,value[2]/ (T)v });
	}


	template<typename U>
	friend bool operator== (const U& c, const Vec& v) {
		return v.value[0] == (T)c && v.value[1] == (T)c &&v.value[2] == (T)c;
	}

	template<typename U>
	bool operator== (const std::array<U, 3>& v)const {
		return value[0] == (T)v[0] && value[1] == (T)v[1] && value[2] == (T)v[2];
	}
	/*
	template<typename U>
	bool operator== (const U& v)const {
		return value[0] == (T)v && value[1] == (T)v &&value[2] == (T)v;
	}

	*/

	template<typename U>
	bool operator== (const std::vector<U>& v)const {
		return value[0] == (T)v[0] && value[1] == (T)v[1] && value[2] == (T)v[2];
	}

	bool operator== (const Vec& v)const {
		return value == v.value;//std::array compares elements
	}

	Vec operator-()const {
		return Vec<T, 3>({ -value[0],-value[1],-value[2] });
	}



	Vec<T, 3>& operator=(const Vec<T, 3>& ot) {
		//self check?
		value = ot.value;
		return *this;
	}

	template<typename U>
	Vec<T, 3>& operator=(const Vec<U, 3>& ot) {
		//self check?
		value[0] = (T)(ot.value[0]);
		value[1] = (T)(ot.value[1]);
		value[2] = (T)(ot.value[2]);
		return *this;
	}

	Vec<T, 3>& operator+=(const Vec<T, 3>& ot) {
		value[0] += ot.value[0];
		value[1] += ot.value[1];
		value[2] += ot.value[2];
		return *this;
	}

	template<typename U>
	Vec<T, 3>& operator+=(const Vec<U, 3>& ot) {
		value[0] += (T)(ot.value[0]);
		value[1] += (T)(ot.value[1]);
		value[2] += (T)(ot.value[2]);
		return *this;
	}

	Vec<T, 3>& operator-=(const Vec<T, 3>& ot) {
		value[0] -= ot.value[0];
		value[1] -= ot.value[1];
		value[2] -= ot.value[2];
		return *this;
	}

	template<typename U>
	Vec<T, 3>& operator-=(const Vec<U, 3>& ot) {
		value[0] -= (T)(ot.value[0]);
		value[1] -= (T)(ot.value[1]);
		value[2] -= (T)(ot.value[2]);
		return *this;
	}

	Vec<T, 3>& operator*=(const Vec<T, 3>& ot) {
		value[0] *= ot.value[0];
		value[1] *= ot.value[1];
		value[2] *= ot.value[2];
		return *this;
	}
	template<typename U>
	Vec<T, 3>& operator*=(const Vec<U, 3>& ot) {
		value[0] *= (T)(ot.value[0]);
		value[1] *= (T)(ot.value[1]);
		value[2] *= (T)(ot.value[2]);
		return *this;
	}

	Vec<T, 3>& operator/=(const Vec<T, 3>& ot) {
		value[0] /= ot.value[0];
		value[1] /= ot.value[1];
		value[2] /= ot.value[2];
		return *this;
	}
	template<typename U>
	Vec<T, 3>& operator/=(const Vec<U, 3>& ot) {
		value[0] /= (T)(ot.value[0]);
		value[1] /= (T)(ot.value[1]);
		value[2] /= (T)(ot.value[2]);
		return *this;
	}

	T sum() {
		return value[0] + value[1] + value[2];
	}

	T normSq() {
		return value[0] * value[0] + value[1] * value[1] + value[2] * value[2];
	}

};
template<typename T>
using Vec3 = Vec<T, 3>;
template<typename T>
class lazy_update {
protected:
	T value;
	std::atomic<T> next;
public:
	
	lazy_update(const lazy_update& init) :value(init.value), next(init.next.load()) {}
	lazy_update(lazy_update&& init) :value(init.value), next(init.next.load()) {}
	lazy_update(T& init) :value(init), next(init) {}
	virtual void update() {
		value = next.load(std::memory_order_relaxed);
	};

	const T& operator()() const {
		return value;
	}

	explicit operator T()const {
		return value;
	}

	lazy_update<T>& operator=(const T& v) {
		next = v;
		return *this;
	}
};

template<typename T>
struct DV :public lazy_update<T>{

public:
    DV(T init) :lazy_update<T>(init) {}


	

	//operator = is not allowed.
	DV<T>& operator=(const DV& v) {
		assert(!"operator = is not allowed.");
		return *this;
	}

	DV<T>& operator=(const T& v) = delete;

	void copy_from(const DV<T>& v) {
        this->value = v.value;
        this->next = v.next;
	}
	DV operator-() {
        return DV<T>(-this->value);
	}
	friend DV<T> operator* (const DV<T>& c, const DV<T>& v) {
		return DV<T>(c()*v());
	}
	friend DV<T> operator+ (const DV<T>& c, const DV<T>& v) {
		return DV<T>(c() + v());
	}
	friend DV operator- (const DV<T>& c, const DV<T>& v) {
		return DV<T>(c() - v());
	}
	friend DV<T> operator/ (const DV<T>& c, const DV<T>& v) {
		return DV<T>(c() / v());
	}
	friend bool operator== (const DV<T>& c, const DV<T>& v) {
		return c() == v();
	}
	friend bool operator!= (const DV<T>& c, const DV<T>& v) {
		return c() != v();
	}

	friend bool operator< (const DV<T>& c, const DV<T>& v) {
		return c() < v();
	}

	friend bool operator> (const DV<T>& c, const DV<T>& v) {
		return c() > v();
	}

	friend bool operator<= (const DV<T>& c, const DV<T>& v) {
		return c() <= v();
	}

	friend bool operator>= (const DV<T>& c, const DV<T>& v) {
		return c() >= v();
	}
	/*
	void add_diff(const T& v) {
	next += v;
	}
	*/

	DV<T>& operator+=(const T& c) {
        auto expected = this->next.load(std::memory_order_relaxed);
		while (
            !this->next.compare_exchange_weak(expected, expected + c));
		return *this;
	}

	DV<T>& operator+=(const DV<T>& c) {
        auto expected = this->next.load(std::memory_order_relaxed);
		while (
            !this->next.compare_exchange_weak(expected, expected + c()));
		return *this;
	}

	DV<T>& operator-=(const T& c) {
        auto expected = this->next.load(std::memory_order_relaxed);
		while (
            !this->next.compare_exchange_weak(expected, expected - c));
		return *this;
	}

	DV<T>& operator-=(const DV<T>& c) {
		//assert(!"operator -= is not allowed.");
        auto expected = this->next.load(std::memory_order_relaxed);
		while (
            !this->next.compare_exchange_weak(expected, expected - c()));
		return *this;
	}


	DV<T>& operator*=(const T& c) {
        auto expected = this->next.load(std::memory_order_relaxed);
		while (
            !this->next.compare_exchange_weak(expected, expected * c));
		return *this;
	}

	DV<T>& operator*=(const DV<T>& c) {
		//assert(!"operator *= is not allowed.");
        auto expected = this->next.load(std::memory_order_relaxed);
		while (
            !this->next.compare_exchange_weak(expected, expected * c()));
		return *this;
	}

	DV<T>& operator/=(const T& c) {
        auto expected = this->next.load(std::memory_order_relaxed);
		while (
            !this->next.compare_exchange_weak(expected, expected / c));
		return *this;
	}

	DV<T>& operator/=(const DV<T>& c) {
		//assert(!"operator /= is not allowed.");
        auto expected = this->next.load(std::memory_order_relaxed);
		while (
            !this->next.compare_exchange_weak(expected, expected / c()));
		return *this;
	}

	void force_set_next_value(const T& v) {
        this->next = v;
	}
};

template<typename T>
struct CAS {
	std::atomic<T> v;
    CAS(T _v) :v(_v) {}
    CAS(CAS&& _v):v(_v.v.load()){}
	CAS& operator=(T c) {
		auto expected = v.load(std::memory_order_relaxed);
		while (
			!v.compare_exchange_weak(expected,c));
		return *this;
	}
	operator T() const{
		return v.load(std::memory_order_relaxed);
	}
	T& operator()() {
		return v.load(std::memory_order_relaxed);
	}
	CAS& operator+=(T c) {
		auto expected = v.load(std::memory_order_relaxed);
		while (
			!v.compare_exchange_weak(expected, expected+c));
		return *this;
	}
	CAS& operator*=(T c) {
		auto expected = v.load(std::memory_order_relaxed);
		while (
			!v.compare_exchange_weak(expected, expected * c));
		return *this;
	}
	CAS& operator-=(T c) {
		auto expected = v.load(std::memory_order_relaxed);
		while (
			!v.compare_exchange_weak(expected, expected - c));
		return *this;
	}
	CAS& operator/=(T c) {
		auto expected = v.load(std::memory_order_relaxed);
		while (
			!v.compare_exchange_weak(expected, expected / c));
		return *this;
	}
};

using real = DV<CAS<double>>;
using real_raw = double;
/*
template<>
struct CAS<double> {
	std::atomic<double> v;
	CAS& operator=(double c) {
		double expected = v.load(std::memory_order_relaxed);
		while (
			!std::atomic_compare_exchange_weak
			(&v, &expected,expected)
			);
		return *this;
	}
	CAS& operator+=(double c) {
		return operator=(v + c);
	}
	CAS& operator*=(double c) {
		return operator=(v * c);
	}
	CAS& operator-=(double c) {
		return operator=(v - c);
	}
	CAS& operator/=(double c) {
		return operator=(v / c);
	}
};
*/
