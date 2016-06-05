#ifndef SWAPDATA_H
#define SWAPDATA_H
#include <unordered_map>
template<typename T>
class SwapData
{
    T data1_internal;
    T data2_internal;
    T* data1;
    T* data2;
public:
    SwapData(){
        data1=&data1_internal;
        data2=&data2_internal;
    }
    SwapData(T&& data):data1_internal(data),data2_internal(data),SwapData(){}
    SwapData(const T& data):data1_internal(data),data2_internal(data),SwapData(){}

    T& first(){
return *data1;
    }
    T& second(){
        return *data2;
    }
    void swap(){
        T* tmp=data2;
        data2=data1;
        data1=tmp;
    }
};
template<typename Key,typename Value,size_t N>
struct SwapUMapArrAccessor {
private:
	size_t idx;
	SwapData<std::unordered_map<Key,Value>[N]>& sw;
public:

	void _emplace(const Key& k,const Value& v) {
		sw.first()[idx].emplace(k, v);
		sw.second()[idx].emplace(k, v);
	}

	void _try_emplace(const Key& k, const Value& v) {
		sw.first()[idx].try_emplace(k, v);
		sw.second()[idx].try_emplace(k, v);
	}
	
	template<class Fn>
	void erase(const Fn& cond) {
		auto& data1 = sw.first()[idx];
		auto& data2 = sw.second()[idx];

		for (auto it1 = data1.begin(), it2 = data2.begin(); it1 != data1.end();) {
			if (cond((const std::unordered_map<Key, Value>::const_iterator)it1)) {
				it1 = data1.erase(it1);
				it2 = data2.erase(it2);
			}
			else {
				++it1;
				++it2;
			}
		}
	}

	bool exist(const Key& k) {
		return sw.first()[idx].find(k) != sw.first()[idx].end();
	}

	

	void _set(const Key& k, const Value& v) {
		sw.first()[idx].at(k)=v;
		sw.second()[idx].at(k) = v;
	}

	auto& at_w(const Key& k) {
		return sw.second()[idx].at(k);
	}


	const auto& operator()(const Key& k)const {
		return sw.first()[idx].at(k);
	}


	SwapUMapArrAccessor(SwapData<std::unordered_map<Key, Value>[N]>& _sw) :sw(_sw) {}

	void _migrate(const size_t dest_idx) {
		sw.first()[dest_idx].swap(sw.first()[idx]);
		sw.second()[dest_idx].swap(sw.second()[idx]);
		idx = dest_idx;
	}
	void init(const size_t _idx) {
		idx = _idx;
	}
};

template<typename T,int subidx>
struct SwapArrAccessor2 {
private:
	size_t idx;
	SwapData<T>& sw;
public:
	const auto& operator()()const {
		return sw.first()[idx][subidx];
	}
	template<typename U>
	SwapArrAccessor2& operator+=(const U& v) {
		sw.second()[idx][subidx] += v;
		return *this;
	}

	template<typename U>
	SwapArrAccessor2& operator-=(const U& v) {
		sw.second()[idx][subidx] -= v;
		return *this;
	}

	template<typename U>
	void _set(const U& v) {
		sw.first()[idx][subidx] = v;
		sw.second()[idx][subidx] = v;
	}

	SwapArrAccessor2(SwapData<T>& _sw):sw(_sw) {}

	void _migrate(const size_t dest_idx) {
		sw.first()[dest_idx][subidx] = sw.first()[idx][subidx];
		sw.second()[dest_idx][subidx] = sw.second()[idx][subidx];
		idx = dest_idx;
	}
	void init(const size_t _idx) {
		idx = _idx;
	}
};

template<typename T>
struct SwapArrAccessor1 {
private:
	size_t idx;
	SwapData<T>& sw;
public:
	const auto& operator()()const {
		return sw.first()[idx];
	}

	/*
	template<typename U>
	SwapArrAccessor1& operator+=(const U& v) {
		sw.second()[idx] += v;
		return *this;
	}

	template<typename U>
	SwapArrAccessor1& operator-=(const U& v) {
		sw.second()[idx] -= v;
		return *this;
	}
	*/

	template<typename U>
	SwapArrAccessor1& operator=(const U& v) {
		sw.second()[idx] = v;
		return *this;
	}
	template<typename U>
	void _set(const U& v) {
		sw.first()[idx] = v;
		sw.second()[idx] = v;
	}

	SwapArrAccessor1(SwapData<T>& _sw):sw(_sw) {}
	void _migrate(const size_t dest_idx) {
		sw.first()[dest_idx] = sw.first()[idx];
		sw.second()[dest_idx] = sw.second()[idx];
		idx = dest_idx;
	}
	void init(const size_t _idx) {
		idx = _idx;
	}
};

#endif // SWAPDATA_H
