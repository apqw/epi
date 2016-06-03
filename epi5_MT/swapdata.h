#ifndef SWAPDATA_H
#define SWAPDATA_H

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

template<typename T, T& sw, int subidx>
struct SwapArrAccessor2 {
private:
	size_t idx;
public:
	const auto& operator()() {
		return sw.first()[idx][subidx];
	}
	template<typename U>
	SwapArrAccessor2& operator+=(const U& v) {
		sw.second()[idx][subidx] += v;
		return *this;
	}

	template<typename U>
	void _set(U& v) {
		sw.first()[idx][subidx] = v;
	}

	SwapArrAccessor2() {}

	void _migrate(const size_t dest_idx) {
		sw.first()[dest_idx][subidx] = sw.first()[idx][subidx];
		sw.second()[dest_idx][subidx] = sw.second()[idx][subidx];
		idx = dest_idx;
	}
	void init(const size_t _idx) {
		idx = _idx;
	}
};

template<typename T, T& sw>
struct SwapArrAccessor1 {
	size_t idx;
	const auto& operator()() {
		return sw.first()[idx];
	}
	template<typename U>
	SwapArrAccessor1& operator+=(const U& v) {
		sw.second()[idx] += v;
		return *this;
	}
	template<typename U>
	void _set(U& v) {
		sw.first()[idx] = v;
	}

	SwapArrAccessor1() {}
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
