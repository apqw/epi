#pragma once
#include "define.h"

namespace _TMP_impl_cpp14{
	// using aliases for cleaner syntax
	template<class T> using Invoke = typename T::type;

	template<unsigned...> struct seq{ using type = seq; };

	template<class S1, class S2> struct concat;

	template<unsigned... I1, unsigned... I2>
	struct concat<seq<I1...>, seq<I2...>>
		: seq<I1..., (sizeof...(I1)+I2)...>{};

	template<class S1, class S2>
	using Concat = Invoke<concat<S1, S2>>;

	template<unsigned N> struct gen_seq;
	template<unsigned N> using GenSeq = Invoke<gen_seq<N>>;

	template<unsigned N>
	struct gen_seq : Concat<GenSeq<N / 2>, GenSeq<N - N / 2>>{};

	template<> struct gen_seq<0> : seq<>{};
	template<> struct gen_seq<1> : seq<0>{};
}

namespace cpp14{
	template<unsigned N> using make_integer_sequence = _TMP_impl_cpp14::GenSeq<N>;
	template<unsigned... I> using integer_sequence = _TMP_impl_cpp14::seq<I...>;
}

namespace cpp11{
	template<typename T> struct type_identity{
		typedef T type;
	};
}

__device__ real p_diff_x(const real x1, const real x2);

__device__ real p_diff_x_old(const real x1, const real x2);

__device__ real p_diff_y(const real y1, const real y2);

__device__ real p_dist_sq(CellPos c1, CellPos c2);

template<unsigned X,unsigned Y,unsigned Z>
struct __midx{
	__host__ __device__ static int idx(int x, int y, int z){
		return x + X*y + X*Y*z;
	}
};

using f3di = __midx<NX + 1, NY + 1, NZ + 1>;

/*
__device__ short atomicAddShort(short* address, short val)

{

	unsigned int *base_address = (unsigned int *)((size_t)address & ~2);

	unsigned int long_val = ((size_t)address & 2) ? ((unsigned int)val << 16) : (unsigned short)val;



	unsigned int long_old = atomicAdd(base_address, long_val);

	if ((size_t)address & 2) {

		return (short)(long_old >> 16);

	}
	else {

		unsigned int overflow = ((long_old & 0xffff) + long_val) & 0xffff0000;

		if (overflow)

			atomicSub(base_address, overflow);

		return (short)(long_old & 0xffff);

	}

}
*/

//__device__ unsigned short atomicIncrementShort_no_overflow(unsigned short* address);

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true);
