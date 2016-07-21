#pragma once
#include <array>
#include <cstdio>
#include <string>
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

__device__ inline real p_diff_x(const real x1, const real x2)
{

	const real diff = x1 - x2;
	return diff - LX*floorr(diff *INV_LX + CDEF(0.5));
}

__device__ inline real p_diff_x_old(const real x1, const real x2)
{
	const real diff = x1 - x2;
	if (diff > CDEF(0.5)*LX)return diff - LX;
	if (diff <= -CDEF(0.5)*LX)return diff + LX;
	return diff;
	
}

__device__ inline real p_diff_y(const real y1, const real y2)
{

	const real diff = y1 - y2;
	return diff - LY*floorr(diff *INV_LY + CDEF(0.5));
}

__device__ inline real p_dist_sq(CellPos c1, CellPos c2){
	const real diffx = p_diff_x(c1.x, c2.x);
	const real diffy = p_diff_y(c1.y, c2.y);
	const real diffz = c1.z - c2.z;
	return diffx*diffx + diffy*diffy + diffz*diffz;
}

template<unsigned X,unsigned Y,unsigned Z>
struct __midx{
	__host__ __device__ static int idx(int x, int y, int z){
		return x + X*y + X*Y*z;
	}
};

using f3di = __midx<NX + 1, NY + 1, NZ + 1>;