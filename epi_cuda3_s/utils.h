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
