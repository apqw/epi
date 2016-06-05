#pragma once
#include <array>
#include <vector>
#include <utility>

#define DEFUFUNC(F,fname) \
template<typename T,class = typename VecValidType1<T>>\
VUnary<T, F, T::size> fname(const T& l){\
	return VUnary<T, F, T::size>(l);\
}

#define DEFUOP(F,op) DEFUFUNC(F,operator##op)

#define DEFBFUNC(F,fname)\
template<typename T,typename U,class = typename VecValidType2<T,U>>\
VExpr<T, F, U, T::size> fname(const T& l, const U& r) {\
	return VExpr<T, F, U, T::size>(l, r);\
}

#define DEFBOP(F,op) DEFBFUNC(F,operator##op)


namespace aml {



	template <std::size_t N>
	class VecBase {
	public:
		static constexpr std::size_t size = N;
	};

	template<class T>
	using VecValidType1 = typename std::enable_if<std::is_base_of<VecBase<T::size>, T>::value>::type;

	template<class T, class U>
	using VecValidType2 = typename std::enable_if<std::is_base_of<VecBase<T::size>, T>::value && std::is_base_of<VecBase<U::size>, U>::value && T::size == U::size>::type;

	
	
	template <class L, template<class,class> class Op, class R,std::size_t N>
	class VExpr:public VecBase<N> {
		const L& l_;
		const R& r_;
	public:
		VExpr(const L& l, const R& r)
			: l_(l), r_(r) {}

		auto operator[](std::size_t i) const
		{
			return Op<decltype(l_[i]),decltype(r_[i])>::Apply(l_[i], r_[i]);
		}


	};

	template <class T, template<class> class UOp, std::size_t N>
	class VUnary :public VecBase<N> {
		const T& v_;
	public:
		VUnary(const T& l)
			: v_(l) {}

		auto operator[](std::size_t i) const
		{
			return UOp<decltype(v_[i])>::Apply(v_[i]);
		}


	};

	template<typename T,std::size_t N>
	class Vec:public VecBase<N>
	{
	private:
		std::array<T, N> _value;

		template<typename itr_begin,std::size_t... I>
		explicit Vec(const itr_begin& arr_itr_begin, std::index_sequence<I...>) :_value{ *(arr_itr_begin + I)... } {};

		template <class L, template<class, class> class Op, class R, std::size_t... I>
		explicit Vec(const VExpr<L, Op, R,N>& vex, std::index_sequence<I...>) :_value{vex[I]...} {};
	public:
		template <class L, template<class, class> class Op, class R>
		Vec(const VExpr<L, Op, R,N>& vex) :Vec(vex, std::make_index_sequence<N>()) {};

		Vec(const std::array<T, N>& arr) :_value(arr) {};

		Vec(const std::initializer_list<T>& iarr) :Vec(iarr.begin(), std::make_index_sequence<N>()) {};

		Vec(const T(&arr)[N]) :Vec(arr, std::make_index_sequence<N>()) {};

		Vec(const std::vector<T>& iarr) :Vec(iarr.cbegin(), std::make_index_sequence<N>()) {};


		T& operator[](std::size_t idx) {
			return _value[idx];
		}

		const T& operator[](std::size_t idx) const {
			return _value[idx];
		}

		template<class E>
		Vec& operator=(const E& r) {
			for (std::size_t i = 0; i < N;i++) {
				_value[i] = r[i];
			}
			return *this;
		}

		template<class E>
		Vec& operator+=(const E& r) {
			*this = *this + r;
			return *this;
		}

		~Vec() {};
	};



	namespace Ops {

		template<typename T, typename U>
		struct _Plus {
			static auto Apply(const T& l, const U& r) {
				return l + r;
			}
		};

		template<typename T, typename U>
		struct _Minus {
			static auto Apply(const T& l, const U& r) {
				return l - r;
			}
		};

		template<typename T, typename U>
		struct _Mul {
			static auto Apply(const T& l, const U& r) {
				return l * r;
			}
		};

		template<typename T, typename U>
		struct _Div {
			static auto Apply(const T& l, const U& r) {
				return l / r;
			}
		};

		template<typename T, typename U>
		struct _Mod {
			static auto Apply(const T& l, const U& r) {
				return l % r;
			}
		};

		/*
			Unary
		*/


		template<typename T>
		struct _Negate {
			static auto Apply(const T& l) {
				return -l;
			}
		};

		template<class CallableClass,typename T>
		struct _UFunc_Impl {
		private:
			static CallableClass f;
		public:
			static auto Apply(const T& l) {
				return f(l);
			}
		};

		template<class CallableClass>
		struct _UFunc {
			template<class T>
			using value = typename _UFunc::template _UFunc_Impl<CallableClass, T>;
		};

	}


	DEFBOP(Ops::_Plus, +);

	DEFBOP(Ops::_Minus, -);

	DEFBOP(Ops::_Mul, *);

	DEFBOP(Ops::_Div, /);

	DEFUOP(Ops::_Negate, -);

	template <typename T, typename U, class = typename VecValidType2<T, U>>
	bool operator==(const T& l, const U& r) {
		bool ret = true;
		for (std::size_t i = 0; i < T::size; i++) {
			ret = ret && (l[i] == r[i]);
		}
		return ret;
	}

	template <typename T, typename U, class = typename VecValidType2<T, U>>
	bool operator!=(const T& l, const U& r) {
		return !(l == r);
	}
	
	template<class CallableClass,typename T,class = typename VecValidType1<T>>
	VUnary<T, (Ops::_UFunc<CallableClass>::value), T::size> Apply(const T& l) {
			return VUnary<T, Ops::_UFunc<CallableClass>::value<T>, T::size>(l);
	}
	
}