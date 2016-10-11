#pragma once
#ifndef UTILS_H_
#define UTILS_H_
#include <string>
#include <memory>
#include <type_traits>
#include "define.h"
#define DIST_SQ(x,y,z) ((x)*(x)+(y)*(y)+(z)*(z))

std::string operator"" _s(const char* ,std::size_t);

real p_diff_x(const real x1, const real x2);
real p_diff_y(const real y1, const real y2);
real p_dist_sq(const real x1, const real y1, const real z1, const real x2, const real y2, const real z2);
real min0(const real a);
std::string state_to_str(CELL_STATE);
template<typename T, typename ...Args>
std::unique_ptr<T> make_unique_c11( Args&& ...args )
{
    return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
}


extern int* __lat_x;
extern int* __lat_y;
extern int* __lat_z;
extern int* per_x_prev_idx;
extern int* per_x_next_idx;
extern int* per_y_prev_idx;
extern int* per_y_next_idx;
extern int* per_z_prev_idx;
extern int* per_z_next_idx;

void init_precalc_lat();
inline int* precalc_lat_x(){
	return __lat_x;
}
inline int* precalc_lat_y(){
	return __lat_y;
}
inline int* precalc_lat_z(){
	return __lat_z;
}

void init_precalc_per();

inline int* precalc_per_next_x(){
	return per_x_next_idx;
}

inline int* precalc_per_prev_x(){
	return per_x_prev_idx;
}
inline int*precalc_per_next_y(){
	return per_y_next_idx;
}

inline int*precalc_per_prev_y(){
	return per_y_prev_idx;
}

inline int*precalc_per_next_z(){
	return per_z_next_idx;
}

inline int*precalc_per_prev_z(){
	return per_z_prev_idx;
}

template<typename T, typename... Rest>
struct is_same_multiple : std::false_type {};

template<typename T, typename First>
struct is_same_multiple<T, First> : std::is_same<T, First> {};

template<typename T, typename First, typename... Rest>
struct is_same_multiple<T, First, Rest...>
    : std::integral_constant<bool, std::is_same<T, First>::value && is_same_multiple<T, Rest...>::value>
{};

template<typename T, typename... U>
struct single_arg_and_eq_ref {
    static constexpr bool value = sizeof...(U) == 1
        &&
        is_same_multiple<
        T,
        typename std::remove_reference<U>::type...
        >::value;
};
#endif
