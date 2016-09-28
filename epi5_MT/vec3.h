#ifndef VEC3_H
#define VEC3_H
#include "utils.h"
#include <cmath>
struct vec3{
    double x,y,z;
    vec3():x(0),y(0),z(0){}
    vec3(double _x,double _y,double _z):x(_x),y(_y),z(_z){}
    vec3 operator+(const vec3& v)const{
        return vec3(x+v.x,y+v.y,z+v.z);
    }
    vec3 operator-(const vec3& v)const{
        return vec3(x-v.x,y-v.y,z-v.z);
    }

    vec3 operator*(double c)const{
        return vec3(c*x,c*y,c*z);
    }

    vec3 operator/(double c)const{
        return vec3(x/c,y/c,z/c);
    }
    static double dot(const vec3& v1,const vec3& v2){
        return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
    }
    static vec3 cross(const vec3& v1,const vec3& v2){
        return vec3(
                v1.y*v2.z-v1.z*v2.y,
                v1.z*v2.x-v1.x*v2.z,
                v1.x*v2.y-v1.y*v2.x
                    );
    }
    double normSq()const{
        return x*x+y*y+z*z;
    }
    double norm()const{
        return sqrt(normSq());
    }
    vec3& normalize(){
        double nm=norm();
        x/=nm;y/=nm;z/=nm;
        return *this;
    }
    vec3& operator+=(const vec3& v){
        x+=v.x;y+=v.y;z+=v.z;
        return *this;
    }
};
inline vec3 operator*(double c,const vec3& v){
    return vec3(c*v.x,c*v.y,c*v.z);
}
inline vec3 vpsub(const vec3& v1,const vec3& v2){
    return vec3(p_diff_x(v1.x,v2.x),p_diff_y(v1.y,v2.y),v1.z-v2.z);
}

#endif // VEC3_H
