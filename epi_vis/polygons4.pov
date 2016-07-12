#version 3.6;
#include "colors.inc"
#include "metals.inc"
#include "textures.inc"

global_settings {
	max_trace_level 64
}

camera {
	location <-120,-240,40>
	sky z
	right 0.24*x*image_width/image_height
	up 0.24*z
	look_at <25,25,25>
}

background{rgb 0}

light_source{<-8,-20,30> color rgb <0.77,0.75,0.75>}
light_source{<25,-12,12> color rgb <0.43,0.45,0.45>}

#declare r=0.1;
#declare rr=0.08;
#declare mth=1.0;
#declare mrad=1.8;
#declare mstr=2.0;
#declare mscl=0.9;
#declare mztr=-1.0;
#declare mztr2=-1.0;
#declare mupz=3.0;
#declare crad=0.8;

#declare f1=finish{reflection 0.15 specular 0.3 ambient 0.42}

#declare t1=texture{pigment{rgbft <1.0,0.6,0.8,0,0.4>} finish{f1}} //dead
#declare t4=texture{pigment{rgbft <0.2,0.0,0.2,0,0.3>} finish{f1}} //fix
#declare t7=texture{pigment{rgbft <0.30,0.0,0.30,0,0.3>} finish{f1}} //musume
#declare t0=texture{pigment{rgbft <0.5,0.122,0.5,0,0.4>} finish{f1}} //al
#declare c1=texture{pigment{rgbft <0.9,0.8,0.9,0.7>} finish{f1}}
#declare r1=texture{pigment{rgb <0.35,0.35,0.35>} finish{f1}}
#declare cell1=texture{pigment{rgb <0.9,0.9,0.9>} finish{f1}}
#declare m1=texture{pigment{rgbft <0.35,0.35,0.35,0.15>} finish{f1}}
union{
#include "outp.pov"
}
//union{
//#include "polygons_d.pov"
//	texture{T_Silver_4B}
//}
