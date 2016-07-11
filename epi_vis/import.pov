#version 3.6;

// Right-handed coordinate system in which the z-axis points upwards
camera {
	location <25,-250,25>
	sky z
	right 0.24*x*image_width/image_height
	up 0.24*z
	look_at <25,50,25>
}

// White background
background{rgb 1}

// Two lights with slightly different colors
light_source{<-8,-20,30> color rgb <0.77,0.75,0.75>}
light_source{<25,-12,12> color rgb <0.38,0.40,0.40>}

// Radius of the Voronoi cell network
#declare r=0.08;

// Radius of the particles
#declare s=1;

// Particles
//union{
//#include "cell_voro_test_p.pov"
//	pigment{rgb 0.95} finish{reflection 0.1 specular 0.3 ambient 0.42}
//}

// Voronoi cells
union{
#include "cell_voro_test.pov"
	pigment{rgbf <1,0.4,0.45,1>} finish{specular 0.5 ambient 0.42}
}
