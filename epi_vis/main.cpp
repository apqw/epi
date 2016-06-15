#include <iostream>
#include <fstream>
#include <string>
#include <voro/voro++.hh>
using namespace std;

static constexpr double LX=50.0;
static constexpr double LY=50.0;
static constexpr double LZ=100.0;

static constexpr int comp_nx=8,comp_ny=8,comp_nz=8;

enum CELL_STATE:int {
    ALIVE = 0,
    DEAD = 1,
    DISA = 2,
    UNUSED = 3, //
    FIX = 4,
    BLANK = 5,
    DER = 6,
    MUSUME = 7,
    AIR = 8,
    MEMB = 9
};

bool count_cond(int state) {
	return state == FIX || state == MUSUME || state == ALIVE || state == DEAD;
}

void draw_polygon(FILE *fp, vector<int> &f_vert, vector<double> &v, int j) {
	static char s[6][128];
	int k, l, n = f_vert[j];

	// Create POV-Ray vector strings for each of the vertices
	for (k = 0; k<n; k++) {
		l = 3 * f_vert[j + k + 1];
		sprintf(s[k], "<%g,%g,%g>", v[l], v[l + 1], v[l + 2]);
	}

	// Draw the interior of the polygon
	fputs("union{\n", fp);
	for (k = 2; k<n; k++) fprintf(fp, "\ttriangle{%s,%s,%s}\n", s[0], s[k - 1], s[k]);
	fputs("\ttexture{t1}\n}\n", fp);

	// Draw the outline of the polygon
	fputs("union{\n", fp);
	for (k = 0; k<n; k++) {
		l = (k + 1) % n;
		fprintf(fp, "\tcylinder{%s,%s,r}\n\tsphere{%s,r}\n",
			s[k], s[l], s[l]);
	}
	fputs("\ttexture{t2}\n}\n", fp);
}

int main(int argc, char *argv[])
{
    if(argc<=1){
        cout << "No cell data specified." << endl;
        exit(1);
    }
    ifstream ifs(argv[1]);
    if(!ifs){
        cout << "Opening cell data failed. Filename:"<<argv[1] << endl;
        exit(1);
    }

    
    string str;
    int count=0;
	double zmax = 0,zmin=LZ;
    while(getline(ifs,str)){
        int index=0;
        int state=0;
        double x,y,z;
        sscanf(str.data(),"%d %d %*f %*f %*f %*f %lf %lf %lf %*f %*d %*f %*f %*d %*f %*d %*d",&index,&state,&x,&y,&z);
        if(count_cond(state)){
			if (zmax < z)zmax = z;
			if (zmin > z)zmin = z;
        }
    }
	ifs.clear();
	ifs.seekg(0, ios::beg);
	voro::container con(0, LX, 0, LY, zmin, zmax, comp_nx, comp_ny, comp_nz, true, true, false, 8);

	while (getline(ifs, str)) {
		int index = 0;
		int state = 0;
		double x, y, z;
		sscanf(str.data(), "%d %d %*f %*f %*f %*f %lf %lf %lf %*f %*d %*f %*f %*d %*f %*d %*d", &index, &state, &x, &y, &z);
		if (count_cond(state)) {//state==MUSUME||state==ALIVE||state==DEAD||
			con.put(count++, x, y, z);
		}
	}
	voro::c_loop_all cl(con);
	vector<int> neigh, f_vert;
	vector<double> v;
	voro::voronoicell_neighbor c;
	FILE* fp4;
	fp4 = fopen("outp.pov", "w");
	if (cl.start()) do if (con.compute_cell(c, cl)) {
		double x, y, z; int id;
		cl.pos(x, y, z); id = cl.pid();

		// Gather information about the computed Voronoi cell
		c.neighbors(neigh);
		c.face_vertices(f_vert);
		c.vertices(x, y, z, v);

		// Loop over all faces of the Voronoi cell
		for (int i = 0, j = 0; i<neigh.size(); i++) {

			// Draw all quadrilaterals, pentagons, and hexagons.
			// Skip if the neighbor information is smaller than
			// this particle's ID, to avoid double counting. This
			// also removes faces that touch the walls, since the
			// neighbor information is set to negative numbers for
			// these cases.
			//if (neigh[i]>id) {
				draw_polygon(fp4, f_vert, v, j);
			//}

			// Skip to the next entry in the face vertex list
			j += f_vert[j] + 1;
		}
	} while (cl.inc());



	
    //con.draw_cells_pov("cell_voro_test.pov");
	//con.draw_particles_pov("cell_voro_test_p.pov");
    //cout << "Hello World!" << endl;
    return 0;
}
