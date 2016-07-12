#include <iostream>
#include <fstream>
#include <string>
#include <voro/voro++.hh>
#include <vector>
#include <cassert>
#define ASSERT assert
#include "Delaunay.h"
using namespace std;

static constexpr double LX=50.0;
static constexpr double LY=50.0;
static constexpr double LZ=100.0;
static constexpr int NMX=100;
static constexpr int NMY=100;
static constexpr double len_th=8.0*1.4142*2.0/4.0;
const double Phi=0.5*(1+sqrt(5.0));
const double phi=0.5*(1-sqrt(5.0));
static constexpr double ccoef=5.0;

static constexpr double grd_e=1.4/2.0;
static constexpr int GNX=LX/grd_e;
static constexpr int GNY=LY/grd_e;
static constexpr double ADHE_CONST=31.3;
static constexpr double ADHE_delta=1.0;
#define b_eps (-0.1)
#define trunc_eps 0.1

static constexpr int comp_nx=8,comp_ny=8,comp_nz=8;
class wall_initial_shape : public voro::wall {
         public:
                 wall_initial_shape() {
 
                         // Create a dodecahedron
                         v.init(-2,2,-2,2,-2,2);
                         v.plane(0,ccoef*Phi,ccoef);v.plane(0,-ccoef*Phi,ccoef);v.plane(0,ccoef*Phi,-ccoef);
                         v.plane(0,-ccoef*Phi,-ccoef);v.plane(ccoef,0,ccoef*Phi);v.plane(-ccoef,0,ccoef*Phi);
                         v.plane(ccoef,0,-ccoef*Phi);v.plane(-ccoef,0,-ccoef*Phi);v.plane(ccoef*Phi,ccoef,0);
                         v.plane(-ccoef*Phi,ccoef,0);v.plane(ccoef*Phi,-ccoef,0);v.plane(-ccoef*Phi,-ccoef,0);
                 };
                 bool point_inside(double x,double y,double z) {return true;}
  bool cut_cell(voro::voronoicell &c,double x,double y,double z) {

                        // Set the cell to be equal to the dodecahedron
                         c=v;
                         return true;
                 }
  bool cut_cell(voro::voronoicell_neighbor &c,double x,double y,double z) {
 
                         // Set the cell to be equal to the dodecahedron
                         c=v;
                         return true;
                 }
         private:
  voro::voronoicell v;
 };
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
struct pos{
double x,y,z;
bool valid=true;
};
       struct dead_p{
	 double z_val;
	 bool registered=false;
	 void set(double z){
	   if(registered){
	     z_val=z>z_val?z:z_val;
	   }else{
	     z_val=z;registered=true;
	   }
	 }
       };
void draw_polygon(FILE *fp, vector<int> &f_vert, vector<double> &v, int j,int state) {
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
    fprintf(fp,"\ttexture{t%d}\n}\n",state);
    //fputs("\ttexture{t1}\n}\n", fp);

    // Draw the outline of the polygon
    fputs("union{\n", fp);
    for (k = 0; k<n; k++) {
        l = (k + 1) % n;
        fprintf(fp, "\tcylinder{%s,%s,r}\n\tsphere{%s,r}\n",
                s[k], s[l], s[l]);
    }
    fputs("\ttexture{c1}\n}\n", fp);
}
double pos_len(const pos& p1,const pos& p2){
double diffx=p1.x-p2.x;
double diffy=p1.y-p2.y;
double diffz=p1.z-p2.z;
return diffx*diffx+diffy*diffy+diffz*diffz;
}
bool pos_len_valid(const pos& p1,const pos& p2){
return pos_len(p1,p2)<len_th;
}
bool memb_tri_valid(const pos& p1,const pos& p2,const pos& p3){
	return pos_len_valid(p1,p2)&&pos_len_valid(p2,p3)&&pos_len_valid(p3,p1);
}

void right_fix(pos& p1,pos& p1r){
	if(p1.x-p1r.x>LX/2.0)p1r.x+=LX;
}

void bot_fix(pos& p1,pos& p1b){
	if(p1.y-p1b.y>LY/2.0)p1b.y+=LY;
}

void trunc_pos(pos& p){
	if(p.x>LX+trunc_eps)p.x=LX+trunc_eps;
	if(p.x<0-trunc_eps)p.x=-trunc_eps;
	if(p.y>LY+trunc_eps)p.y=LY+trunc_eps;
	if(p.y<0-trunc_eps)p.y=-trunc_eps;
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
    std::vector<int> state_arr;
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
    voro::container con(0, LX, 0, LY, zmin, zmax, comp_nx, comp_ny, comp_nz, false, false, true, 8);
    //wall_initial_shape(wis);
    // con.add_wall(wis);
   FILE* fp4;
    fp4 = fopen("outp.pov", "w");
    //fprintf(fp4,"intersection{box{\n<%lf,%lf,-1000>,<%lf,%lf,%lf> texture{r1}\n}\n",b_eps,b_eps,LX-b_eps,LY-b_eps,LZ+1000);
    //fputs("blob{\n\tthreshold mth\n",fp4);
    bool memb_flg=false;
       std::vector<pos> memb_pos_store;
       std::vector<double> agek_c;

       dead_p dead_grid[GNX][GNX];
       
    while (getline(ifs, str)) {
        int index = 0;
        int state = 0;
        double x, y, z,agek;
        sscanf(str.data(), "%d %d %*f %*f %f %*f %lf %lf %lf %*f %*d %*f %*f %*d %*f %*d %*d", &index, &state,&agek, &x, &y, &z);
	if(memb_flg==false&&state!=MEMB){
	  memb_flg=true;
	//  fprintf(fp4,"\ttexture{r1}\n}texture{r1}}");
	    }
	if(memb_flg==false){
	memb_pos_store.push_back({x,y,z});
	//  fprintf(fp4,"\tsphere{<%lf,%lf,%lf> ,mrad,mstr scale z*mscl translate <0,0,mztr>}\n",x,y,z);
	}
        if (count_cond(state)) {//state==MUSUME||state==ALIVE||state==DEAD||
            state_arr.push_back(state);
	    if(state==DEAD&&agek>=ADHE_CONST-ADHE_delta){
	      int grx=x/grd_e;grx=(grx+GNX)%GNX;
	      int gry=y/grd_e;gry=(gry+GNY)%GNY;
	      dead_grid[grx][gry].set(z);
	    }
	    agek_c.push_back(agek);
            con.put(count++, x, y, z);
        }
    }
	 Delaunay delaunay;
	 vertexSet dvset;
	 for(int gix=-GNX;gix<2*GNX;gix++){
	   for(int giy=-GNY;giy<2*GNY;giy++){
	     dead_p& dp_ref=dead_grid[(gix+GNX)%GNX][(giy+GNY)%GNY];
	     if(dp_ref.registered){
	       vertex tmpv(gix,giy);
	       tmpv.user_data=&dp_ref;
	       dvset.insert(tmpv);
	     }

	   }

	 }
	 triangleSet trset;
	 delaunay.Triangulate(dvset,trset);

printf("memb_pos_size:%d\n",memb_pos_store.size());
assert(memb_pos_store.size()==NMX*NMY);
/*
if(memb_pos_store[0].x>LX/2.0){
memb_pos_store[0].x-=LX;
}
if(memb_pos_store[0].y>LY/2.0){
memb_pos_store[0].y-=LY;
}
*/
pos* memb_d=new pos[NMX*3*NMY*3];
for(int my=0;my<NMY-1;my++){
for(int mx=0;mx<NMX-1;mx++){
pos& mpos=memb_pos_store[my*NMX+mx];
pos& mposr=memb_pos_store[my*NMX+mx+1];
pos& mposb=memb_pos_store[(my+1)*NMX+mx];
pos& mposrb=memb_pos_store[(my+1)*NMX+mx+1];

right_fix(mpos,mposr);
bot_fix(mpos,mposb);
right_fix(mposb,mposrb);
bot_fix(mposr,mposrb);


}}
#define MEMB_COPY(ix,iy) \
for(int my=0;my<NMY;my++){\
for(int mx=0;mx<NMX;mx++){\
memb_d[(my+iy*NMY)*NMX*3+mx+ix*NMX]=memb_pos_store[my*NMX+mx];\
pos& tm=memb_d[(my+iy*NMY)*NMX*3+mx+ix*NMX];\
tm.x+=(ix-1)*LX;\
tm.y+=(iy-1)*LY;\
}}\
0

MEMB_COPY(0,0);
MEMB_COPY(0,1);
MEMB_COPY(0,2);
MEMB_COPY(1,0);
MEMB_COPY(1,1);
MEMB_COPY(1,2);
MEMB_COPY(2,0);
MEMB_COPY(2,1);
MEMB_COPY(2,2);
int mstartx=NMX-8;
int mendx=NMX*2+8;
int mstarty=NMY-8;
int mendy=NMY*2+8;
double upz=5.0;
	fprintf(fp4,"mesh{\n");
for(int my=mstarty;my<=mendy;my++){
for(int mx=mstartx;mx<=mendx;mx++){
pos& mpos=memb_d[my*NMX*3+mx];
pos& mposr=memb_d[my*NMX*3+mx+1];
pos& mposb=memb_d[(my+1)*NMX*3+mx];
pos& mposrb=memb_d[(my+1)*NMX*3+mx+1];
trunc_pos(mpos);
trunc_pos(mposr);
trunc_pos(mposb);
trunc_pos(mposrb);
	fprintf(fp4,"\ttriangle{<%lf,%lf,%lf>,<%lf,%lf,%lf>,<%lf,%lf,%lf>}\n",mpos.x,mpos.y,mpos.z,mposr.x,mposr.y,mposr.z,mposb.x,mposb.y,mposb.z);
fprintf(fp4,"\ttriangle{<%lf,%lf,%lf>,<%lf,%lf,%lf>,<%lf,%lf,%lf>}\n\n",mposr.x,mposr.y,mposr.z,mposb.x,mposb.y,mposb.z,mposrb.x,mposrb.y,mposrb.z);

	fprintf(fp4,"\ttriangle{<%lf,%lf,%lf+mupz>,<%lf,%lf,%lf+mupz>,<%lf,%lf,%lf+mupz>}\n",mpos.x,mpos.y,mpos.z,mposr.x,mposr.y,mposr.z,mposb.x,mposb.y,mposb.z);
fprintf(fp4,"\ttriangle{<%lf,%lf,%lf+mupz>,<%lf,%lf,%lf+mupz>,<%lf,%lf,%lf+mupz>}\n\n",mposr.x,mposr.y,mposr.z,mposb.x,mposb.y,mposb.z,mposrb.x,mposrb.y,mposrb.z);

if(my==mstarty){
	fprintf(fp4,"\ttriangle{<%lf,%lf,%lf+mupz>,<%lf,%lf,%lf>,<%lf,%lf,%lf>}\n",mpos.x,mpos.y,mpos.z,mpos.x,mpos.y,mpos.z,mposr.x,mposr.y,mposr.z);
	fprintf(fp4,"\ttriangle{<%lf,%lf,%lf+mupz>,<%lf,%lf,%lf+mupz>,<%lf,%lf,%lf>}\n",mpos.x,mpos.y,mpos.z,mposr.x,mposr.y,mposr.z,mposr.x,mposr.y,mposr.z);
}

if(my==mendy){
	fprintf(fp4,"\ttriangle{<%lf,%lf,%lf+mupz>,<%lf,%lf,%lf>,<%lf,%lf,%lf>}\n",mposb.x,mposb.y,mposb.z,mposb.x,mposb.y,mposb.z,mposrb.x,mposrb.y,mposrb.z);
	fprintf(fp4,"\ttriangle{<%lf,%lf,%lf+mupz>,<%lf,%lf,%lf+mupz>,<%lf,%lf,%lf>}\n",mposb.x,mposb.y,mposb.z,mposrb.x,mposrb.y,mposrb.z,mposrb.x,mposrb.y,mposrb.z);
}

if(mx==mstartx){
	fprintf(fp4,"\ttriangle{<%lf,%lf,%lf+mupz>,<%lf,%lf,%lf>,<%lf,%lf,%lf>}\n",mpos.x,mpos.y,mpos.z,mpos.x,mpos.y,mpos.z,mposb.x,mposb.y,mposb.z);
	fprintf(fp4,"\ttriangle{<%lf,%lf,%lf+mupz>,<%lf,%lf,%lf+mupz>,<%lf,%lf,%lf>}\n",mpos.x,mpos.y,mpos.z,mposb.x,mposb.y,mposb.z,mposb.x,mposb.y,mposb.z);
}
if(mx==mendx){
	fprintf(fp4,"\ttriangle{<%lf,%lf,%lf+mupz>,<%lf,%lf,%lf>,<%lf,%lf,%lf>}\n",mposr.x,mposr.y,mposr.z,mposr.x,mposr.y,mposr.z,mposrb.x,mposrb.y,mposrb.z);
	fprintf(fp4,"\ttriangle{<%lf,%lf,%lf+mupz>,<%lf,%lf,%lf+mupz>,<%lf,%lf,%lf>}\n",mposr.x,mposr.y,mposr.z,mposrb.x,mposrb.y,mposrb.z,mposrb.x,mposrb.y,mposrb.z);
}
}
}

fprintf(fp4,"translate <0,0,mztr>\ntexture{m1}}\n");
    voro::c_loop_all cl(con);
    vector<int> neigh, f_vert;
    vector<double> v;
    voro::voronoicell_neighbor c;
    //fprintf(fp4,"difference{\n");
    fprintf(fp4,"union{\n");
    vector<pos> cell_pos;
    if (cl.start()) do if (con.compute_cell(c, cl)) {
                double x, y, z; int id;
                cl.pos(x, y, z); id = cl.pid();

		//fprintf(fp4,"\tsphere{<%lf,%lf,%lf>,crad texture{cell1}}\n",x,y,z);
                // Gather information about the computed Voronoi cell
                c.neighbors(neigh);
                c.face_vertices(f_vert);
                c.vertices(x, y, z, v);
		bool added=false;
                // Loop over all faces of the Voronoi cell
                for (int i = 0, j = 0; i<neigh.size(); i++) {
		  /*
		  bool overz=true;
    for (int k = 0; k<f_vert[j]; k++) {
        int l = 3 * f_vert[j + k + 1];
	
	overz=overz&&v[l+2]>=zmax-4.0;
    }
		  */
		  
		  
                    // Draw all quadrilaterals, pentagons, and hexagons.
                    // Skip if the neighbor information is smaller than
                    // this particle's ID, to avoid double counting. This
                    // also removes faces that touch the walls, since the
                    // neighbor information is set to negative numbers for
                    // these cases.
		  // if ((neigh[i]>id||neigh[i]<0)) {
                    draw_polygon(fp4, f_vert, v, j,state_arr[id]);
		    // }
		  if(!added&&neigh[i]<0){
		    cell_pos.push_back({x,y,z});added=true;
		  }

                    // Skip to the next entry in the face vertex list
                    j += f_vert[j] + 1;
                }
		if(!added)fprintf(fp4,"\tsphere{<%lf,%lf,%lf>,crad texture{cell1}}\n",x,y,z);
            } while (cl.inc());
    fprintf(fp4,"}\n");
    /*
fprintf(fp4,"intersection{box{\n<%lf,%lf,-1000>,<%lf,%lf,%lf> texture{r1}\n}\n",b_eps,b_eps,LX-b_eps,LY-b_eps,LZ+1000);
 fprintf(fp4,"union{\n");
 for(pos& p:cell_pos){
   fprintf(fp4,"\tsphere{<%lf,%lf,%lf>,crad texture{cell1}}\n",p.x,p.y,p.z);
 }
 fprintf(fp4,"}}\n");
    */
    /*
fputs("blob{\n\tthreshold mth\n",fp4);
    for(int ii=0;ii<posstore.size();ii++){
      fprintf(fp4,"\tsphere{<%lf,%lf,%lf> ,mrad,mstr scale z*mscl translate <0,0,mztr2>}\n",posstore[ii].x,posstore[ii].y,posstore[ii].z);
	}
fprintf(fp4,"\ttexture{r1}\n}texture{r1}}");
    fprintf(fp4,"}\n");
    */
    fclose(fp4);

    //con.draw_cells_pov("cell_voro_test.pov");
    //con.draw_particles_pov("cell_voro_test_p.pov");
    //cout << "Hello World!" << endl;
delete[] memb_d;
    return 0;
}
