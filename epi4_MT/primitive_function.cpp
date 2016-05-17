#include "define.h"
#include "primitive_function.h"



/*
p_diff_x:ŽüŠú‹«ŠEðŒ‚É‚¨‚¯‚éxÀ•W‚Ì· x1-x2
*/
double p_diff_x(double x1, double x2) {

	/*
	ŽüŠú‹«ŠEðŒ‚É‚¨‚¯‚éÀ•W‚Ì·‚ÌŽæ‚è•û:
	1.•’Ê‚É·‚ð‚Æ‚é
	2.‚»‚Ì·‚Ìâ‘Î’l‚ªA(·‚ð‚Æ‚Á‚½Ž²‚Ì)ŒvŽZ—Ìˆæ‚Ì’·‚³‚Ì”¼•ª‚æ‚è‘å‚«‚¯‚ê‚ÎŽŸ‚Ì‚æ‚¤‚ÉC³
		2.1. ·‚ª•‰‚È‚çŒvŽZ—Ìˆæ‚Ì’·‚³‚ð‘«‚·
		2.2. ·‚ª³‚È‚çŒvŽZ—Ìˆæ‚Ì’·‚³‚ðˆø‚­

	(ŽüŠú‹«ŠE->’[‚ª‚Â‚È‚ª‚Á‚Ä‚é->‚à‚Á‚Æ‚à—£‚ê‚½‚Æ‚µ‚Ä‚à‚»‚Ì‹——£‚Í—Ìˆæ‚Ì’·‚³‚Ì”¼•ª)
	*/

	using namespace cont;
	double diff = x1 - x2;
	if (diff > 0.5*LX)return diff - LX;
	if (diff <= -0.5*LX)return diff + LX;
	return diff;
}

/*
p_diff_y:ŽüŠú‹«ŠEðŒ‚É‚¨‚¯‚éyÀ•W‚Ì· y1-y2
*/
double p_diff_y(double y1, double y2) {
	using namespace cont;
	double diff = y1 - y2;
	if (diff > 0.5*LY)return diff - LY;
	if (diff <= -0.5*LY)return diff + LY;
	return diff;
}

double p_dist_sq(double x1, double y1, double z1, double x2, double y2, double z2) {
	double diffx = p_diff_x(x1, x2);
	double diffy = p_diff_y(y1, y2);
	double diffz = z1 - z2;
	return diffx*diffx + diffy*diffy + diffz*diffz;
}