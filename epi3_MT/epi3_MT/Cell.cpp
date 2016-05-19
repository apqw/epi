#include "Cell.h"
#include "util_func.h"
#include "primitive_func.h"

std::atomic<uint32_t> Cell::construction_count = 0;

void Cell::do_interact_w_connected() {

}

double Cell::c_der_to_der::operator()(Cell* me, Cell* oppo) {
	if (is_near_delta(me, oppo, cont::delta_R)) {
		return ljmain_der_near(me, oppo);
	}
	else if (is_near(me, oppo)) {
		return cont::para_ljp2;
	}
	else {
		return ljmain_der_far(me, oppo);
	}
}


double Cell::c_al_air_de_to_al_air_de_fix_mu::operator()(Cell* me, Cell* oppo) {
	if (is_near(me, oppo)) {
		return ljmain(me, oppo);
	}
	else {
		double spf = cont::K_DESMOSOME;
		if (me->agek() > cont::THRESH_SP && oppo->agek() > cont::THRESH_SP) {
			spf = cont::K_TOTAL;
		}
		return adhesion(me, oppo,spf);
	}
}

double Cell::c_fix_mu_to_fix_mu::operator()(Cell* me, Cell* oppo) {
    if (me->pair.get() == oppo) {
		return 0;
	}
	else {
		if (is_near(me, oppo)) {
			double tmp = ljmain(me, oppo);
			if (fabs(tmp) > 100) {
				printf("warn: c_fix_mu_to_fix_mu ljm too strong: ljm=%lf\n", tmp);
			}
			return tmp;
		}
		else {
			double tmp = adhesion(me, oppo, cont::K_DESMOSOME);
			if (fabs(tmp) > 100) {
				printf("warn: c_fix_mu_to_fix_mu adhe too strong: adhe=%lf\n", tmp);
			}
			return tmp;
		}
	}
}

double Cell::c_mu_to_memb::operator()(Cell* me, Cell* oppo) {
	if (is_near(me, oppo)) {
		return ljmain(me, oppo);
	}
	else {
		double spf = 0;
		if (paired_with_fix(me)) {
			spf= cont::Kspring;
		}
		else if (me->rest_div_times > 0) {
			spf = cont::Kspring_d;
		}
		return adhesion(me,oppo, spf);
	}
}

double Cell::c_fix_to_memb::operator()(Cell* me, Cell* oppo) {
	if (is_near(me, oppo)) {
		return ljmain(me, oppo);
	}
	else {
		double distlj = sqrt(cellDistSq(me,oppo));
		double LJ2 = distlj / (me->radius() + oppo->radius());
		return  -(cont::Kspring / distlj) * (LJ2 - 1.0);
	}
}

double Cell::c_memb_to_memb::operator()(Cell* me, Cell* oppo) {
	using namespace cont;
	double rad_sum = me->radius() + oppo->radius();
	double cr_dist = rad_sum*P_MEMB;
	double lambda_dist = (1.0 + P_MEMB)*rad_sum;
	double distSq = cellDistSq(me, oppo);
	if ( distSq< cr_dist*cr_dist) {
		//assert(fabs(ljmain(me, oppo)) < 1000);
		double LJ6 = cr_dist*cr_dist / distSq;
		LJ6 = LJ6*LJ6*LJ6;
		return 4.0*cont::eps_m*LJ6*(LJ6 - 1) / distSq;
	}
	else if (distSq < rad_sum*rad_sum) {
		double distlj = sqrt(distSq);
		double tmp = -(DER_DER_CONST / distlj) * (distlj / cr_dist - 1.0);
		assert(fabs(tmp) < 1000);
		return tmp;

		
	}
	else {
		double distlj = sqrt(distSq);
		double LJ6 = cr_dist / (lambda_dist - distlj);
		LJ6 = LJ6*LJ6;
		LJ6 = LJ6*LJ6*LJ6;
		double tmp= -(DER_DER_CONST / rad_sum)*((1.0 - P_MEMB) / P_MEMB)
			- 4 * eps_m*(LJ6*(LJ6 - 1.0)) / ((lambda_dist - distlj)*distlj);

		assert(fabs(tmp) < 1000);
		return tmp;
	}
}

double Cell::c_other::operator()(Cell* me, Cell* oppo) {
	return ljmain(me, oppo);
}

void Cell::wall_interact() {
	double distlj = 2.0*pos[2]();
	double LJ6 = radius() / pos[2]();
	LJ6 = LJ6*LJ6;
	LJ6 = LJ6*LJ6*LJ6;
	double ljm = 4.0*cont::eps_m*LJ6*(LJ6 - 1.0) / (distlj*distlj);
	pos[2] += cont::DT_Cell* ljm*2.0*pos[2]();
}

void Cell::DER_interact() {
	assert(state() == DER);
	connected_cell.foreach([&](Cell* conn) {
		switch (conn->state())
		{
		case DER:
			if (this>conn)_interact<c_der_to_der>(this, conn);
			break;
		default:
			_interact<c_other>(this, conn);
			break;
		}
	});
}

void Cell::AL_AIR_DE_interact() {
	assert(get_state_mask(state())&(ALIVE_M | AIR_M | DEAD_M));
	connected_cell.foreach([&](Cell* conn) {
		switch (conn->state())
		{
		case ALIVE:case AIR:case DEAD:
			if (this>conn)_interact<c_al_air_de_to_al_air_de_fix_mu>(this, conn);
			break;
		case DER:
			break;
		default:
			_interact<c_al_air_de_to_al_air_de_fix_mu>(this, conn);
			break;
		}
	});
}

void Cell::FIX_interact() {
	assert(state() == FIX);
	set_dermis();
	connected_cell.foreach([&](Cell* conn) {
		
		switch (conn->state())
		{
		case FIX:
			if (this>conn)_interact<c_fix_mu_to_fix_mu>(this, conn);
			break;
		case MUSUME:
			_interact<c_fix_mu_to_fix_mu>(this, conn);
			break;
		case MEMB:
			
			if (dermis == conn) {
				_interact<c_fix_to_memb>(this, conn);
			}
			else {
				_interact<c_other>(this, conn);
			}
			break;
		default:
			//c_other()(this, conn);
			break;
		}
	});
}

void Cell::MUSUME_interact() {
	assert(state() == MUSUME);
	set_dermis();
	connected_cell.foreach([&](Cell* conn) {
		switch (conn->state())
		{
		case MUSUME:
			if (this>conn)_interact<c_fix_mu_to_fix_mu>(this, conn);
			break;
		case MEMB:
			
			if (dermis == conn) {
				_interact<c_mu_to_memb>(this, conn); //fix_to_memb‚É‚È‚Á‚Ä‚½
			}
			else {
				_interact<c_other>(this, conn);
			}
			break;
		default:
			break;
		}
	});
}

void Cell::MEMB_interact() {
	assert(state() == MEMB);
	connected_cell.foreach([&](Cell* conn) {
		switch (conn->state())
		{
		case MEMB:
			
			if(this>conn)_interact<c_memb_to_memb>(this, conn);
			break;
		default:
			break;
		}
	});
}

void Cell::memb_bend_calc1() {
	mbd.dn = sqrt(cellDistSq(this, mbd.memb_r.get()));
	assert(mbd.dn >= 1.0e-5);
	mbd.nv = p_diff_v3(mbd.memb_r->pos, pos) / mbd.dn;

	mbd.dm = sqrt(cellDistSq(this,mbd.memb_u.get()));
	assert(mbd.dm >= 1.0e-5);
	mbd.mv = p_diff_v3(mbd.memb_u->pos, pos) / mbd.dm;
}

void Cell::memb_bend_calc2() {
	mbd.ipn = vec_sum(mbd.nv*mbd.memb_r->mbd.nv);
	mbd.ipm = vec_sum(mbd.mv*mbd.memb_u->mbd.mv);
}

void Cell::update() {
	state.update();
	pos[0].update();
	pos[1].update();
	pos[2].update();
	radius.update();
	//ca2p.update();
	//ca2p_avg.update();
	//IP3.update();
	//ex_inert.update();

	agek.update();
	ageb.update();
	ex_fat.update();
	in_fat.update();
	spring_nat_len.update();
	//div_age_thresh.update();
	rest_div_times.update();
}

Vec3<double> MEMB_bend_data::memb_bend_force_sqr() {
	auto& nvr = memb_r->mbd.nv;
	auto& nvl = memb_l->mbd.nv;
	auto& nvll = memb_ll->mbd.nv;

	auto& mvu = memb_u->mbd.mv;
	auto& mvb = memb_b->mbd.mv;
	auto& mvbb = memb_bb->mbd.mv;

	auto& ipnl = memb_l->mbd.ipn;
	auto& ipnll = memb_ll->mbd.ipn;

	auto& ipmb = memb_b->mbd.ipm;
	auto& ipmbb = memb_bb->mbd.ipn;

	auto& dnl = memb_l->mbd.dn;
	auto& dmb = memb_b->mbd.dm;
	
	double x = cont::KBEND*(
		-(1.0 - ipn)*(nvr[0] - ipn*nv[0]) / dn
		+ (1.0 - ipnl)*((nv[0] - ipnl*nvl[0]) / dnl - (nvl[0] - ipnl*nv[0]) / dn)
		+ (1.0 - ipnll)*(nvll[0] - ipnll*nvl[0]) / dnl

		- (1.0 - ipm)*(mvu[0] - ipm*mv[0]) / dm
		+ (1.0 - ipmb)*((mv[0] - ipmb*mvb[0]) / dmb - (mvb[0] - ipmb*mv[0]) / dm)
		+ (1.0 - ipmbb)*(mvbb[0] - ipmbb*mvb[0]) / dmb);

	double y = cont::KBEND*(
		-(1.0 - ipn)*(nvr[1] - ipn*nv[1]) / dn
		+ (1.0 - ipnl)*((nv[1] - ipnl*nvl[1]) / dnl - (nvl[1] - ipnl*nv[1]) / dn)
		+ (1.0 - ipnll)*(nvll[1] - ipnll*nvl[1]) / dnl

		- (1.0 - ipm)*(mvu[1] - ipm*mv[1]) / dm
		+ (1.0 - ipmb)*((mv[1] - ipmb*mvb[1]) / dmb - (mvb[1] - ipmb*mv[1]) / dm)
		+ (1.0 - ipmbb)*(mvbb[1] - ipmbb*mvb[1]) / dmb);

	double z = cont::KBEND*(
		-(1.0 - ipn)*(nvr[2] - ipn*nv[2]) / dn
		+ (1.0 - ipnl)*((nv[2] - ipnl*nvl[2]) / dnl - (nvl[2] - ipnl*nv[2]) / dn)
		+ (1.0 - ipnll)*(nvll[2] - ipnll*nvl[2]) / dnl

		- (1.0 - ipm)*(mvu[2] - ipm*mv[2]) / dm
		+ (1.0 - ipmb)*((mv[2] - ipmb*mvb[2]) / dmb - (mvb[2] - ipmb*mv[2]) / dm)
		+ (1.0 - ipmbb)*(mvbb[2] - ipmbb*mvb[2]) / dmb);

	return Vec3<double>({ x,y,z });
	
	/*
	return cont::KBEND*(
		-(1.0 - ipn)*(nvr - ipn*nv) / dn
		+ (1.0 - ipnl)*( (nv - ipnl*nvl) / dnl - (nvl - ipnl*nv) / dn)
		+ (1.0 - ipnll)*(nvll - ipnll*nvl) / dnl

		- (1.0 - ipm)*(mvu - ipm*mv) / dm
		+ (1.0 - ipmb)*((mv - ipmb*mvb) / dmb - (mvb - ipmb*mv) / dm)
		+ (1.0 - ipmbb)*(mvbb - ipmbb*mvb) / dmb);
	*/
}
/*
	do memb_bend_calc1() and memb_bend_calc2() before doing this
*/
void Cell::memb_bend_interact() {
	pos += cont::DT_Cell*mbd.memb_bend_force_sqr();
}

void Cell::pair_interact() {
	if (pair!=nullptr&&this > pair.get()) {
		double dist = sqrt(cellDistSq(this, pair.get()));
		double force = cont::Kspring_division*(dist - spring_nat_len());
        auto dum = (cont::DT_Cell*force)*p_diff_v3(pos, pair->pos)/dist;
		/*
        printf("pair interact\n");
		printf("addr:%d\n", this);
        printf("cx:%lf,cy:%lf,cz:%lf\n",pos[0](),pos[1](),pos[2]());
        printf("diffx:%lf,diffy:%lf,diffz:%lf\n",dum[0](),dum[1](),dum[2]());
		*/
		//•„†‹t‚É‚È‚Á‚Ä‚½
		//‚Î‚Ë‚Ì—Í‚ÌŒü‚«l‚¦‚ë
        pos -= dum;
		pair->pos += dum;
	}
}

void Cell::set_dermis() {
	double d1Sq = cont::LX*cont::LX;
	double distanceSq=0;
	dermis = nullptr;
	connected_cell.foreach([&](Cell* cptr) {
		if (cptr->state() == MEMB) {
			distanceSq = cellDistSq(this, cptr);
			if (distanceSq < d1Sq) {//2æ‚É‚È‚Á‚Ä‚È‚©‚Á‚½
				d1Sq = distanceSq;
				dermis = cptr;
			}
		}
	});
}

bool Cell::divide_try()
{
	using namespace cont;
	if (pair != nullptr)return false;

    double div_gamma = DT_Cell*(is_malignant ? accel_div:1)*eps_kb*u0 / (div_age_thresh*stoch_div_time_ratio);
	if (STOCHASTIC&&genrand_real()>div_gamma) {
		return false;
	}
	pair = std::make_shared<Cell>(
		MUSUME,
		pos,
		radius,
		ca2p,
        ca2p,//this is avg value,do not use orignal avg
		IP3,
		ex_inert,
		0,0,//set ages 0
		0,0,//set fats 0
		delta_L,//init
		agki_max, //musume thresh
		rest_div_times,
		is_malignant,
		is_touch
		);
	//do pairing pair->pair
	spring_nat_len.force_set_next_value(delta_L);
	ageb.force_set_next_value(0);
	if (state() == MUSUME) {
		pair->rest_div_times -= 1;
		rest_div_times -= 1;
	}
	//FIX divs infinitely
	assert(dermis != nullptr);
	auto divdir = div_direction(this, dermis);
    pos += divdir*(0.5*delta_L);
    pair->pos-= divdir*(0.5*delta_L);



	return true;
}
double Cell::ageb_const() {
	using namespace cont;
	return (is_malignant ? accel_div : 1)*eps_kb*(u0 + alpha_b*min0(ca2p_avg() - u0));
}

double Cell::agek_const() {
	using namespace cont;
	assert(state() == DEAD || state() == AIR || state() == ALIVE);
	if (state() == DEAD || state() == AIR) {
		return eps_kk*u0;
	}
	else {
		return (is_malignant ? accel_diff : 1.0)*eps_ks*
			(S0 + alpha_k*min0(ca2p_avg() - u0));
	}
}
void Cell::FIX_state_renew() {
	assert(state() == FIX);
	set_dermis();
	if (pos[2]() > 10) {
		//printf("warn: too high z value of FIX. z=%lf\n",pos[2]());
	}
    if(dermis==nullptr){
        printf("err\n");
        printf("x:%lf,y:%lf,z:%lf\n",pos[0](),pos[1](),pos[2]());
        printf("connected_num:%d\n",connected_cell.count());
        assert(dermis != nullptr);
    }


	//if (state() == MUSUME&&rest_div_times <= 0)return false;
	if (ageb() >= div_age_thresh*(1.0 - cont::stoch_div_time_ratio)) {
		pair_generated=divide_try();
	}
	else {
        ageb += cont::DT_Cell*ageb_const(); //ok
	}

	
}
void Cell::set_as_deleted() {
	pending_kill = true;
}

bool Cell::should_deleted() {
	return pending_kill;
}

void Cell::set_as_pair_gen() {
	pair_generated = true;
}

bool Cell::has_new_pair() {
	return pair_generated;
}

void Cell::set_as_no_more_new_pair()
{
	pair_generated = false;
}

void Cell::MUSUME_state_renew() {
	assert(state() == MUSUME);
	set_dermis();
#ifdef  UNPAIR_DBG
	dermis = nullptr;
#endif //  UNPAIR_DBG

	
	if (dermis == nullptr && pair == nullptr) {
		if (SYSTEM == WHOLE) {
			printf("ALIVE detected\n");
			state = ALIVE;
		}
		else if (SYSTEM == BASAL) {
			state = DISA;
			set_as_deleted();
		}
		return;
	}

	if (rest_div_times>0&&ageb() >= div_age_thresh*(1.0 - cont::stoch_div_time_ratio)) {
		pair_generated = divide_try();
	}
	else {
		ageb += cont::DT_Cell*ageb_const();
	}

}

void Cell::DEAD_AIR_state_renew() {
	using namespace cont;
	assert(state() == DEAD || state() == AIR);
	if (agek() >= ADHE_CONST && connected_cell.count() <= DISA_conn_num_thresh) {
		state = DISA;
		set_as_deleted();
	}
	else {
		agek += DT_Cell*agek_const();
	}
}
double Cell::k_lipid_release() {
	using namespace cont;
	return 0.25*lipid_rel*(1 + tanh((ca2p_avg() - ubar) / 0.01))*(1 + tanh((agek() - THRESH_SP) / delta_lipid));
}

double Cell::k_lipid() {
	using namespace cont;
	return 0.25*lipid*(1 + tanh(( ubar- ca2p_avg()) / 0.01))*(1 + tanh((agek() - THRESH_SP) / delta_lipid));
}
void Cell::ALIVE_state_renew() {
	using namespace cont;
	assert(state() == ALIVE);

	if (agek() >= THRESH_DEAD) {
		state = DEAD;
	}
	else {
		agek += DT_Cell*agek_const();
		double tmp = k_lipid_release()*in_fat();
		in_fat += DT_Cell*(k_lipid()*(1.0 - in_fat()) - tmp);
		ex_fat += DT_Cell*tmp;
	}
}

void Cell::pair_disperse() {
	using namespace cont;
	assert(pair != nullptr);
	assert(pair->pair.get() == this);
	double rad_sum = radius() + pair->radius();
	double unpair_th = unpair_dist_coef*rad_sum;
	double distSq = 0;
	if (spring_nat_len < 2.0*radius()) {
		spring_nat_len += DT_Cell*eps_L;
		pair->spring_nat_len.force_set_next_value(spring_nat_len() + DT_Cell*eps_L);
	}
	else if ((distSq=cellDistSq(this, pair.get()))>unpair_th*unpair_th) {
		spring_nat_len.force_set_next_value(0);
		pair->spring_nat_len.force_set_next_value(0);
		pair->pair = nullptr;
		pair = nullptr;
		printf("unpaired. distSq:%lf\n",distSq);
	}
}

void Cell::set_lattice()
{
	using namespace cont;
	lat[0] = (int)(pos[0]() / dx);
	lat[1] = (int)(pos[1]() / dy);
	lat[2] = (int)(pos[2]() / dz);
	
	if (lat[0] < 0) {
		lat[0] +=NX;
	}
	else if (lat[0] >= NX) {
		lat[0] -=NX;
	}

	if (lat[1] < 0) {
		lat[1] += NY;
	}
	else if (lat[1] >= NY) {
		lat[1] -= NY;
	}

	if (lat[2] < 0) {
		lat[2] += NZ;
	}
	else if (lat[2] >= NZ) {
		lat[2] -= NZ;
	}
}

Cell::~Cell()
{
}
