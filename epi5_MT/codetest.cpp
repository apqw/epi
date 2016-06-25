#include "codetest.h"
#include "atomics.h"
#include "cell.h"
#include "cellmanager.h"
#include <cstdio>
#include <cassert>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <memory>
#include "utils.h"
bool atomic_double_test(){
    double dval=0.1;
atomic_double a(dval);
atomic_double b(a); //atomic_double c(atomic_double(dval));
assert(dval==a&&dval==b&&(double)a==(double)b);
printf("define test:a %lf == %lf b %lf == %lf\n",(double)a,dval,(double)b,dval);
a=0.2;
assert(0.2==a);
printf("subst test ok\n");
a=dval;
a+=0.1;
assert(dval+0.1==a);
printf("increment test ok\n");
a=dval;
a-=0.1;
assert(dval-0.1==a);
printf("decrement test ok\n");
a=dval;
a*=0.1;
assert(dval*0.1==a);
printf("mul subst test ok\n");
a=dval;
a/=0.1;
assert(dval/0.1==a);
printf("div subst test ok\n");

a=dval;
double tmp=dval;
double ill_par=dval;
constexpr size_t ITR = 1000000;
tbb::parallel_for(tbb::blocked_range<size_t>(0,ITR),[&](const tbb::blocked_range<size_t>& range){
    for(size_t i=range.begin();i!=range.end();++i){
        a+=0.1;
        ill_par+=0.1;
    }
});
for(size_t i=0;i<ITR;i++){
    tmp+=0.1;
}
assert(tmp==a);
printf("parallel increment test ok, atomic_double:%lf non_atomic_double:%lf non_parallel:%lf\n",(double)a,ill_par,tmp);
return true;
}


bool lfpstack_test(){
Lockfree_push_stack<atomic_double,100> a;
printf("lfps init test ok\n");
a.push_back(0.1);
a.push_back(0.1);
a.push_back(0.1);
a.push_back(0.1);
printf("lfps size:%zd pushed:%lf ok\n",a.size(),(double)a[0]);

a.clear();

printf("lfps size:%zd pushed:%lf ok\n",a.size(),(double)a[0]);

return true;

}

bool cell_test() {
    auto cells = make_unique_c11<CellManager>();
    auto p = cells ->create(ALIVE,0);
    auto q = cells->create(DEAD,0); //auto r = cells->create(DEAD);
	p->x += 0.1;
    printf("p idx:%zd q idx:%zd p prex:%lf\n", p->get_index(), q->get_index(),p->x());
	pos_copy(*cells);
	printf("pos swapped p pos:%lf\n", p->x());
    printf("cell man size:%zd\n", cells->size());
	cells ->add_remove_queue(1);
	cells ->remove_exec();
    printf("cell man size:%zd\n", cells->size());
	cells->remove_exec();
    printf("cell man size:%zd\n", cells->size());
	return true;
}
