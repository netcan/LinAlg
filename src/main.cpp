/****************************************************************************
 > File Name: main.cpp
 > Author: Netcan
 > Blog: http://www.netcan666.com/
 > Mail: netcan1996@gmail.com
 > Created Time: 2018-05-25 -- 17:17
 ****************************************************************************/

#include "LinAlg.h"
using namespace LinAlg;

int main() {
	Matrix m({
		 {1,2,3,4},
		 {5,6,7,8},
		 {9,10,11,12},
	});
	Matrix n({
		{1,2,3},
		{4,5,6},
		{7,8,10},
	});
	Matrix r({
		{43, 63, 57, 35},
		{63, 26, 32, 35},
		{57, 32, 78, 76},
		{35, 35, 76, 25},
	});

	r.show();
	printf("det(r) = %lf\n", r.det());

	Vector y({16,22,32,44});
	Vector x = r.LUsolve(y);
	x.show();
	puts("");

	n.inv().show();

	double det = 1;
	for(auto x: r.Jacobi()) {
		printf("%lf ", x);
		det *= x;
	}
	puts("");
	printf("%f\n", det);
	printf("%f\n", r.det());

	return 0;
}
