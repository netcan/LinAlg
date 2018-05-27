/****************************************************************************
 > File Name: main.cpp
 > Author: Netcan
 > Blog: http://www.netcan666.com/
 > Mail: netcan1996@gmail.com
 > Created Time: 2018-05-25 -- 17:17
 ****************************************************************************/

#include "Matrix.h"
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
		{7,8,9},
	});
	Matrix r({
		{43, 63, 57, 35},
		{96, 26, 32, 35},
		{80, 29, 78, 76},
		{5, 97, 10, 25},
	});

	r.show();
	printf("det(r) = %lf\n", r.det());

	Vector y{16,22,32,44};
	Vector x = r.LUsolve(y);
	x.show();

	puts("");

	return 0;
}
