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
	Vector x({1,2,3}, VecType::Row);
	Vector y({4,5,6,7}, VecType::Row);
	y.setType(VecType::Col);

	Matrix m({
		 {1,2,3,4},
		 {5,6,7,8},
		 {9,10,11,12},
	});
	Matrix n({
		 {5,6,7,8},
		 {1,2,3,4},
		 {9,10,11,12},
	 });

	(m * n.T()).show();

	return 0;
}
