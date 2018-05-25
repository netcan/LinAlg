/****************************************************************************
 > File Name: main.cpp
 > Author: Netcan
 > Blog: http://www.netcan666.com/
 > Mail: netcan1996@gmail.com
 > Created Time: 2018-05-25 -- 17:17
 ****************************************************************************/

#include "Matrix.h"
using namespace Matrix;

int main() {
	Vector x({1,2,3}, VecType::Row);
	x.T().show();
	return 0;
}
