/****************************************************************************
 > File Name: Matrix.cpp
 > Author: Netcan
 > Blog: http://www.netcan666.com/
 > Mail: netcan1996@gmail.com
 > Created Time: 2018-05-25 -- 17:09
 ****************************************************************************/

#include "Matrix.h"
using namespace Matrix;

void Vector::show() {
	printf(type == VecType::Row ? "R":"C");
	printf("(");
	for(size_t i = 0; i < data.size(); ++i) printf("%lf%c", data[i], i == data.size() - 1 ? ')':' ');
}

Vector Vector::T() {
	Vector ret = *this;
	ret.type = this->type == VecType::Row ?
	           VecType::Col: VecType::Row;
	return ret;
}

