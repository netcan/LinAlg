/****************************************************************************
 > File Name: Matrix.h
 > Author: Netcan
 > Blog: http://www.netcan666.com/
 > Mail: netcan1996@gmail.com
 > Created Time: 2018-05-25 -- 17:08
 ****************************************************************************/

#ifndef MATRIX_MATRIX_H
#define MATRIX_MATRIX_H

#include <vector>
#include "stdio.h"
using std::vector;
using std::initializer_list;

namespace Matrix {
	enum class VecType {
		Row, Col
	};
	class Vector {
	private:
		vector<double> data;
		VecType type;
	public:
		Vector(size_t n = 1, VecType type = VecType::Col): data(n, 0), type(type) {};
		Vector(initializer_list<double> data, VecType type = VecType::Col): data(data) {}

		VecType getType() const { return type; }
		void setType(VecType type) { Vector::type = type; }
		void show();
		Vector T();

	};
}

#endif //MATRIX_MATRIX_H
