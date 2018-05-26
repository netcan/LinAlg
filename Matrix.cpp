/****************************************************************************
 > File Name: Matrix.cpp
 > Author: Netcan
 > Blog: http://www.netcan666.com/
 > Mail: netcan1996@gmail.com
 > Created Time: 2018-05-25 -- 17:09
 ****************************************************************************/

#include "Matrix.h"
using namespace LinAlg;

void Vector::show() {
	printf(type == VecType::Row ? "R":"C");
	printf("(");
	for(size_t i = 0; i < data.size(); ++i) printf("%lf%c", data[i], i == data.size() - 1 ? ')':' ');
	puts("");
	return;
}

Vector Vector::T() {
	Vector ret = *this;
	ret.type = this->type == VecType::Row ?
	           VecType::Col: VecType::Row;
	return ret;
}

Vector& Vector::operator+=(const Vector &rhs) {
	assert(getSize() == rhs.getSize());
	for(int i = 0; i < getSize(); ++i) data[i] += rhs.data[i];
	return *this;
}

Vector& Vector::operator-=(const Vector &rhs) {
	assert(getSize() == rhs.getSize());
	*this += rhs * -1;
	return *this;
}

Vector& Vector::operator*=(double c) {
	for(int i = 0; i < getSize(); ++i) data[i] *= c;
	return *this;
}


RetType LinAlg::operator*(const Vector &lhs, const Vector &rhs) {
	RetType ret;
	if(lhs.type != VecType::Col || rhs.type != VecType::Row) { // 内积
		assert(lhs.getSize() == rhs.getSize());
		ret._data.val = 0;
		for(int i = 0; i < lhs.getSize(); ++i)
			ret._data.val += lhs.data[i] * rhs.data[i];
	} else {
		ret.dtype = RetType::DType::MATRIX;
		ret._data.matVal = new Matrix(lhs.getSize(), rhs.getSize());
		for(size_t i = 0; i < lhs.getSize(); ++i)
			for(size_t j = 0; j < rhs.getSize(); ++j)
				ret._data.matVal->data[i][j] = lhs.data[i] * rhs.data[j];
	}
	return ret;
}



Vector Matrix::getNRowVec(size_t n) const {
	assert(n < getRowSize());
	return Vector(data[n], VecType::Row);
}

Vector Matrix::getNColVec(size_t n) const {
	assert(n < getColSize());
	Vector ret(getRowSize());
	for(int i = getRowSize() - 1; i >= 0; --i)
		ret.data[i] = data[i][n];
	return ret;
}

void Matrix::show() {
	for(size_t i = 0; i < getRowSize(); ++i) {
		for (size_t j = 0; j < getColSize(); ++j)
			printf("%10.6f ", data[i][j]);
		puts("");
	}
	return;
}

Matrix Matrix::T() {
	Matrix ret(getColSize(), getRowSize());
	for(size_t i = 0; i < getRowSize(); ++i)
		for (size_t j = 0; j < getColSize(); ++j)
			ret.data[j][i] = data[i][j];
	return ret;
}

Matrix& Matrix::operator+=(const Matrix &rhs) {
	assert(getRowSize() == rhs.getRowSize() && getColSize() == rhs.getColSize());
	for(size_t i = 0; i < getRowSize(); ++i)
		for (size_t j = 0; j < getColSize(); ++j)
			data[i][j] += rhs.data[i][j];
	return *this;
}

Matrix& Matrix::operator-=(const Matrix &rhs) {
	assert(getRowSize() == rhs.getRowSize() && getColSize() == rhs.getColSize());
	*this += rhs * -1;
	return *this;
}

Matrix& Matrix::operator*=(double c) {
	for(size_t i = 0; i < getRowSize(); ++i)
		for (size_t j = 0; j < getColSize(); ++j)
			data[i][j] *= c;
	return *this;
}



Matrix LinAlg::operator*(Matrix lhs, const Matrix &rhs) {
	assert(lhs.getColSize() == rhs.getRowSize());
	Matrix ret(lhs.getRowSize(), rhs.getColSize());
	for(size_t i = 0; i < lhs.getRowSize(); ++i)
		for (size_t j = 0; j < rhs.getColSize(); ++j)
			ret.data[i][j] = (lhs.getNRowVec(i) * rhs.getNColVec(j))._data.val;
	return ret;
}
