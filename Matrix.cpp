/****************************************************************************
 > File Name: Matrix.cpp
 > Author: Netcan
 > Blog: http://www.netcan666.com/
 > Mail: netcan1996@gmail.com
 > Created Time: 2018-05-25 -- 17:09
 ****************************************************************************/

#include "Matrix.h"
#include <cmath>
#include <algorithm>
using namespace LinAlg;
using std::make_pair;

bool LinAlg::Eq(double a, double b) {
	return fabs(a-b) < eps;
}

bool LinAlg::isZero(double a) {
	return fabs(a) < eps;
}

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
			ret._data.val += lhs[i] * rhs[i];
	} else {
		ret.dtype = RetType::DType::MATRIX;
		ret._data.matVal = new Matrix(lhs.getSize(), rhs.getSize());
		for(size_t i = 0; i < lhs.getSize(); ++i)
			for(size_t j = 0; j < rhs.getSize(); ++j)
				ret._data.matVal->data[i][j] = lhs[i] * rhs[j];
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

Matrix Matrix::LUdecomp() const {
	// A = LU
	Matrix LU(*this); // 最终结果存放到一个矩阵中
	size_t s = std::min(getRowSize(), getColSize());
	for(size_t k = 0; k < s; ++k) {
		double x = 1.0 / LU.data[k][k];
		// l1 = 1/a_{11} * ...
		for(size_t i = k + 1; i < getRowSize(); ++i)
			LU.data[i][k] *= x;
		// A -= lu^t,
		for(size_t i = k + 1; i < getRowSize(); ++i)
			for(size_t j = k + 1; j < getColSize(); ++j)
				LU.data[i][j] -= LU.data[i][k] * LU.data[k][j];
	}
	return LU;
}

void Matrix::PLU_P_update(vector<size_t> &P, const Matrix &mat, size_t k) const {
	// 选主元
	size_t row_idx = k;
	double max_val = fabs(mat.data[P[k]][k]);
	for(size_t i = k; i < mat.getRowSize(); ++i) {
		double x = fabs(mat.data[P[i]][k]);
		if(x > max_val) {
			x = max_val;
			row_idx = i;
		}
	}
	std::swap(P[k], P[row_idx]);
}

pair<vector<size_t>, Matrix> Matrix::PLUdecomp() const {
	// PA = A' = LU
	vector<size_t> P(getRowSize());
	Matrix LU(*this);
	size_t s = std::min(getRowSize(), getColSize());
	for(size_t i = 0; i < getRowSize(); ++i) P[i] = i;

	for(size_t k = 0; k < s; ++k) {
		// 选主元
		PLU_P_update(P, LU, k);
		double x = 1.0 / LU.data[P[k]][k];
		// l1 = 1/a_{11} * ...
		for(size_t i = k + 1; i < getRowSize(); ++i)
			LU.data[P[i]][k] *= x;
		// A -= lu^t,
		for(size_t i = k + 1; i < getRowSize(); ++i)
			for(size_t j = k + 1; j < getColSize(); ++j)
				LU.data[P[i]][j] -= LU.data[P[i]][k] * LU.data[P[k]][j];
	}
	return make_pair(P, LU);
}

double Matrix::det() const {
	// det(A) = det(PLU) = det(P)det(L)det(U) = det(P)det(U)
	assert(getRowSize() == getColSize() && "size of matrix should be square");
	PLUType PLU = PLUdecomp();
	const auto &P = PLU.first;
	const auto &LU = PLU.second;

	bool sign = false; // 负号
	for(size_t i = 0; i < PLU.first.size(); ++i)
		for(size_t j = i + 1; j < PLU.first.size(); ++j)
			if(PLU.first[i] > PLU.first[j]) sign = !sign;
	double ret = sign?-1.0:1.0;
	for(size_t i = 0; i < LU.data.size(); ++i)
		ret *= LU.data[P[i]][i];
	return ret;
}

Vector Matrix::LUsolve_L(const Matrix &LU, const Vector &y) const {
	// 求出Lz = y的z
	Vector z(y.getSize());
	for(size_t i = 0; i < z.getSize(); ++i) {
		z[i] = y[i];
		for (size_t j = 0; j < i; ++j)
			z[i] -= LU.data[i][j] * z[j];
	}
	return z;
}

Vector Matrix::LUsolve_U(const Matrix &LU, const Vector &z) const {
	// 求出Ux = z的x
	Vector x(z.getSize());
	for(int i = z.getSize() - 1; i >=0; --i) {
		x[i] = z[i];
		for(int j = i + 1; j < z.getSize(); ++j)
			x[i] -= LU.data[i][j] * x[j];
		x[i] /= LU.data[i][i];
	}
	return x;
}


Vector Matrix::LUsolve(const Vector &y) const {
	// Ax = y
	assert(y.getSize() == getRowSize());
	assert(! isZero(det()) && "det is equal 0");
	Matrix LU = std::move(LUdecomp());
	Vector x(y.getSize());
	x = std::move(LUsolve_L(LU, y));
	x = std::move(LUsolve_U(LU, x));
	return x;
}

Matrix Matrix::inv() const {
	// A(x_1, ..., x_n) = (e_1, ..., e_n)
	assert(! isZero(det()) && "matrix is not invertible");
	size_t n = getRowSize();
	Matrix LU = std::move(LUdecomp());
	Matrix ret(n, n);
	Vector e(n), x(n);
	e[0] = 1.0;
	for(size_t j = 0; j < n; ++j) {
		if(j >= 1) std::swap(e[j], e[j-1]);
		// 解出Ax = e
		x = std::move(LUsolve_L(LU, e));
		x = std::move(LUsolve_U(LU, x));
		for(size_t i = 0; i < n; ++i)
			ret.data[i][j] = x[i];
	}
	return ret;
}

void Matrix::R(double cosphi, size_t p, size_t q, bool T) {
	// 通过平面旋转进行相似变换
	double sinphi = sqrt(1 - cosphi * cosphi);
	size_t n = getRowSize();

	for(size_t k = 0; k < n; ++k) {
		double  ap =  data[T?p:k][T?k:p] * cosphi + data[T?q:k][T?k:q] * sinphi,
				aq = -data[T?p:k][T?k:p] * sinphi + data[T?q:k][T?k:q] * cosphi;
		data[T?p:k][T?k:p] = ap;
		data[T?q:k][T?k:q] = aq;
	}

	return;
}

vector<double> Matrix::Jacobi() const {
	assert(getRowSize() == getColSize());
	size_t n = getRowSize();
	for(size_t i = 0; i < n; ++i)
		for(size_t j = 0; j < n; ++j)
			assert(Eq(data[i][j], data[j][i]));

	vector<double> eigenvals(n);
	Matrix D(*this); // 对角化

	// g(A) = \sum a[i][i]^2
	auto g = [](const Matrix &A)->double {
		double ret = 0.0;
		for(size_t i = 0; i < A.getRowSize(); ++i)
			ret += A.data[i][i] * A.data[i][i];
		return ret;
	};
	auto f = [](const Matrix &A)->double {
		double ret = 0.0;
		for(size_t i = 0; i < A.getRowSize(); ++i)
			for(size_t j = 0; j < A.getRowSize(); ++j)
				if(i!=j) ret += A.data[i][j] * A.data[i][j];
		return ret;
	};

	double cur_g = g(D), last_g = 0.0;

	while (! Eq(cur_g, last_g)) { // 当g(D)不变时，对角化完成
		for (size_t p = 0; p < n; ++p) { // 逐个处理a_{pq}
			for (size_t q = p + 1; q < n; ++q) {
				if (isZero(D.data[p][q])) continue;
				double  tan2phi = 2 * D.data[p][q] / (D.data[p][p] - D.data[q][q]),
						cosphi = sqrt(0.5 * (1.0 + 1.0/sqrt(1+tan2phi * tan2phi)));
				// D = R^T A R
				D.R(cosphi, p, q, true);
				D.R(cosphi, p, q, false);
			}
		}
		last_g = cur_g;
		cur_g = g(D);
	}

	for(size_t i = 0; i < n; ++i) eigenvals[i] = D.data[i][i];
	return eigenvals;
}
