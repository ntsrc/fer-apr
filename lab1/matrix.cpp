#include "matrix.h"
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <functional>

namespace apr
{

matrix::matrix(const char *filename)
{
	std::ifstream ifs(filename);
	if (!ifs)
		throw std::runtime_error("file error");

	std::string line;
	ifs = std::getline(ifs, line);
	std::istringstream iss(line);

	std::istream_iterator<value_type> beg(iss), end;
	elem_.insert(elem_.end(), beg, end);

	cols_ = elem_.size();

	beg = std::istream_iterator<value_type>(ifs);
	elem_.insert(elem_.end(), beg, end);

	if (elem_.size() % cols_ != 0)
		throw std::runtime_error("file error");

	rows_ = elem_.size() / cols_;
}


auto matrix::at(size_type row, size_type col) const
{
	if (row < 0 || row >= rows_ || col < 0 || col >= cols_)
		throw std::runtime_error("invalid element");

	return elem[idx(row, col)];
}

auto &matrix::at(size_type row, size_type col)
{
	if (row < 0 || row >= rows_ || col < 0 || col >= cols_)
		throw std::runtime_error("invalid element");

	return elem[idx(row, col)];
}

matrix &matrix::operator+=(const matrix &m)
{
	if (rows_ != m.rows_ || cols_ != m.cols_)
		throw std::runtime_error("incompatible matrices");

	std::transform(elem_.begin(), elem_.end(), m.elem_.begin(), elem_.begin(), std::plus<value_type>());

	return *this;
}

matrix &matrix::operator-=(const matrix &m)
{
	if (rows_ != m.rows_ || cols_ != m.cols_)
		throw std::runtime_error("incompatible matrices");

	std::transform(elem_.begin(), elem_.end(), m.elem_.begin(), elem_.begin(), std::minus<value_type>());

	return *this;
}

matrix &matrix::operator*=(const matrix &m)
{
	if (cols_ != m.rows_)
		throw std::runtime_error("incompatible matrices");

	std::vector<value_type> new_elem(rows_ * m.cols_);

	for (auto i = 0; i != rows_; ++i)
	{
		for (auto j = 0; j != m.cols_; ++j)
		{
			for (auto k = 0; k != cols_; ++k)
				new_elem[idx(i, j)] += elem_[idx(i, k)] * m.elem_[idx(k, j)];
		}
	}

	cols_ = m.cols_;
	elem_ = new_elem;

	return *this;
}

auto &matrix::operator*=(value_type v)
{
	std::transform(elem_.begin(), elem_.end(), elem_.begin(), [s](auto e){ return e * s; })

	return *this;
}

void matrix::LU_decomposition()
{
	if (rows_ != cols_)
		throw std::runtime_error("non-square matrix");

	for (auto i = 0; i != rows_ - 1; ++i)
	{
		if (std::abs(elem_[idx(i, i)]) < eps)
			throw std::runtime_error("LU error");

		for (auto j = i + 1; j != rows_; ++j)
		{
			elem_[idx(j, i)] /= elem_[idx(i, i)];
			for (auto k = i + 1; k != rows_; ++k)
				elem_[idx(j, k)] -= elem_[idx(j, i)] * elem_[idx(i, k)];
		}
	}
}

Matrix Matrix::LUP_decomposition()
{
	if (row != col)
		throw std::runtime_error("non-square matrix");

	std::vector<unsigned> Pv(row);
	for (int i = 0; i < row; ++i)
		Pv[i] = i;

	for (int k = 0; k < row-1; ++k) {
		double pivot = 0;
		int l = k;

		for (int i = k; i < row; ++i) {
			if (std::abs(elem[i*col+k]) > pivot) {
				pivot = std::abs(elem[i*col+k]);
				l = i;
			}
		}

		if (pivot < EPSILON)
			throw std::runtime_error("singular matrix");

		std::swap(Pv[k], Pv[l]);

		for (int j = 0; j < row; ++j)
			std::swap(elem[k*col+j], elem[l*col+j]);

		for (int i = k+1; i < row; ++i) {
			elem[i*col+k] /= elem[k*col+k];
			for (int j = k+1; j < row; ++j)
				elem[i*col+j] -= elem[i*col+k]*elem[k*col+j];
		}
	}

	Matrix P(row);
	for (int i = 0; i < row; ++i) {
		P(i, Pv[i]) = 1;
	}

	return P;
}

void Matrix::print_to_file(const std::string& filename) const
{
	std::ofstream ofs(filename);
	if (!ofs)
		throw std::runtime_error("file error");

	ofs << *this;
}

bool operator==(const Matrix& m1, const Matrix& m2)
{
	if (&m1 == &m2)
		return true;

	if (m1.rows() != m2.rows() || m1.columns() != m2.columns())
		return false;

	for (int i = 0; i < m1.rows(); ++i)
		for (int j = 0; j < m1.columns(); ++j)
			if (std::abs(m1(i,j) - m2(i,j)) > EPSILON)
				return false;

	return true;
}

Matrix operator~(const Matrix& m)
{
	Matrix trans(m.columns(), m.rows());

	for (int i = 0; i < trans.rows(); ++i)
		for (int j = 0; j < trans.columns(); ++j)
			trans(i,j) = m(j,i);

	return trans;
}

Matrix operator!(const Matrix& m)
{
	if (m.rows() != m.columns())
		throw std::runtime_error("singular matrix");

	Matrix inverse(m.rows());

	std::pair<Matrix, Matrix> LUP_P = LUP_decomposition(m);

	for (int i = 0; i < m.columns(); ++i) {
		Matrix b(m.rows(), 1);
		
		for (int j = 0; j < b.rows(); ++j)
			b(j,0) = LUP_P.second(j,i);

		Matrix y = forward_supstitution(LUP_P.first, b);
		Matrix x = backward_supstitution(LUP_P.first, y);

		for (int j = 0; j < inverse.rows(); ++j) {
			inverse(j,i) = x(j,0);
		}
	}

	return inverse;
}

double determinant(const Matrix& m)
{
	if (m.rows() != m.columns())
		throw std::runtime_error("non-square matrix");

	std::pair<Matrix,Matrix> LUP_P = LUP_decomposition(m);

	int n = 0;
	for (int i = 0; i < LUP_P.second.rows(); ++i)
		if (std::abs(LUP_P.second(i,i)) > EPSILON)
			++n;

	int det_P = n % 2 == 0 ? 1 : -1;

	double det_U = LUP_P.first(0,0);
	for (int i = 1; i < LUP_P.first.rows(); ++i)
		det_U *= LUP_P.first(i,i);

	return det_P*det_U;
}

Matrix operator+(const Matrix& m1, const Matrix& m2)
{
	Matrix sum = m1;

	return sum += m2;
}

Matrix operator-(const Matrix& m1, const Matrix& m2)
{
	Matrix dif = m1;

	return dif -= m2;
}

Matrix operator*(const Matrix& m1, const Matrix& m2)
{
	Matrix prod = m1;

	return prod *= m2;
}

Matrix operator*(double s, const Matrix& m)
{
	Matrix scaled = m;

	return scaled *= s;
}

Matrix operator*(const Matrix& m, double s)
{
	return s*m;
}

Matrix operator-(const Matrix& m)
{
	return -1*m;
}

Matrix LU_decomposition(const Matrix& m)
{
	Matrix LU = m;

	LU.LU_decomposition();

	return LU;
}

std::pair<Matrix,Matrix> LUP_decomposition(const Matrix& m)
{
	Matrix LUP = m;

	Matrix P = LUP.LUP_decomposition();

	return std::make_pair(LUP, P);
}

Matrix forward_supstitution(const Matrix& L, const Matrix& b)
{
	if (L.rows() != L.columns() || L.columns() != b.rows() || b.columns() != 1)
		throw std::runtime_error("forward supstitution error");

	Matrix y = b;

	for (int i = 0; i < y.rows()-1; ++i)
		for (int j = i+1; j < y.rows(); ++j)
			y(j,0) -= L(j,i)*y(i,0);

	return y;
}

Matrix backward_supstitution(const Matrix& U, const Matrix& y)
{
	if (U.rows() != U.columns() || U.columns() != y.rows() || y.columns() != 1)
		throw std::runtime_error("backward supstitution error");

	Matrix x = y;

	for (int i = x.rows()-1; i >= 0; --i) {
		if (std::abs(U(i,i)) < EPSILON)
			throw std::runtime_error("backward supstitution error");

		x(i,0) /= U(i,i);
		for (int j = 0; j < i; ++j)
			x(j,0) -= U(j,i)*x(i,0);
	}

	return x;
}

std::ostream& operator<<(std::ostream& os, const Matrix& m)
{
	os << '\n';

	for (int i = 0; i < m.rows(); ++i) {
		os << m(i,0);
		for (int j = 1; j < m.columns(); ++j)
			os << ' ' << m(i,j);
		os << '\n';
	}

	return os;
}

}
